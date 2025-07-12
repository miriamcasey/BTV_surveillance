# packages ------------------------------------------------------------------------------------
library(foreign)
library(janitor)
library(tidyverse)
library(epiR)
library(data.table)
library(parallel)
library(actuar)
library(tidyverse)
library(ggdist)
library(cowplot)
library(ggbeeswarm)
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d3 <-read_csv("d3_data.csv")
wanted_hex_ref <- unique(d3$hex_id)

lookup1 <- read_csv("lookup1.csv")
# functions -----------------------------------------------------------------------------------------

ppv_function <- function(prevalence, sensitivity, specificity) {ppv <- 
  (sensitivity * prevalence) / ( (sensitivity * prevalence) + ((1 - specificity) * (1 - prevalence)) )
return(ppv)
}
npv_function <- function(prevalence, sensitivity, specificity) {npv <-
  (specificity * (1 - prevalence))/(((1 - sensitivity) * prevalence) + ((specificity) * (1 - prevalence)))
return(npv)
}

ppv_function(0.015, 0.892, 0.984)
npv_function(0.015, 0.892, 0.984)

# zero truncated binomial sampling for within-herd infecteds
## if herd infected must have at least one infected animal
## function for zero truncated binomial sampling from "actuar" package
repeat_rbinom <- function(xc, within_herd_prevalence1) 
{
  y <- rztbinom(1, xc, within_herd_prevalence1) 
  return(y)
}

# epiR sample size calculation retruning count herds to test
ss_fun <- function(n_herds, vector_herd_sizes, between_herd_prevalence, within_herd_prevalence, test_sens)
{
  ss <- rsu.sssep.rs2st(H = n_herds, N = vector_herd_sizes, pstar.c = between_herd_prevalence, 
                        pstar.u = within_herd_prevalence, se.p = 0.95, se.c = 0.95, se.u = test_sens)
  
  return(ss)
}

# sample animals within a herd to test
sample_animals_function <- function(herd_size, n_infected_cattle, n_uninfected_cattle, n_cattle_to_sample)
{
  seq1 <- c(rep(0, n_uninfected_cattle), rep(1, n_infected_cattle))
  seq2 <- sample(seq1) # randomise
  samp <- sample(1:herd_size, n_cattle_to_sample, replace = FALSE)
  infected_for_testing <- sum(seq2[samp])
  return(infected_for_testing)
}

# simulate BT surveillance in a hexagon

make_first_table <- function(between_herd_prevalence, 
                                                               within_herd_prevalence,
                                                               hexagon_reference,
                                                               test_sensitivity,
                                                               test_specificity) 
{
  test_data <- d3[d3$hex_id == hexagon_reference ,  ] 
  
  big_herds <- test_data[test_data$hs22 > 100, ]
  
  ### how many herds and animals to sample? ------------------------------------------------------
  
  sampling_advised <- ss_fun(n_herds = nrow(big_herds), vector_herd_sizes = big_herds$hs22,
                             between_herd_prevalence = between_herd_prevalence,
                             within_herd_prevalence = within_herd_prevalence,
                             test_sens = test_sensitivity
  )
  
  herds_to_sample1 <- sampling_advised[[1]]$nsample # count of herds to sample

  
   overall_animal_prevalence <-  within_herd_prevalence * between_herd_prevalence
  ppv1 <- ppv_function(prevalence = overall_animal_prevalence,
                       sensitivity = test_sensitivity,
                       specificity = test_specificity
                       )
  npv1 <- npv_function(prevalence = overall_animal_prevalence,
                       sensitivity = test_sensitivity,
                       specificity = test_specificity)

  
  summary_df1 <- 
    data.frame(
      hex_id = unique(test_data$hex_id),
      herds_to_sample = herds_to_sample1,
      test_sensitivity = test_sensitivity,
      test_specificity = test_specificity,
      between_herd_prevalence = between_herd_prevalence,
      within_herd_prevalence = within_herd_prevalence,
      overall_animal_prevalence = overall_animal_prevalence,
      ppv = ppv1,
      npv = npv1,
      count_all_herds = nrow(test_data),
      count_all_herds_gr20 = nrow(test_data[test_data$hs22 > 20, ]),
      count_all_herds_gr100 = nrow(test_data[test_data$hs22 > 100, ]),
      count_all_animals = sum(test_data$hs22)
    )
  
  
  return(summary_df1)

}



simulate_btv_sampling_function <- function(between_herd_prevalence, 
          within_herd_prevalence,
          hexagon_reference,
          test_sensitivity,
          test_specificity) 
{
   library(actuar) 
  test_data <- d3[d3$hex_id == hexagon_reference , ]#& d3$hs22 > 20, ] 
  herd_sizes <- test_data$hs22
    n_herds <- length(herd_sizes)
    
    ## simulated infected herds based on between herd prevalence -----------------------------------------------------
    test_data$herd_inf_or_not <-  rbinom(n_herds, 1, between_herd_prevalence)
    
    ## if herds infected simulate infected cattle within herd based on within herd prevalence -----------------------
    herd_sizes <- test_data$hs22[test_data$herd_inf_or_not == 1]
    test_data$n_infected_cattle <- 0 # default zero
    n_inf <- mapply(repeat_rbinom, herd_sizes, within_herd_prevalence)
    test_data$n_infected_cattle[test_data$herd_inf_or_not == 1] <- n_inf
    test_data$n_uninfected_cattle <- test_data$hs22 - test_data$n_infected_cattle
    
    ## now sample herds and cattle for testing ---------------------------------------------------------------------
    
    ### select herds with greater than 100 cattle present ---------------------------------------
    big_herds <- test_data[test_data$hs22 > 100, ]
    
    ### how many herds and animals to sample? ------------------------------------------------------
    
    sampling_advised <- ss_fun(n_herds = nrow(big_herds), vector_herd_sizes = big_herds$hs22,
                               between_herd_prevalence = between_herd_prevalence,
                               within_herd_prevalence = within_herd_prevalence,
                               test_sens = test_sensitivity
    )
    
    herds_to_sample <- sampling_advised[[1]]$nsample # count of herds to sample
    big_herds$cattle_to_sample <- sampling_advised[[2]]$nsample # count of animals in herds to sample
    
    ###  sample  herds ----------------------------------------------------------------------------
    test <- seq(1:nrow(big_herds))
    ref_sampled <- sample(test, herds_to_sample, replace = FALSE)
    ref_sampled
    big_herds$tested <- 0
    big_herds$tested[ref_sampled] <- 1
    test_data1 <- big_herds[big_herds$tested == 1, ] ## pull out tested herds
    
    ### sample cattle within selected herds -----------------------------------------------------------
    test_data1$n_infected_cattle_tested <- mapply(sample_animals_function, test_data1$hs22, 
                                                  test_data1$n_infected_cattle,
                                                  test_data1$n_uninfected_cattle,
                                                  test_data1$cattle_to_sample 
    )
    
    test_data1$n_uninfected_cattle_tested <- test_data1$cattle_to_sample - test_data1$n_infected_cattle_tested
    
    # simulate test results ---------------------------------------------------------------------------------------
    prob_false_pos <- 1 - test_specificity

    test_data1$true_positives <- rbinom(nrow(test_data1), test_data1$n_infected_cattle_tested, test_sensitivity)
    test_data1$false_positives <- rbinom(nrow(test_data1), test_data1$n_uninfected_cattle_tested, prob_false_pos)
    test_data1$false_negatives <- test_data1$n_infected_cattle_tested - test_data1$true_positives
    test_data1$true_negatives <- test_data1$n_uninfected_cattle_tested - test_data1$false_positives

    test_data1$sum_pos_in_herd <- test_data1$true_positives + test_data1$false_positives
    test_data1$any_pos_in_herd <- test_data1$sum_pos_in_herd > 0 
    test_data1$true_pos_herd <- test_data1$true_positives > 0
    test_data1$false_pos_herd <- test_data1$true_positives == 0 & test_data1$false_positives > 0
    
    # test_data1 = tested herds
    # test_data = all herds > 20 animals in hexagon
    summary_df <- 
      data.frame(
        count_all_herds = nrow(test_data),
        count_infected_herds = sum(test_data$herd_inf_or_not),
        count_infected_animals = sum(test_data$n_infected_cattle),
        mean_prop_infected_animals_in_infected_herds = mean(test_data$n_infected_cattle[test_data$herd_inf_or_not == 1] /
                                      test_data$hs22[test_data$herd_inf_or_not == 1] ),
        
        count_all_herds_tested = nrow(test_data1), ## check
        count_infected_herds_tested = sum(test_data1$herd_inf_or_not),
        count_all_animals_in_tested_herds = sum(test_data1$hs22),
        count_infected_animals_in_tested_herds = sum(test_data1$n_infected_cattle),
        
        count_infected_cattle_tested = sum(test_data1$n_infected_cattle_tested),
        count_uninfected_cattle_tested = sum(test_data1$n_uninfected_cattle_tested),
        
        min_cattle_in_herd_sampled = min(test_data1$cattle_to_sample),
        median_cattle_in_herd_sampled = median(test_data1$cattle_to_sample),
        max_cattle_in_herd_sampled = max(test_data1$cattle_to_sample),
        sum_cattle_in_herd_sampled = sum(test_data1$cattle_to_sample),
        
        count_true_pos_cattle = sum(test_data1$true_positives),
        count_false_pos_cattle = sum(test_data1$false_positives),
        count_true_neg_cattle = sum(test_data1$true_negatives),
        count_false_neg_cattle = sum(test_data1$false_negatives),
        
        mean_prop_infected_animals_in_tested_herds = mean(test_data1$n_infected_cattle[test_data1$herd_inf_or_not == 1] /
                                                      test_data1$hs22[test_data1$herd_inf_or_not == 1] ),
        count_apparent_pos_herds = sum(test_data1$any_pos_in_herd),
        count_apparent_pos_animals = sum(test_data1$true_positives + test_data1$false_positives),
        count_false_pos_herds = sum(test_data1$any_pos_in_herd & test_data1$herd_inf_or_not == 0),
        count_true_pos_herds = sum(test_data1$any_pos_in_herd & test_data1$herd_inf_or_not == 1),
        count_false_neg_herds = sum(!test_data1$any_pos_in_herd & test_data1$herd_inf_or_not == 1),
        count_true_neg_herds = sum(!test_data1$any_pos_in_herd & test_data1$herd_inf_or_not == 0),
        mean_prop_apparent_positive_animals_in_tested_herds = mean(test_data1$sum_pos_in_herd[test_data1$any_pos_in_herd] / 
                                                      test_data1$cattle_to_sample[test_data1$any_pos_in_herd])
      )
    
  # if zero test positives
    summary_df$mean_prop_apparent_positive_animals_in_tested_herds[summary_df$count_apparent_pos_herds == 0] <- 0
  # if zero infected tested
    summary_df$mean_prop_infected_animals_in_tested_herds[summary_df$count_infected_herds_tested == 0] <- 0
    summary_df$prop_infected_herds <- summary_df$count_infected_herds / summary_df$count_all_herds
    summary_df$prop_apparent_pos_herds <- summary_df$count_apparent_pos_herds / summary_df$count_all_herds_tested
    ### this is a problem above
    tab2 <- reshape2::melt(summary_df)
    tab2$value <- round(tab2$value, digits = 2)
   tab2$variable <- as.character(tab2$variable)
   tab3 <- merge(tab2, lookup1, by.x = "variable", by.y = "variable")
   tab3 <- tab3[order(tab3$order), ]
   tab3 <- tab3[tab3$include == "y",]
   tab3 <- tab3[ , c("variable1", "value")]
   
   rownames(tab3) <- NULL
   colnames(tab3) <- c("variable", "value")
    return(tab3)
  }

