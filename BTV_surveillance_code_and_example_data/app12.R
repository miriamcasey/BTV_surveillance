library(shiny)
library(kableExtra)

source("code_for_IVJ_BTV_functions3.R")
ui <- fluidPage(
  
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1")
  ),
  
  # Ensure mobile responsiveness and scalable UI
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      body {
        font-size: 16px;
      }
      .form-control {
        min-width: 100%;
      }
    "))
  ),
  
  
  titlePanel("Two-stage surveillance in 20km radius (Ireland)"),
  
  # Clear explanatory text
  tags$div(style = "margin-bottom: 20px;",
           p("This app supports the Irish Veterinary Journal article:"),
           tags$strong("Preparation for a potential outbreak of bluetongue virus in Ireland."),
           p("It allows users to explore how diagnostic test sensitivity, specificity, and infection prevalence affect the required sample size and expected surveillance outcomes."),
           p("The model uses cattle populations from 20km-radius hypothetical temporary control zones in Ireland, applying a two-stage sampling design. Sample sizes are calculated using the ",
             tags$code("rsu.sssep.rs2st"),
             " function from the ",
             tags$code("Epitools"),
             " package. Monte Carlo simulation is used to estimate outcomes. Please see our Methods tab for a longer description and our paper for a full description.")
  ),
  sidebarLayout(
    sidebarPanel(
      
      sliderInput(
        inputId = "between_herd_prevalence",
        "Between herd prevalence",
        min = 0.025,
        max = 0.7,
        value = 0.05,
        step = 0.01
      ),
      
      sliderInput(
        inputId = "within_herd_prevalence",
        "Within herd prevalence",
        min = 0.05,
        max = 1,
        value = 0.3,
        step = 0.05
      ),

      sliderInput(
        inputId = "test_sensitivity",
        "Diagnostic test sensitivity",
        min = 0.8,
        max = 1,
        value = 0.99,
        step = 0.01
      ),
      
      sliderInput(
        inputId = "test_specificity",
        "Diagnostic test specificity",
        min = 0.95,
        max = 1,
        value = 1,
        step = 0.001
      ),
      numericInput(
        inputId = "hex_id",
        label = "Temporary Control Zone (20km radius) 1 - 8",
        value = 1,
        min = 1,
        max = 8,
        step = 1
      )   
,

numericInput(
  inputId = "chosen_reps",
  label = "Times to repeat simulation",
  min = 1,
  max = 100,
  value = 100,
  step = 1
)   
,
      actionButton("goButton", "Generate Results")
    ),
    
mainPanel(
  tabsetPanel(
    tabPanel("Results",
             plotOutput("tp_fp_plot", height = "500px"),
             htmlOutput("tableOutput1"),
             htmlOutput("tableOutput")
    ),
    tabPanel("Methods",
             fluidPage(
               h3("Methods"),
               
               h4("Surveillance design"),
               p("We designed targeted surveillance for a 20km radius temporary control zone (TCZ) after the detection of an initial case of bluetongue (BT). This surveillance would be in addition to (1) enhanced passive surveillance for clinical and postmortem signs of BT (2) post import testing and (3) targeted surveillance based on the likely sites of windborne midge dispersion."),
               p("We design a two stage (herd and animal level) surveillance plan, considering expected between-herd and within-herd prevalence of BT. For logistical reasons, we aim to minimise the numbers of herds requiring visits. Our geographical sampling unit is the 20km radius TCZ."),
               
               h4("Irish cattle data"),
               p("Herd profiles were estimated using the DAFM’s Animal Identification and Movements (AIM) database for Jan, May, and Sept 2022. Herd sizes were averaged over these dates. DAFM’s LPIS was used to delineate land area and herd centroids. 5% of herds not recorded in LPIS were mapped using Electoral Division data."),
               
               h4("Definition of TCZs"),
               p("A total of 1,064 overlapping 20km radius circles were created to cover the island of Ireland. These TCZs allow for localized surveillance and estimation of herds present based on centroid data."),
               
               h4("Sampling logistics"),
               p("Given that minimising herds visited is logistically more beneficial to the Irish state veterinary service than minimising animals sampled, and also given that larger herds could more easily provide facilities for sampling, we selected herds of greater than 100 cattle in size. Because our surveillance deign was based on herds comprising >100 cattle, to enable the simulation to run with low between herd prevalence (requiring more herds to be sampled), we selected TCZs (n = 862) with > 70 herds comprising more than 100 cattle."),
               
               h4("Sample size estimation"),
               p("We used the “EpiR” package (Stevenson et al., 2024) to estimate the sample size required in each TCZ. The two-stage representative survey design tool in the “EpiR” package allowed consideration of both within- and between- herd prevalence. We calculated how many herds and how many animals from within each herd need to be sampled to be 95% confident of detecting disease at the herd and individual animal level."),
               
               h4("Scenario simulation"),
               p("To explore the effectiveness of our planned surveillance, including levels of false negative and false positive results under different test interpretation conditions, we conducted a simulation study. Static between- and within herd prevalence was assumed. The steps in the simulation were as follows.
1.	For each TCZ, infected herds were simulated based on the expected between herd prevalence and the total count of cattle herds.
2.	For each potential TCZ, herds with 100 cattle or more were selected.
3.	The number of herds required to be sampled, based on the “EpiR” based sample size calculation was randomly selected from the herds defined in point 2.
4.	If any of the selected herds had been simulated as infected herds, infected animals within that herd were simulated based on expected within herd prevalence and herd size.
5.	From each selected herd, the count of animals required to be sampled, based on sample size calculation, was randomly selected.
6.	Test results in selected animals were simulated based on the animal’s simulated true infection status, and test sensitivity and specificity. 
7.	Simulated true and apparent prevalence was compared, and the levels of false positives and false negatives under different test and prevalence conditions were reviewed.
"),
               
               h4("References"),
               tags$ul(
                 tags$li("Stevenson et al., 2024. Package ‘epiR’."),
                 tags$li("Tratalos et al., 2020. Spatial and network characteristics of Irish cattle movements."),
                 tags$li("Zimmermann et al., 2016. The Irish land-parcels identification system (LPIS).")
               )
             )
    )
    
  )
)

  )
)

server <- function(input, output, session) {
  
  # Trigger this only when goButton is clicked
  data_reactive <- eventReactive(input$goButton, {

    xa <- replicate(input$chosen_reps, simulate_btv_sampling_function (
      between_herd_prevalence = input$between_herd_prevalence,
      within_herd_prevalence = input$within_herd_prevalence,
      hexagon_reference = input$hex_id,
      test_sensitivity = input$test_sensitivity,
      test_specificity = input$test_specificity),
                    simplify = FALSE)
    
    x1 <- rbindlist(xa)
    
    x1$reps <- rep(1:input$chosen_reps, each = 27)
    
    test <- data.table(x1)
    test3 <- test[,   .(
      .N,
      min = round(min(value), digits = 2),
      lq = round(quantile(value, 0.25), digits = 2),
      mean = round(mean(value), digits = 2),
      median = round(median(value), digits = 2),
      uq = round(quantile(value, 0.75), digits = 2),
      max = round(max(value), digits = 2),
      count_zeros = sum(value == 0)
    ),
    by = .(variable)]
    test3
    lookup1
    test4 <- merge(test3, lookup1, by.x = "variable", by.y = "variable1")
    test4 <- test4[order(test4$order), ]
    test4$text <- paste0(test4$median, " (", test4$lq, "-", test4$uq, ")")
    test4$round
    test4$text[test4$round == "y"] <- paste0(round(test4$median[test4$round == "y"] ),
                                             " (",
                                             round(test4$lq[test4$round == "y"] ), 
                                             "-", round(test4$uq[test4$round == "y"] ), ")")
    
    
    
    test5 <- test4[ , c("variable", "text")]
    colnames(test5) <- c("Variable", "Median (IQR)")
    
    
    tab1 <- make_first_table(between_herd_prevalence = input$between_herd_prevalence,
                             within_herd_prevalence = input$within_herd_prevalence,
                             hexagon_reference = input$hex_id,
                             test_sensitivity = input$test_sensitivity,
                             test_specificity = input$test_specificity)
    
    tab2 <- reshape2::melt(tab1)
   # tab2$value <- round(tab2$value, digits = 4)
    lookup1
    colnames(tab2)
    colnames(lookup1)
    tab3 <- merge(tab2, lookup1, by.x = "variable", by.y = "variable")
    tab3 <- tab3[order(tab3$order), ]
    tab3 <- tab3[ , c("variable1", "value")]
    colnames(tab3) <- c("Variable", "Value")
    rownames(tab3) <- NULL

    
    list(test5, x1, tab3) ## both dataframes in list

  })
  
  output$tableOutput1 <- renderUI({
    req(data_reactive())  # Ensure data is available
    
    df <- data_reactive()[[3]]
    
    # Identify the column to format (change "Value" to your actual column name)
    column_to_format <- "Value"  # <-- CHANGE THIS to your column name

    # Apply conditional formatting: integers as-is, decimals with 3 decimal places
    df[[column_to_format]] <- sapply(df[[column_to_format]], function(x) {
      if (is.na(x)) return(NA)  # keep NA values
      if (x == floor(x)) {
        format(x, nsmall = 0, scientific = FALSE)
      } else {
        format(round(x, 3), nsmall = 3, scientific = FALSE)
      }
    })
    
    # Now convert to HTML table with kable
    kable_html <- knitr::kable(df,  format = "html",
                               caption = "Table 1: Inputs") %>%
      kable_styling("striped", full_width = F)
    
    HTML(kable_html)
  })
  
  

  # Render the HTML table
  output$tableOutput <- renderUI({
    req(data_reactive())  # ensure it's available

    # Convert to HTML table
    df <- data_reactive()[[1]]

    kable_html <- knitr::kable(df, format = "html", 
                               caption = "Table 2: Simulation outputs"
                               
                               ) %>%
      kable_styling("striped", full_width = F)
    
     HTML(kable_html)



    #HTML(knitr::kable(df, format = "html", table.attr = "class='table table-bordered'"))
  })
  
  output$tp_fp_plot <- renderPlot({
    df <- data_reactive()[[2]]
    req(nrow(df) > 0)
    
    x2 <- df[df$variable %in% c("True positive herds", "False positive herds", "False negative herds", "Infected herds tested"), ]
    
    x2$variable <- factor(x2$variable, levels = c("False positive herds", "True positive herds", "False negative herds","Infected herds tested"))
    max_value <- max(x2$value)
    max_value_breaks <- round(max_value / 10)
    
    # Define color palette if not already
    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
    
    ggplot(x2) +
      geom_bar(
        aes(y = value,
            fill = variable),
        dodge.width = 0.9,
        shape = 21,
        size = 3,
        alpha = 0.7
      ) + facet_wrap(variable ~ .) + coord_flip() +
      xlab("Count simulations") + ylab("Count herds") +
  #    scale_y_continuous(breaks = seq(0, max_value, by = max_value_breaks)) +
      scale_colour_manual(values = cbPalette) +
      scale_fill_manual(values = cbPalette) +
      theme_minimal(base_size = 25) +
      theme(legend.position = "none") +
      coord_flip()
  })
  

}

shinyApp(ui = ui, server = server)