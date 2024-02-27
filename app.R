# set working directory
# setwd("~/Documents/DATA3888/Capstone Assignment/shinyETHAN")

# Load packages ----
load("data/KA7_app.RData")
source("extractHelper.R")

library(shiny)
library(shinyjs)
library(DT)
library(caret)
library(survival)
library(glmnet)


# User interface ----
# User interface ----
ui <- fluidPage(useShinyjs(),
  
  titlePanel("A7 KT Optimizer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      p("Please select a valid CSV file containing gene expression data sourced from a renal allograft biopsy of a patient."),
      helpText("The shiny app folder contains 2 sample patient csv files to demonstrate the functionalities of our app."),
      br(),
      fileInput("patient_file", label = "Select Patient Data", accept = ".csv", 
                placeholder = "patient_biopsy.csv"),
      em("Example of valid input:"),
      tableOutput("demoData"),
      helpText("The data must contain 2 columns - the gene ID/probe ID (or any other specific identifier)
               and the respective gene expression value for that gene.")
    ),
    
    mainPanel(
      
      tabsetPanel(
        selected = "Info", type = "pills",
        tabPanel(
          "Info", 
          h1("Welcome to A7 KT Optimizer"),
          helpText("Warning: This app is designed to be used by medical practitioners"),
          h4("The A7 KT Optimizer is designed to provide information for doctors in order to make an informed judgement
          on whether to go through with a kidney transplant or not."),
          em("This app has been designed to combat the scarcity of avaliable kidneys for transplantation, 
          and ensure that all transplants have long term success."),
          h4("The app has three goals:"),
          br(),
          h4(div("1. Minimize risk of graft failure", style = "color:DarkRed"), align = "left"),
          h4(div("2. Maximize long term survival", style = "color:SeaGreen"), align = "center"),
          h4(div("3. Optimize kidney allocation", style = "color:DarkGoldenRod"), align = "right"),
          br(),
          br(),
          h5(span("Goal 1", style = "color:DarkRed"), "is achieved through the 'Survival curve' tab."),
          actionButton("p1Details", "Click to show more details"), align = "center",
          hidden(
            div(id='textdiv1',
                textOutput("p1text"))
            ),
          br(),
          h5(span("Goal 2", style = "color:SeaGreen"), "is achieved through the 'Dosage Classifier' tab."),
          actionButton("p2Details", "Click to show more details"),
          hidden(
            div(id='textdiv2',
                textOutput("p2text"))
            ),
          br(),
          br(),
          h5("You will now be able to infer if this kidney transplant is worth preceding with."),
          h5("Your judgement, along with the information gained through this app, is vital in", 
             span("optimizing the process of kidney allocation", style = "color:DarkGoldenRod"), 
             "and ensuring that all kidneys find a suitable long term home.")
          ),
        
                  
        tabPanel(
          "Survival Curve",
          br(),
          helpText("Below is the survival curve for graft failure of sample population."),
          p("Once you enter a new patient's biopsy data, the predicted survival curve for the new patient
                   will be displayed (in red) on top of the survival curve of the population."),
          p("This should allow you to visualise the predicted graft survival of your patient, and compare it with
            a broader population"),
          checkboxInput("hidePop", "Hide population survival curve", value = FALSE),
          plotOutput("plotFinal"),
          textOutput("testText")),
        
        tabPanel(
          "Dosage Classifier",
          br(),
          p("This classifier aims to make predictions on the dosage of immunosuppressants required by the patient
            to retain stable graft function. It will classify the patient as either needing <10mg/day or >10mg/day."),
          helpText("As stated on the 'info' page, patients in the <10mg/day class 
          are argued to have a higher chance of long term survival compared to the >10mg/day class,
          due to the implications of a higher dose of these drugs, which drastically increase the chance for
                   infection (Duncan & Wilkes, 2005)." ),
          strong(htmlOutput("immuneClass")),
          
          br(),
          br(),
          br(),
          helpText("A more thorough look at the effect of immunosupressants on survival can be seen in the figure below."),
          actionButton("imageBut", "Show figure"),
          hidden(
            div(id='imagediv',
                img(src = "figure3.png", height = 320, width = 360),
                p("As you can see, the survival of the group that was weaned off of the drugs is far better than
                  that of the group that maintained the drugs. This illustrates the effect of immunosupressants on survival"),
                helpText("Kaplan-Meier curve of the outcomes in immunosuppressant weaning and maintaining groups. 
                Adapted from “Weaning Immunosuppressant in Patients with Failing Kidney Grafts and The Outcomes: 
                         A Single-Center Retrospective Cohort Study” by Ryu et al, 2020,", 
                         em("Scientific Reports, 10, p. 6425."))
            )
          )
          
        )
                  
      )
      
    )
    
  )
  
)


# Server logic
server <- function(input, output) {
  
  output$demoData <- renderTable({ ## displays example input
    gene_symbol <- c("ABCD9", "AKAP2", "WOW202", "GENE1", "etc...")
    expression_value <- c("4.53", "2.13", "1.23", "8.75", "etc...")
    df <- data.frame("Gene Symbol" = gene_symbol, "Expression Value" = expression_value, 
                     check.names = FALSE)
    df
  })
  
  return_patient_data <- reactive({ ## turns input into data frame :)
    if (is.null(input$patient_file)) return(NULL)
    read.csv(file = input$patient_file$datapath)
  })
  
  patient_survival_curve <- reactive({ ## gets input data frame and tries to plot survival curve
    req(return_patient_data())
    patient_data <- return_patient_data()
    patient_data <- extractData(patient_data, "s")
    survfit(surv_coxnet, x=survivalX, y=survivalY.surv, s=survivalS, newx=patient_data$Expression.Value, data=survdata)
  })
  
  output$plotFinal <- renderPlot({ ## plots both the pop survival curve and the patient
    if (is.null(return_patient_data()) & !input$hidePop){
      plot(model_survival_curve, xlab = "Years post transplant", ylab = "Probability of survival",
           xaxs = "S", xscale = 365.25, yscale = 100)
    }else if (!is.null(return_patient_data()) & !input$hidePop){
      plot(model_survival_curve, xlab = "Years post transplant", ylab = "Probability of survival",  
           xaxs = "S", xscale = 365.25, yscale = 100)
      lines(patient_survival_curve(), col = "red")
    }else if (!is.null(return_patient_data()) & input$hidePop){
      plot(patient_survival_curve(), xlab = "Years post transplant", ylab = "Probability of survival", col = "red",
           xaxs = "S", xscale = 365.25, yscale = 100)
    }else if (is.null(return_patient_data()) & !input$hidePop){
      
    }
  })
  
  immuneClassFac <- reactive({
    req(return_patient_data())
    patient_data <- return_patient_data()
    patient_data <- extractData(patient_data, "i", immuneIDX)
    # predict(immuneKNN, t(patient_data))
    icNum <- as.numeric(predict(immuneKNN, t(patient_data)))
    if(icNum == 1){ #<10
      "<10mg/day"
    }else{
      ">10mg/day"
    }
  })
  
  output$immuneClass <- renderText({
      paste("This patient requires ","<font color=\"#FF0000\"><b>", immuneClassFac(), 
            " </b></font> of immunosuppressants to retain stable graft function.")
  })
  
  
  observeEvent(input$p1Details, {
    toggle('textdiv1')
    output$p1text <- renderText({"Using an elastic net cox regression model, the app
      plots the predicted survival curve for the graft failure of a patient. Through visualising the
      survival curve of a patients graft, you will be able to infer the success of the the
      transplant over time."})
  })
  
  observeEvent(input$p2Details, {
    toggle('textdiv2')
    output$p2text <- renderText({"By creating a binary predictor for the required immunosuppressant dosage the doctor will be
            able to see whether the patient requires >/< 10mg of immunosuppressants per day to retain stable 
            graft function. Patients that fall into the '<' category, have a higher chance of long term survival,
            due to the implications of a higher dose of these drugs, which drastically increase the chance for
            infection (Duncan & Wilkes, 2005)."})
  })
  
  observeEvent(input$imageBut, {
    toggle('imagediv')
  })
  
}

# Run the app
shinyApp(ui, server)