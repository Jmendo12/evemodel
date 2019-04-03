library(shiny)
library(tidyverse)
library(ggplot2)
library(ape)
library(ggtree)

source("../scripts/EVEcore.R")
source("../scripts/twoThetaTest.R")

# load "dataSets"
load("data.RData")

dataSets$salmon20$shiftSpecies <- c("Salp","Okis","Omyk","Ssal")
dataSets$sim$shiftSpecies <- c("D","E")

lapply(dataSets, function(dataSet){
  dataSet$g_phylo <- ggtree(dataSet$tree)
  dataSet$orderedSpcs <- 
    dataSet$g_phylo$data %>% filter(isTip) %>% arrange(y) %>% with(label)
  dataSet$initParams <- initialParamsTwoTheta(dataSet$gene.data,colSpecies = colnames(dataSet$gene.data),shiftSpecies = dataSet$shiftSpecies)
  return(dataSet)
}) -> dataSets


getMeanSigmaTwoTheta <- function(theta1, theta2, sigma2, alpha, beta, tree, shiftSpecies, colSpecies)
{
  isThetaShiftEdge <- 1:Nedge(tree) %in% getEdgesFromMRCA(tree, tips = shiftSpecies)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  N <- Nedge(tree)
  expression.var <- calcExpVarOU(tree, thetas = ifelse(isThetaShiftEdge,theta2,theta1), 
                                 alphas = rep(alpha,N), sigma2s = rep(sigma2,N),
                                 rootVar = sigma2/(2*alpha), rootE = theta1)
  
  covar.matrix <- calcCovMatOU(tree, alphas = alpha, evol.variance = expression.var$evol.variance)
  
  expanded.matrix <- expandECovMatrix(expected.mean = expression.var$expected.mean,covar.matrix,
                                      sigma2, alpha, beta, index.expand)
  
  return(list(mean = expanded.matrix$expected.mean, sigma = expanded.matrix$cov.matr))
}

fitTwoThetasVarFix <- function(tree, gene.data.row, initParams, shiftSpecies, colSpecies, isParamFixed){
  isThetaShiftEdge <- 1:Nedge(tree) %in% getEdgesFromMRCA(tree, tips = shiftSpecies)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  # Calculate max alpha based on the proportion of variance from covariance
  alphaMax <- -log(.01) / min(tree$edge.length[tree$edge[,2] <= Ntip(tree)])
  
  paramIterations <- numeric()
  
  # function to optimize
  LLPerGeneTwoTheta <- function(par)
  {
    params <- initParams
    params[!isParamFixed] <- par
    ll <- logLikTwoTheta( theta1 = params[1], theta2 = params[2], sigma2 = params[3],
                          alpha = params[4], beta = params[5],
                          tree = tree, isThetaShiftEdge = isThetaShiftEdge, 
                          gene.data.row = gene.data.row, index.expand = index.expand)
    paramIterations <<- rbind(paramIterations,c(params,ll=ll))
    return(-ll)
  }
  
  res <- optim( par = initParams[!isParamFixed], fn = LLPerGeneTwoTheta, method = "L-BFGS-B",
           lower = c(-Inf, -Inf, 1e-10, 1e-10, 1e-10)[!isParamFixed], 
           upper = c(Inf, Inf, Inf, alphaMax, Inf)[!isParamFixed])
  
  return(list(optimRes = res, paramIterations=paramIterations))
}



ui <- fixedPage(
  fixedRow(
    column( 3,
      selectInput("selData",label = "Dataset:",choices = names(dataSets)),
      selectInput("selGene",label = "Gene:",choices = rownames(dataSets[[1]]$gene.data)),
      wellPanel( # parameter inputs
        tags$table(
          tags$tr(
            tags$td("theta1"), tags$td( numericInput("theta1",NULL,1)), tags$td( checkboxInput("theta1Fix", NULL))
          ),
          tags$tr(
            tags$td("theta2"), tags$td( numericInput("theta2",NULL,1)), tags$td( checkboxInput("theta2Fix", NULL))
          ),
          tags$tr(
            tags$td("alpha"), tags$td( numericInput("alpha",NULL,0.5)), tags$td( checkboxInput("alphaFix", NULL))
          ),
          tags$tr(
            tags$td("sigma2"), tags$td( numericInput("sigma2",NULL,1)), tags$td( checkboxInput("sigma2Fix", NULL))
          ),
          tags$tr(
            tags$td("beta"), tags$td( numericInput("beta",NULL,1)), tags$td( checkboxInput("betaFix", NULL))
          )
        )
      ),
      actionButton("runOptim", "Optimize!")
    ),
    
    column(9,
      tags$table(
        tags$tr(
          tags$td( # phylogeny
            plotOutput("treePlot",width = 250,height = 250)
          ), 
          tags$td( # expression data
            plotOutput("exprPlot",width = 250,height = 250)
          ),
          tags$td(# covariance matrix
            plotOutput("covMatPlot",width = 250,height = 250)
          )
        )
      ),
      rbokehOutput("optimPlot")
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  optimReact <- eventReactive(input$runOptim, {
    dataSet <-  reactDataset()
    gene.data.row <- reactDataRow()
    initParams <- unlist(reactParams())
    
    fitTwoThetasVarFix(tree = dataSet$tree,
                       gene.data.row = gene.data.row,
                       initParams = initParams,
                       shiftSpecies = dataSet$shiftSpecies,
                       colSpecies = colnames(dataSet$gene.data),
                       isParamFixed = reactParamFix())
  })
  
  # update gene selection when changing dataset
  observe({
    dataSet = reactDataset()
    geneIDs <- rownames(dataSet$gene.data)
    updateSelectInput(session, "selGene", choices = geneIDs, selected = geneIDs[1])
  })
  
  # update initial params when changing gene
  observe({
    dataSet <- reactDataset()
    geneID <- input$selGene

    # TODO: fix so that this never happens:
    if(!(geneID %in% rownames(dataSet$gene.data)))
      geneID <- rownames(dataSet$gene.data)[1]

    for( param in c("theta1","theta2","alpha","sigma2","beta")){
      updateNumericInput(session, param, value = dataSet$initParams[geneID,param])
    }
  })
  
  reactParams <- reactive({
    list(
      theta1 = input$theta1,
      theta2 = input$theta2,
      alpha = input$alpha,
      sigma2 = input$sigma2,
      beta = input$beta
    )
  })
  
  reactParamFix <- reactive({
    c(
      theta1 = input$theta1Fix,
      theta2 = input$theta2Fix,
      alpha = input$alphaFix,
      sigma2 = input$sigma2Fix,
      beta = input$betaFix
    )
  })
  
  reactMeanSigma <- reactive({
    params <- reactParams() # params = as.list(dataSet$initParams[1,])
    dataSet <- reactDataset()
    meanSigma <- do.call(getMeanSigmaTwoTheta, 
                         args = c(params, list( tree = dataSet$tree, shiftSpecies = dataSet$shiftSpecies, 
                                                colSpecies = colnames(dataSet$gene.data)) ) )
    
    idx <- order(match(colnames(dataSet$gene.data),dataSet$orderedSpcs))
    meanSigma$mean <- meanSigma$mean[idx]
    meanSigma$sigma <- meanSigma$sigma[idx,idx]
    return(meanSigma)
  })
  
  reactDataset <- reactive({
    dataSet <- dataSets[[input$selData]]
    return(dataSet)
  })
  
  reactDataRow <- reactive({
    dataSet <- reactDataset()
    geneID <- input$selGene
    
    # TODO: fix so that this never happens:
    if(!(geneID %in% rownames(dataSet$gene.data))) 
      geneID <- rownames(dataSet$gene.data)[1]
    
    gene.data.row <- dataSet$gene.data[geneID, ]
    return(gene.data.row)
  })
  
  output$treePlot <- renderPlot({
    dataSet <- reactDataset()
    dataSet$g_phylo + 
      geom_tiplab(align = T) +
      scale_y_continuous(limits = c(1,Ntip(dataSet$tree)),expand=expand_scale(add=1)) +
      scale_x_continuous(expand=expand_scale(mult=c(0.01,0.21))) + 
      xlab(" ") + 
      theme_classic() +
      theme(axis.line.y = element_blank(),axis.text.y = element_blank(), axis.ticks.y = element_blank() )
  })
  
  output$exprPlot <- renderPlot({
    dataSet <- reactDataset()
    gene.data.row <- reactDataRow()
    tibble( x=gene.data.row, spc=names(gene.data.row) ) %>% 
      left_join( select(dataSet$g_phylo$data,label,y), by=c("spc"="label")) %>% 
      ggplot( aes(x=x,y=y) ) +
      geom_point() +
      scale_y_continuous(limits = c(1,Ntip(dataSet$tree)),expand=expand_scale(add=1)) +
      xlab("expression level") +
      theme_bw()+
      theme(axis.line.y = element_blank(),axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), axis.title.y = element_blank() )
  })
  
  output$covMatPlot <- renderPlot({
    meanSigma <- reactMeanSigma()
    sigma <- meanSigma$sigma
    colnames(sigma) <- paste0(sprintf("%03i",1:nrow(sigma)),"_",colnames(sigma))
    rownames(sigma) <- colnames(sigma)
    reshape2::melt(sigma) %>% 
      ggplot( aes(x=Var1,y=Var2,fill=value)) + 
      geom_tile() + 
      theme_bw()+
      theme(axis.line.y = element_blank(),axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), axis.title.y = element_blank() )
    
  })
  
  # output$optimPlot <- renderPlot({
  #   x <- optimReact()
  #   as.tibble(x$paramIterations) %>%
  #     mutate(iteration=1:n()) %>% 
  #     gather(key = "param",value = "value", -iteration) %>% 
  #     ggplot( aes(x=iteration,y=value)) + geom_point() + geom_line() + facet_grid( param ~ .,scales = "free_y")
  # })
  
  output$optimPlot <- renderRbokeh({
    x <- optimReact()
    myData <- as.tibble(x$paramIterations) %>% mutate(iteration=1:n())
    
    tools <- c("pan", "wheel_zoom", "box_zoom", "box_select", "reset")
    
    lapply(colnames(x$paramIterations), function(par){
      figure(width = 800, height = 100, tools = tools,
             xlab = "iteration", ylab = par) %>%
        ly_points("iteration", par, data = myData,
                  color = par, size = 5, legend = FALSE)
      
    }) -> splom_list
    
    grid_plot(splom_list, ncol = 1, same_axes = c(TRUE,FALSE))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

