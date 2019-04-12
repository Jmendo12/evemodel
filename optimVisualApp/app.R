library(shiny)
library(tidyverse)
library(ape)
library(ggtree)
library(rbokeh)

source("../scripts/EVEcore.R")
source("../scripts/twoThetaTest.R")

source("dyOptIterPlot.R")

# load "dataSets"
load("data.RData")

# initialize a table to store optimization results
resTblInit <- bind_rows(
  tibble(geneID=character(),theta1=numeric(),theta2=numeric(),sigma2=numeric(),alpha=numeric(),
         beta=numeric(),ll=numeric(), iterations=integer(),convCode=integer(), transform=character()),
  rownames_to_column(as.data.frame(dataSets$salmon20$initParams),var = "geneID"),
  rownames_to_column(as.data.frame(dataSets$sim$initParams),var = "geneID"))

# define the names of the parameters here
paramNames = c("theta1", "theta2", "sigma2", "alpha", "beta")

parTransFuns <- list(
  none = list(
    trans = function(par){par},
    untrans = function(par){par}
  ),
  log = list(
    trans = function(par){
      setNames(ifelse( names(par) %in% c("sigma2", "alpha", "beta"), log(par), par), names(par))
    },
    untrans = function(par){
      setNames(ifelse( names(par) %in% c("sigma2", "alpha", "beta"), exp(par), par), names(par))
    }),
  rho = list(
    trans = function(par){
      par["sigma2"] <- par["sigma2"]/(2*par["alpha"])
      par
    },
    untrans = function(par){
      par["sigma2"] <- par["sigma2"]*(2*par["alpha"])
      par
    })
)

# dataSetID <- "salmon20"
# params <- as.list(dataSets[[dataSetID]]$initParams[3,])
getMeanSigmaTwoTheta <- function(params, dataSetID)
{
  tree <- dataSets[[dataSetID]]$tree
  shiftSpecies <- dataSets[[dataSetID]]$shiftSpecies
  colSpecies <- colnames(dataSets[[dataSetID]]$gene.data)
  
  isThetaShiftEdge <- 1:Nedge(tree) %in% getEdgesFromMRCA(tree, tips = shiftSpecies)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  N <- Nedge(tree)
  expression.var <- calcExpVarOU(tree, thetas = ifelse(isThetaShiftEdge,params$theta2,params$theta1), 
                                 alphas = rep(params$alpha,N), sigma2s = rep(params$sigma2,N),
                                 rootVar = params$sigma2/(2*params$alpha), rootE = params$theta1)
  
  covar.matrix <- calcCovMatOU(tree, alphas = params$alpha, evol.variance = expression.var$evol.variance)
  
  expanded.matrix <- expandECovMatrix(expected.mean = expression.var$expected.mean,covar.matrix,
                                      params$sigma2, params$alpha, params$beta, index.expand)
  
  return(list(nodeMean = expression.var$expected.mean,
              nodeVar = expression.var$evol.variance,
              expandedMean = expanded.matrix$expected.mean,
              sigma = expanded.matrix$cov.matr))
}

# curData <- list()
# curData$dataSetID <- "salmon20"
# curData$geneID <- "OG0000053_2"
# curData$params <- as.list(dataSets[[curData$dataSetID]]$initParams[curData$geneID, ])
# curData$gene.data.row = dataSets[[curData$dataSetID]]$gene.data[curData$geneID,]
# 
# dataSetID = curData$dataSetID
# gene.data.row = curData$gene.data.row
# initParams <- unlist(curData$params)
# isParamFixed = c(F,F,F,F,F)
fitTwoThetasVarFix <- function(dataSetID, gene.data.row, initParams, isParamFixed, transFun = parTransFuns$none){
  tree <- dataSets[[dataSetID]]$tree
  shiftSpecies <- dataSets[[dataSetID]]$shiftSpecies
  colSpecies <- colnames(dataSets[[dataSetID]]$gene.data)
  
  isThetaShiftEdge <- 1:Nedge(tree) %in% getEdgesFromMRCA(tree, tips = shiftSpecies)
  
  # match the column species with the phylogeny tip labels
  index.expand <- match(colSpecies, tree$tip.label)
  
  # Calculate max alpha based on the proportion of variance from covariance
  alphaMax <- -log(.01) / min(tree$edge.length[tree$edge[,2] <= Ntip(tree)])
  
  paramIterations <- numeric()
  
  initParams <- transFun$trans(initParams)
  
  # function to optimize
  LLPerGeneTwoTheta <- function(par)
  {
    params <- initParams
    params[!isParamFixed] <- par
    params <- transFun$untrans(params)
    ll <- logLikTwoTheta( theta1 = params[1], theta2 = params[2], sigma2 = params[3],
                          alpha = params[4], beta = params[5],
                          tree = tree, isThetaShiftEdge = isThetaShiftEdge, 
                          gene.data.row = gene.data.row, index.expand = index.expand)
    paramIterations <<- rbind(paramIterations,c(params,ll=ll))
    return(-ll)
  }
  
  res <- optim( par = initParams[!isParamFixed], fn = LLPerGeneTwoTheta, method = "L-BFGS-B",
           lower = transFun$trans(c(-Inf, -Inf, 1e-10, 1e-10, 1e-3))[!isParamFixed], 
           upper = transFun$trans(c(Inf, Inf, Inf, alphaMax, Inf))[!isParamFixed])
  
  res$par <- transFun$untrans(res$par)
  
  return(list(optimRes = res, paramIterations=paramIterations))
}

doSimulate <- function(params, dataSetID){
  meanSigma <- getMeanSigmaTwoTheta(params, dataSetID)
  rmvnorm(n=1, mean=meanSigma$expandedMean, sigma = meanSigma$sigma)
}



ui <- fixedPage(
  fixedRow(
    column( 3,
      selectInput("selData",label = "Dataset:",choices = names(dataSets)),
      selectInput("selGene",label = "Gene:",choices = rownames(dataSets[[1]]$gene.data)),
      selectInput("selTrans",label = "Transform:",choices = names(parTransFuns)),
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
      actionButton("runOptim", "Optimize!"),
      actionButton("updateParams", "Update params"),
      actionButton("simulate", "Simulate")
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
      dyOptIterPlotUI(id = "optimPlot",paramNames = c(paramNames,"ll")),
      dataTableOutput('restable')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  curData <- reactiveValues( # define reactive values and set defaults
    dataSetID = names(dataSets)[1],
    geneID = rownames(dataSets[[1]]$gene.data)[1],
    gene.data.row = dataSets[[1]]$gene.data[1,],
    params = as.list(dataSets[[1]]$initParams[1,])
  )
  
  optRes <- reactiveValues(
    all = resTblInit # currentOpt
  )
  
  observeEvent( input$updateParams, {
    res <- optRes$currentOpt
    if( is.null(res) ) return()
    # update param Inputs
    for( param in names(res$optimRes$par)){
      cat('updateNumericInput(session, inputID="',param,'", value = ',res$optimRes$par[param],'\n',sep = "")
      updateNumericInput(session, param, value = as.vector(res$optimRes$par[param]))
    }
  })
  
  observeEvent( input$simulate, {
    dataSet <- dataSets[[curData$dataSetID]] # dataSet = dataSets$salmon20
    params <- curData$params # params = as.list(dataSet$initParams[1,])
    
    simData <- doSimulate(params, curData$dataSetID)

    # add simdata to current dataSet
    
    # add empty row if simulated data not already exists
    if(!("simulated" %in% rownames(dataSet$gene.data))){
      dataSet$gene.data <- rbind(dataSet$gene.data, simulated=0)
      dataSet$initParams <- rbind(dataSet$initParams, simulated=0)
    }
    # set the data
    dataSet$gene.data["simulated",] <- simData
    dataSet$initParams["simulated",] <- unlist(params)[colnames(dataSet$initParams)]
    
    # save it to the global dataset
    dataSets[[input$selData]] <<- dataSet
    
    # update curData
    curData$geneID <- "simulated"
    curData$gene.data.row <- dataSet$gene.data["simulated",]
    
    # update selGene
    geneIDs <- rownames(dataSet$gene.data)
    updateSelectInput(session, "selGene", choices = geneIDs, selected = "simulated")
  })
  
  observeEvent( input$runOptim, {
    initParams <- unlist(curData$params)
    isParamFixed <- reactParamFix()
    transFun <- parTransFuns[[input$selTrans]]
    
    res <- fitTwoThetasVarFix(
      dataSetID = curData$dataSetID,
      gene.data.row = curData$gene.data.row,
      initParams = initParams,
      isParamFixed = isParamFixed,
      transFun = transFun)
    
    # add result to list
    optimResult <- initParams
    optimResult[!isParamFixed] <- res$optimRes$par
    optimResult <- bind_cols(as.data.frame(t(optimResult)),
                             geneID=curData$geneID, ll= -res$optimRes$value,
                             iterations=res$optimRes$counts[1],convCode=res$optimRes$convergence,
                             transform=input$selTrans)

    # update reactive values
    optRes$all <- bind_rows(optRes$all,optimResult) %>% distinct()
    optRes$currentOpt <- res
  })
  
  # update gene selection when changing dataset
  observe({
    dataSet <- dataSets[[input$selData]]
    geneIDs <- rownames(dataSet$gene.data)
    updateSelectInput(session, "selGene", choices = geneIDs, selected = geneIDs[1])
  })
  
  # update initial params and curData when changing gene
  observe({
    dataSetID <- isolate(input$selData)
    dataSet <- dataSets[[dataSetID]]
    geneID <- input$selGene
    
    if(!(geneID %in% rownames(dataSet$gene.data))){
      cat("!!! selGene and selData not synced !!!\n")
      geneID <- rownames(dataSet$gene.data)[1]
    }

    params <- dataSet$initParams[geneID, ]
    
    for( param in paramNames){
      # NOTE: for some reason the value cannot be named!
      updateNumericInput(session, param, value = as.vector(params[param]))
    }
    
    curData$dataSetID = dataSetID
    curData$geneID = geneID
    curData$gene.data.row = dataSet$gene.data[geneID,]
    curData$params = as.list(params)
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
  throttledParams <- reactParams %>% throttle(500)
  # update curData$params with throttled params
  observe({
    curData$params <- reactParams()
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
    params <- curData$params # params = as.list(dataSet$initParams[1,])
    dataSetID <- curData$dataSetID
    
    meanSigma <- getMeanSigmaTwoTheta( params, dataSetID)
   
    # rearrange covariance matrix to match the species order in the tree
    idx <- order(match(colnames(dataSets[[dataSetID]]$gene.data),dataSets[[dataSetID]]$orderedSpcs))
    meanSigma$sigma <- meanSigma$sigma[idx,idx]
    
    return(meanSigma)
  })
  
  

  output$treePlot <- renderPlot({
    tree <- dataSets[[curData$dataSetID]]$tree
    g_phylo <- dataSets[[curData$dataSetID]]$g_phylo
    alpha <- curData$params$alpha

    x <- seq(from=min(g_phylo$data$x),to=max(g_phylo$data$x),length.out = 200)
    alphaStrip <- tibble( x=x,z=exp(-x*alpha))

    g_phylo + 
      geom_tiplab(align = T) +
      scale_y_continuous(limits = c(0,Ntip(tree)+1),expand=expand_scale(add=0)) +
      scale_x_continuous(expand=expand_scale(mult=c(0.01,0.21))) + 
      xlab(" ") + 
      geom_tile( data = alphaStrip, aes(x=x,fill=z,y=0.3), height=0.4) +
      scale_fill_gradient(low="white",high = "red",limits=c(0,1),guide = F) +
      theme_classic() +
      theme(axis.line.y = element_blank(),axis.text.y = element_blank(), axis.ticks.y = element_blank() )
  })
  
  output$exprPlot <- renderPlot({
    dataSet <- dataSets[[curData$dataSetID]]
    gene.data.row <- curData$gene.data.row
    meanSigma <- reactMeanSigma()
    params <- curData$params
    
    # calculate measured species means
    spcMean <- 
      tibble( x=gene.data.row, spc=names(gene.data.row) ) %>% 
        group_by(spc) %>% summarise(spcMean = mean(x))
    
    spcMeanVar <-
      tibble( spc = dataSet$tree$tip.label,
            spcExpMean=meanSigma$nodeMean[1:Ntip(dataSet$tree)],
            spcExpVar=meanSigma$nodeVar[1:Ntip(dataSet$tree)]) %>%
      left_join( select(dataSet$g_phylo$data,label,y), by=c("spc"="label")) %>%
      left_join( spcMean, by="spc" ) %>% 
      mutate( x1 = spcExpMean - 2*sqrt(spcExpVar), x2 = spcExpMean + 2*sqrt(spcExpVar)) %>% 
      mutate( tau = spcExpVar*params$beta ) %>% 
      mutate( tau_x1 = spcMean - 2*sqrt(tau), tau_x2 = spcMean + 2*sqrt(tau))

    
    tibble( x=gene.data.row, spc=names(gene.data.row) ) %>% 
      left_join( select(dataSet$g_phylo$data,label,y), by=c("spc"="label")) %>% 
      ggplot() +
      geom_segment( data=spcMeanVar, aes( y=y,yend=y,x=x1,xend=x2), size=3, color="red", alpha=0.5) +
      geom_point( data=spcMeanVar, aes(x=spcExpMean,y=y), color="red", size=4,alpha=0.5) +
      geom_segment( data=spcMeanVar, aes( y=y,yend=y,x=tau_x1,xend=tau_x2), size=1, color="blue", alpha=0.5) +
      geom_point( data=spcMeanVar, aes(x=spcMean,y=y), color="blue", size=4,shape=18, alpha=0.5) +
      geom_point( aes(x=x,y=y) ) +
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
  
  callModule(module = dyOptIterPlot, id= "optimPlot", 
             paramNames = c(paramNames,"ll"),
             paramListReact = reactive({
               x <- optRes$currentOpt
               if(is.null(x))
                 return(NULL)
               else
                 apply(x$paramIterations,2,function(col){data.frame(idx=seq_along(col),val=col)})
              }),
             group = "ho")
  
  output$restable <- renderDataTable({
    optRes$all %>% filter(geneID == curData$geneID)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

