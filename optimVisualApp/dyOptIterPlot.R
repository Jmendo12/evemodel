# install.packages("dygraphs")
library(dygraphs)

dyOptIterPlotUI <- function(id, paramNames){
  ns <- NS(id)
  tagList(
    lapply(paramNames, function(parName)
      dygraphOutput(ns(paste0("dygraph_",parName)),height = ifelse(rev(paramNames)[1]==parName,"100px","80px"))
    )
  )
}

dyOptIterPlot <- function(input, output, session, paramNames, paramListReact, group) {
  cat("Module call\n")
  # generte output graph for each
  lapply(paramNames, function(parName){
    cat("  generating renderDygraph",parName,"\n")
    output[[paste0("dygraph_",parName)]] <- 
      renderDygraph({
        paramList <- paramListReact()
        if( is.null(paramList) ) {
          cat("renderDygraph",parName," is NULL\n")
          return(NULL)
        }
        cat("renderDygraph",parName,"\n")
        dygraph(paramList[[parName]],ylab = parName, group = group) %>% 
          dyOptions(drawPoints = TRUE, pointSize = 2, stepPlot = TRUE) %>% 
          { if( parName == tail(paramNames,1) )
              dyRangeSelector(., height = 10)
            else
              dyAxis(., "x",axisLabelFontSize = 0) 
          }
      })  
  })
}





##### testing ####
# 
# x <- 
# fitTwoThetasVarFix(isParamFixed = c(F,F,F,F,F),
#                    shiftSpecies = dataSets$salmon20$shiftSpecies,
#                    tree = dataSets$salmon20$tree, 
#                    colSpecies = colnames(dataSets$salmon20$gene.data),
#                    initParams = dataSets$salmon20$initParams[1,],
#                    gene.data.row = dataSets$salmon20$gene.data[1,])
# 
# y <- apply(x$paramIterations,2,function(col){data.frame(idx=seq_along(col),val=col)})
# 
# shinyApp(
#   ui = shinyUI(fluidPage(
#     checkboxInput(inputId = "checkit",label = "activate",value = F),
#     dyOptIterPlotUI(id = "hei",paramNames = names(y))
#   )
#   ),
#   server = shinyServer(function(input, output) {
#     callModule(module = dyOptIterPlot, id= "hei", paramNames = names(y),
#                paramListReact = reactive({if(input$checkit) return(y) else return(NULL)}),group = "ho")
#   })
# )