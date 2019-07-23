#' Plot positive and negative predictive values across different prevalences
#' 
#' @description Positive and negative predictive values (PPV and NPV) are two diagnostic statistics that change depending on the prevalence, so if you don't have a discrete prevalence to work with
#' this function can create a plot that shows the positive and negative predictive values across all possible prevalences (as long as you have already calculated the sensitivity and specificity).
#' @usage predvalPlot(sens, spec, nsteps=1000, title=NULL, rounding=2)
#' 
#' @param sens the sensitivity of the prediction
#' @param spec the specificity of the prediction
#' @param nsteps the number of steps between prevalence 0\% and 100\% (i.e. the number of steps in the X-axis) (default: 1000)
#' @param title title of the plot (if left blank, it will just indicate the input sensitivity and specificity)
#' @param rounding number of significant digits for displaying the sensitivity, specificity, PPV, and NPV (default: 2)
#' @return Plotly plot of predictive values vs. prevalence
#' @author Lara Murphy, Aditya M. Rao
#' @examples
#' predvalPlot(sens = 0.9, spec = 0.8)
#' @importFrom magrittr %>%
#' @import httpuv
#' @export
predvalPlot <- function(sens, spec, nsteps=1000, title=NULL, rounding=2) {
  if(is.null(sens) || is.null(spec)){stop("Sensitivity and specificity must be provided.")}
  i = 1:nsteps
  k = i/nsteps
  ppv = ((k*sens)/((k*sens) + (1-spec)*(1-k)))*100
  npv = ((spec * (1-k)) / (((1-sens)*k) + (spec*(1-k))))*100
  dt = data.table::as.data.table(k, keep.rownames = F)
  dt$ppv = ppv
  dt$npv = npv
  setnames(dt, "k", "prevalence")
  df = as.data.frame(dt)
  #df$text = paste0("Prevalence: ",df$prevalence*100,"%\nPredictive Value: ",round(sens*100,rounding),"%")
  if(is.null(title)){
    title = paste0("Positive and Negative Predictive Value\n",round(sens*100,rounding),"% Sensitivity and ",round(spec*100,rounding),"% Specificity")
  }
  
  #NOTE: eventually someone should make a ggplot version of this too
  
  df %>% plotly::plot_ly(.,type="scatter",mode="markers") %>%
    plotly::add_trace(y = ~ppv,
                      x = ~prevalence*100,
                      name = "PPV (%)",
                      hoverinfo="text",
                      text = paste0("Prevalence: ",df$prevalence*100,"%\nPPV: ",round(df$ppv,rounding),"%"),
                      marker = list(
                        color = 'rgb(71, 150, 215)',
                        opacity = 0.8,
                        size = 8)) %>%
    plotly::add_trace(y = ~npv,
                      x = ~prevalence*100,
                      name = "NPV (%)",
                      hoverinfo="text",
                      text = paste0("Prevalence: ",df$prevalence*100,"%\nNPV: ",round(df$npv,rounding),"%"),
                      marker = list(
                        color = 'rgb(249, 111, 88)',
                        opacity = 0.8,
                        size = 8)) %>%
    plotly::layout(title=title,
                   xaxis = list(title = "Prevalence (%)"),
                   yaxis = list(title = "Predictive Value (%)"))
}

#Note: I think this is the only function that depends on httpuv
