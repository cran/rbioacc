#' Plotting method for \code{fitTK} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{fitTK}.  It plots the fit obtained for each
#' variable in the original dataset.
#' 
#' @param x And object returned by fitTK
#' @param \dots Additional arguments
#' 
#' @return a plot of class \code{ggplot}
#' 
#' @export
#' 
#' @import ggplot2
#' 
plot.fitTK <- function(x, ...){
  
  fit <- x

  df <- .df_for_plot(fit)
  
  # HACK TO BE > 0
  df$q50 <- ifelse(df$q50<0,0,df$q50)
  df$qinf95 <- ifelse(df$qinf95<0,0,df$qinf95)
  df$qsup95 <- ifelse(df$qsup95<0,0,df$qsup95)

  plt <- ggplot(data = df) + 
    theme_classic() +
    labs(x = "Time", y = "Concentration") +
    # scale_y_continuous(limits = c(0,NA)) +
    geom_ribbon(
      aes_string(x = 'time', ymin = 'qinf95', ymax = 'qsup95'), fill = "grey80") +
    geom_line(aes_string(x = 'time', y = 'q50'), color = "orange") +
    geom_point(aes_string(x = 'time', y = 'observation' )) + 
    facet_wrap(~variable, scales = "free")

  return(plt)
}

