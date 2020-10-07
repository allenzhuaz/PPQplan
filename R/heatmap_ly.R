#' A General Heatmap for Dynamically Assessing Power of the Sampling Plan Using a General Specification Limit.
#'
#' The function for dynamically plotting (ggplot) the heatmap to evaluate the sampling plan based on a general lower and/or upper specification limits.
#'
#' @usage heatmap_ly(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, test.point, dynamic)
#' @param attr.name (optional) user-defined attribute name for sampling plan assessment
#' @param attr.unit (optional) user-defined attribute unit
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param test.point (optional) actual process data points for testing whether the processes pass PPQ
#' @param dynamic logical; if \code{TRUE}, then convert the plain heatmap to dynamic graph using plotly.
#' @return
#' A Plain or Dynamic Heatmap for Sampling Plan Assessment.
#' @seealso \code{pp} and \code{PPQ.occurve}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' heatmap_ly(attr.name = "Thickness", attr.unit = "%",Llim = -0.2, Ulim = 0.2,
#' mu = seq(-0.2, 0.2, 0.001), sigma = seq(0,0.2, 0.001),
#' test.point=data.frame(c(0.1,-0.05),c(0.15,0.05)), n=2, dynamic = T)
#' }
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @export
heatmap_ly <- function(attr.name="", attr.unit="", Llim, Ulim, mu, sigma, n=1, test.point =c(), dynamic=TRUE){

  para <- expand.grid(mu,sigma)
  ct.df <- data.frame(para, sapply(1:nrow(para), function(i) pp(mu = para[i, 1], sigma = para[i, 2], n = n, Llim = Llim, Ulim = Ulim)))
  colnames(ct.df) <- c("Mean", "Std.Dev", "Pass.Prob")
  p <- ggplot2::ggplot(ct.df, aes(x = Mean, y = Std.Dev, z = Pass.Prob)) +
    theme_bw() +
    theme(panel.background = element_rect(fill = NA), panel.ontop = TRUE) +
    theme(panel.grid.major=element_line(linetype='dashed'),
          panel.grid.minor=element_line(linetype='dashed')) +
    geom_tile(aes(fill=Pass.Prob)) +
    geom_contour(color="white", breaks=c(0.90, 0.95, 0.99)) +
    scale_fill_distiller(palette = "RdYlGn", direction = 1, na.value = "red",
                         limits = c(0.8, 1.0), breaks = c(0.80, 0.90, 0.95, 0.99)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = "Mean Value", y = "Standard Deviation", fill = "Passing\nProbability \n") +
    ggtitle(paste0("Heatmap for ", attr.name, "\nLSL = ", Llim, attr.unit, ", USL = ", Ulim, attr.unit, ", ", n, ifelse(n==1, " Batch", " Batches")))

  if(dynamic==TRUE){
    if(is.null(test.point)){
      ggplotly(p, tooltip = c("x","y", "level", "fill"))
    } else {
      test.point <- data.frame(test.point)
      colnames(test.point) <- c("Mean.Value", "Std.Dev")
      ggplotly(p + geom_point(data = test.point, mapping = aes(x = Mean.Value, y=Std.Dev, z=NULL), shape=8, size=2), tooltip = c("x","y", "level", "fill") )
    }
  } else{
    ct.mat <- matrix(ct.df$Pass.Prob, nrow = length(mu), ncol = length(sigma), byrow = FALSE)
    filled.contour(mu, sigma, ct.mat, levels = c(0, 0.8, 0.9, 0.95, 0.99, 1),
                   col = c("red", "orange", "yellow", "light green", "green"),
                   ylab = paste0("Standard Deviation", " (", attr.unit, ")"),
                   xlab = paste0("Mean for ", attr.name, " (", attr.unit,")"),
                   main = paste0("Heatmap for ", attr.name, "\nLSL = ", Llim, attr.unit, ", USL = ", Ulim, attr.unit, ", for ", n, ifelse(n==1, " Batch", " Batches")),
                   plot.axes = {axis(1)
                     axis(2)
                     points(x = test.point[, 1], y = test.point[, 2], pch = 8)
                   }, key.axes = axis(4, at = c(0.8, 0.9, 0.95, 0.99)))
  }
}
