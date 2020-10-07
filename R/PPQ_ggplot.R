#' Heatmap/Contour GGPlot for Dynamically Assessing Power of the CQA PPQ Plan Using General Multiplier.
#'
#' The function for dynamically plotting (ggplot) the heatmap to evaluate the PPQ plan based on the specification test, given lower and upper specification limits.
#'
#' @usage PPQ_ggplot(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch, k,
#' test.point, dynamic)
#' @param attr.name (optional) user-defined attribute name for PPQ assessment
#' @param attr.unit (optional) user-defined attribute unit
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param k general mulipler for constructing the specific interval
#' @param test.point (optional) actual process data points for testing whether the processes pass PPQ
#' @param dynamic logical; if \code{TRUE}, then convert the heatmap ggplot to dynamic graph using plotly.
#' @return
#' Dynamic Heatmap (or Countour Plot) for PPQ Assessment.
#' @seealso \code{PPQ_pp} and \code{PPQ_occurve}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' mu <- seq(95, 105, 0.1)
#' sigma <- seq(0.1,1.7,0.1)
#' PPQ_ggplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373, dynamic = FALSE)
#' test <- data.frame(mu=c(97,98.3,102.5), sd=c(0.55, 1.5, 0.2))
#' PPQ_ggplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373, test.point = test)
#' }
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @export
PPQ_ggplot <- function(attr.name="", attr.unit="", Llim, Ulim, mu, sigma, n=10, n.batch=1, k, test.point =c(), dynamic=TRUE){

  para <- expand.grid(mu,sigma)
  ct.df <- data.frame(para, sapply(1:nrow(para), function(i) PPQ_pp(mu = para[i,1], sigma=para[i,2], n = n, Llim = Llim, Ulim = Ulim, k=k)))
  colnames(ct.df) <- c("Mean.Value", "Std.Dev", "Pass.Prob")
  p <- ggplot2::ggplot(ct.df, aes(x = Mean.Value, y = Std.Dev, z = Pass.Prob)) +
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
    ggtitle(paste0("Heatmap for ", attr.name, "\nLSL = ", Llim, attr.unit, ", USL = ", Ulim, attr.unit,  ", k = ", k))
  if(is.null(test.point)){
    if(dynamic==TRUE){
      ggplotly(p, tooltip = c("x","y", "level", "fill"))
    } else return(p)
  } else{
    test.point <- data.frame(test.point)
    colnames(test.point) <- c("Mean.Value", "Std.Dev")
    if(dynamic==TRUE){
      ggplotly(p + geom_point(data = test.point, mapping = aes(x = Mean.Value, y=Std.Dev, z=NULL), shape=8, size=2), tooltip = c("x","y", "level", "fill") )
    } else {
      return(p + geom_point(data = test.point, mapping = aes(x = Mean.Value, y=Std.Dev, z=NULL), shape=8, size=2) )
    }
  }
}
