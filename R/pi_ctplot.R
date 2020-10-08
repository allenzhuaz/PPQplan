#' Heatmap/Contour Plot for Assessing Power of the CQA PPQ Plan Using Prediction Interval.
#'
#' The function for plotting the heatmap to evaluate the PPQ plan based on the specification test, given lower and upper specification limits.
#'
#' @usage pi_ctplot(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch, alpha, test.point)
#' @param attr.name user-defined attribute name for PPQ assessment
#' @param attr.unit user-defined attribute unit
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param alpha significant level for constructing the prediction interval.
#' @param test.point (optional) actual process data points for testing whether the processes pass PPQ
#' @return
#' Heatmap (or Contour Plot) for PPQ Assessment.
#' @seealso \code{pi_pp} and \code{pi_occurve}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' ## Example verifying simulation resutls in the textbook page 249
#' mu <- seq(95, 105, 0.1)
#' sigma <- seq(0.2, 3.5, 0.1)
#' pi_ctplot(attr.name = "Composite Assay", attr.unit = "%LC",
#' mu = mu, sigma = sigma, Llim=95, Ulim=105)
#' mu <- seq(90, 110, 0.5)
#' pi_ctplot(attr.name = "Composite Assay", attr.unit = "%LC",
#' mu = mu, sigma = sigma, Llim=90, Ulim=110)
#'
#' mu <- seq(95,105,0.1)
#' sigma <- seq(0.1,2.5,0.1)
#' pi_ctplot(attr.name = "Sterile Concentration Assay", attr.unit = "%",
#' mu = mu, sigma = sigma, Llim=95, Ulim=105)
#' test <- data.frame(mean=c(97,98.3,102.5), sd=c(0.55, 1.5, 1.2))
#' pi_ctplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, test.point=test)
#' }
#' @export

pi_ctplot <- function(attr.name="", attr.unit="", Llim=1.5, Ulim=3.5, mu, sigma, n=10, n.batch=1, alpha=0.05, test.point=c()){

  para <- expand.grid(mu,sigma)
  ct.df <- data.frame(para, sapply(1:nrow(para), function(i) pi_pp(mu = para[i,1], sigma=para[i,2], n = n, n.batch = n.batch, Llim = Llim, Ulim = Ulim, alpha = alpha)))
  colnames(ct.df) <- c("Mean", "Std Dev", "Passing Probability")

  ct.mat<-matrix(ct.df$`Passing Probability`, nrow=length(mu), ncol=length(sigma), byrow=FALSE)

  filled.contour(mu,sigma,ct.mat,levels=c(0, 0.8, 0.9, 0.95, 0.99,1),
                 col=c("red", "orange", "yellow", "light green", "green"), ylab=paste0("Standard Deviation",  " (", attr.unit, ")"),
                 xlab=paste0("Mean for ", attr.name, " (", attr.unit, ")"),
                 main=paste0("Probability of Passing ", round(100*(1-alpha),2), "% Confidence Criteria \n", attr.name, " Spec: ", Llim, "-", Ulim, attr.unit),
                 plot.axes={axis(1); axis(2); points(x=test.point[,1], y=test.point[,2], pch=8)},
                 key.axes=axis(4, at=c(0.80, 0.90, 0.95, 0.99)))
}

