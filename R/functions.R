utils::globalVariables(c("Mean", "Mean.Value", "Std.Dev", "Pass.Prob"))
#' Probability of Passing Specification Test for a Release Batch
#'
#' The function for calculating the probability of passing critical quality attributes (CQA) specification test .
#'
#' @usage rl.pp(Llim, Ulim, mu, sigma, NV)
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param NV nominal volume for the specification test.
#' @return
#' A numeric value of the passing/acceptance probability
#' @seealso \code{PPQ.pp}, \code{pi.pp} and \code{ti.pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' rl.pp(Llim=1.5, Ulim=3.5, mu=2.5, sigma=0.8)
#' @export
rl.pp <- function(Llim, Ulim, mu, sigma, NV=10){
  if (NV>=10){
  pnorm(Ulim, mean = mu, sd = sigma)-pnorm(Llim, mean = mu, sd = sigma)
  } else if (NV>3 & NV<10){
  (pnorm(Ulim, mean = mu, sd = sigma)-pnorm(Llim, mean = mu, sd = sigma))^3
  } else if (NV>2 & NV<=3){
  (pnorm(Ulim, mean = mu, sd = sigma)-pnorm(Llim, mean = mu, sd = sigma))^5
  } else {
  pnorm(Ulim, mean = mu, sd = sigma/sqrt(5))-pnorm(Llim, mean = mu, sd = sigma/sqrt(5))
  }
}


#' Probability of Passing PPQ Test Using General Multiplier
#'
#' The function for calculating the probability of passing critical quality attributes (CQA) PPQ test .
#'
#' @usage PPQ.pp(Llim, Ulim, mu, sigma, n, n.batch, k)
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param k general mulipler for constructing the specific interval
#' @return
#' A numeric value of the passing/acceptance probability
#' @seealso \code{rl.pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' PPQ.pp(Llim = 90, Ulim = 110, mu=105, sigma=1.5, n=10, k=3.1034)
#'
#' # One-sided tolerance interval with k=0.753 (95/67.5 one-sided tolerance interval LTL)
#' PPQ.pp(sigma=0.03, mu=1.025, n=40, Llim=1, Ulim=Inf, k=0.753)
#'
#' sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN = PPQ.pp, mu=97, n=10, Llim=95, Ulim=105, k=2.373)
#' sapply(X=seq(0.1,10,0.1), FUN = PPQ.pp, mu=97, n=10, Llim=95, Ulim=105, k=2.373)
#'
#' sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN =  PPQ.pp, mu=100, n=10, Llim=95, Ulim=105, k=2.373)
#'
#' sigma <- seq(0.1, 4, 0.1)
#' pp1 <- sapply(X=sigma, FUN =  PPQ.pp, mu=97, n=10, Llim=95, Ulim=105, k=2.373)
#' pp2 <- sapply(X=sigma, FUN =  PPQ.pp, mu=98, n=10, Llim=95, Ulim=105, k=2.373)
#' pp3 <- sapply(X=sigma, FUN =  PPQ.pp, mu=99, n=10, Llim=95, Ulim=105, k=2.373)
#' pp4 <- sapply(X=sigma, FUN =  PPQ.pp, mu=100, n=10, Llim=95, Ulim=105, k=2.373)
#' plot(sigma, pp1, xlab="Standard Deviation", main="LSL=95, USL=105, k=2.373, n=10",
#' ylab="Probability of Passing", type="o", pch=1, col=1, lwd=1, ylim=c(0,1))
#' lines(sigma, pp2, type="o", pch=2, col=2)
#' lines(sigma, pp3, type="o", pch=3, col=3)
#' lines(sigma, pp4, type="o", pch=4, col=4)
#' legend("topright", legend=paste0(rep("mu=",4),c(97,98,99,100)), bg="white",
#' col=c(1,2,3,4), pch=c(1,2,3,4), lty=1, cex=0.8)
#'
#' mu <- seq(95, 105, 0.1)
#' pp5 <- sapply(X=mu, FUN =  PPQ.pp, sigma=0.5, n=10, Llim=95, Ulim=105, k=2.373)
#' pp6 <- sapply(X=mu, FUN =  PPQ.pp, sigma=1, n=10, Llim=95, Ulim=105, k=2.373)
#' pp7 <- sapply(X=mu, FUN =  PPQ.pp, sigma=1.5, n=10, Llim=95, Ulim=105, k=2.373)
#' pp8 <- sapply(X=mu, FUN =  PPQ.pp, sigma=2, n=10, Llim=95, Ulim=105, k=2.373)
#' pp9 <- sapply(X=mu, FUN =  PPQ.pp, sigma=2.5, n=10, Llim=95, Ulim=105, k=2.373)
#' plot(mu, pp5, xlab="Mean Value", main="LSL=95, USL=105, k=2.373, n=10",
#' ylab="Probability of Passing", type="o", pch=1, col=1, lwd=1, ylim=c(0,1))
#' lines(mu, pp6, type="o", pch=2, col=2)
#' lines(mu, pp7, type="o", pch=3, col=3)
#' lines(mu, pp8, type="o", pch=4, col=4)
#' lines(mu, pp9, type="o", pch=5, col=5)
#' legend("topright", legend=paste0(rep("sigma=",5),seq(0.5,2.5,0.5)), bg="white",
#' col=c(1,2,3,4,5), pch=c(1,2,3,4,5), lty=1, cex=0.8)
#' }
#' @export
PPQ.pp <- function(Llim, Ulim, mu, sigma, n=10, n.batch=1, k=2.373){
  Func <- function(V){
    (pnorm(q = Ulim-k*sqrt(V), mean = mu, sd = sigma/sqrt(n))-pnorm(q = Llim+k*sqrt(V), mean = mu, sd = sigma/sqrt(n)))*dchisq(x = (n-1)*V/sigma^2, df = n-1)
  }
  (min(1, integrate(Func, lower = 0, upper = ((Ulim-Llim)/(2*k))^2, rel.tol = 1e-10)$value*(n-1)/sigma^2))^n.batch
}

#' Operating Characteristic (OC) Curves for the CQA PPQ Plan Using General Multiplier.
#'
#' The function for  plotting the OC curve to show the PPQ plan, given lower and upper specification limits.
#'
#' @usage PPQ.occurve(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch, k, add.reference)
#' @param attr.name (optional) user-defined attribute name
#' @param attr.unit (optional) user-defined attribute unit
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param k general mulipler for constructing the specific interval
#' @param add.reference logical; if \code{TRUE}, then add reference OC curves (Baseline and High Performance) in the plot.
#' @return
#' OC curves for specification test and PPQ plan.
#' @seealso \code{PPQ.pp} and \code{rl.pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' PPQ.occurve(attr.name = "Sterile Concentration Assay", attr.unit="%", Llim=95, Ulim=105,
#' mu=97, sigma=seq(0.1, 10, 0.1), n=10, k=2.373, add.reference=TRUE)
#' PPQ.occurve(attr.name = "Sterile Concentration Assay", attr.unit="%", Llim=95, Ulim=105,
#' mu=100, sigma=seq(0.1, 10, 0.1), n=10, k=2.373, add.reference=TRUE)
#' PPQ.occurve(attr.name = "Sterile Concentration Assay", attr.unit="%", Llim=95, Ulim=105,
#' mu=seq(95,105,0.1), sigma=1, n=10, k=2.373)
#' PPQ.occurve(attr.name = "Sterile Concentration Assay", attr.unit="%", Llim=95, Ulim=105,
#' mu=seq(95,105,0.1), sigma=1, n=10, k=2.373, add.reference=TRUE)
#'
#' PPQ.occurve(attr.name = "Protein Concentration", attr.unit="%", Llim=90, Ulim=110,
#' mu=seq(90, 110, 0.1), sigma=1.25, k=2.373)
#' 
#' ## Only display referece curves, leave k as NULL by default
#' PPQ.occurve(attr.name = "Sterile Concentration Assay", attr.unit="%LC", Llim=95, Ulim=105, 
#' mu=98, sigma=seq(0.1, 10, 0.1), n=10, add.reference=TRUE)
#' }
#' @author Yalin Zhu
#' @export
PPQ.occurve <- function(attr.name = "", attr.unit = "", Llim = 1.5, Ulim = 3.5, mu = 2.5, 
                        sigma = 0.1, n = 10, n.batch = 1, k = NULL, add.reference = FALSE) {
  rlpp <- rl.pp(Llim = Llim, Ulim = Ulim, mu = mu, sigma = sigma)
  if (length(mu) == 1 & length(sigma) > 1) {
    if (!is.null(k)) {
      PPQpp <- sapply(X = sigma, FUN = PPQ.pp, mu = mu, n = n, Llim = Llim, 
                      Ulim = Ulim, n.batch = n.batch, k = k)
    }
    plot(sigma, rlpp, xlab = paste0("Standard Deviation (", attr.unit, 
                                    ")"), ylab = "Probability of Passing", type = "l", lty = 1, 
         col = 1, lwd = 2, ylim = c(0, 1))
    if (!is.null(k)) {
      lines(sigma, PPQpp, lty = 1, col = 2, lwd = 2)
    }
    bl.sigma <- sigma[which.min(abs(rlpp - 0.9))]
    hp.sigma <- sigma[which.min(abs(rlpp - 0.95))]
    abline(h = 0.95, lty = 2, col = 16)
    abline(h = 0.9, lty = 2, col = 16)
    abline(h = 0.2, lty = 2, col = 16)
    abline(h = 0.05, lty = 2, col = 16)
    abline(v = bl.sigma, lty = 2, col = 16)
    abline(v = hp.sigma, lty = 2, col = 16)
    points(bl.sigma, 0.2, pch = "*", col = 4, cex = 3)
    points(hp.sigma, 0.05, pch = "*", col = 3, cex = 3)
    if (add.reference == TRUE) {
      bl.k <- optimize(f = function(x) {
        abs(PPQ.pp(sigma = bl.sigma, mu = mu, n = n, Llim = Llim, 
                   Ulim = Ulim, n.batch = n.batch, k = x) - 0.2)
      }, interval = c(0, 10), tol = 1e-06)$minimum
      bl.PPQpp <- sapply(X = sigma, FUN = PPQ.pp, mu = mu, n = n, 
                         Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = bl.k)
      lines(sigma, bl.PPQpp, lty = 2, col = 4, lwd = 2)
      hp.k <- optimize(f = function(x) {
        abs(PPQ.pp(sigma = hp.sigma, mu = mu, n = n, Llim = Llim, 
                   Ulim = Ulim, n.batch = n.batch, k = x) - 0.05)
      }, interval = c(0, 10), tol = 1e-06)$minimum
      hp.PPQpp <- sapply(X = sigma, FUN = PPQ.pp, mu = mu, n = n, 
                         Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = hp.k)
      lines(sigma, hp.PPQpp, lty = 2, col = 3, lwd = 2)
      if (is.null(k)) {
        legend("topright", legend = c(paste0("Spec: ", Llim, "-", 
                                             Ulim, attr.unit), paste0("Interval with k = ", round(bl.k, 
                                                                                                  3), " Baseline"), paste0("Interval with k = ", round(hp.k, 
                                                                                                                                                       3), " High Perform")), bg = "white", col = c(1, 4, 3), 
               lty = c(1, 2, 2), cex = 0.6)
        
      } else {
        legend("topright", legend = c(paste0("Spec: ", Llim, "-", 
                                             Ulim, attr.unit), paste0("Interval with k = ", k), paste0("Interval with k = ", 
                                                                                                       round(bl.k, 3), " Baseline"), paste0("Interval with k = ", 
                                                                                                                                            round(hp.k, 3), " High Perform")), bg = "white", col = c(1, 
                                                                                                                                                                                                     2, 4, 3), lty = c(1, 1, 2, 2), cex = 0.6)
      }
    } else {
      legend("topright", legend = c(paste0("Spec: ", Llim, "-", Ulim, 
                                           attr.unit), paste0("Interval with k = ", k)), bg = "white", 
             col = c("1", "2"), lty = c(1, 1), cex = 0.6)
    }
    title(paste0(attr.name, " OC curves for ", n.batch, " future result \n Assume process mean = ", 
                 mu, attr.unit), cex = 0.9)
  } else if (length(mu) != 1 & length(sigma) == 1) {
    if (!is.null(k)) {
      PPQpp <- sapply(X = mu, FUN = PPQ.pp, sigma = sigma, n = n, 
                      Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = k)
    }
    plot(mu, rlpp, xlab = paste0("Mean Value (", attr.unit, ")"), ylab = "Probability of Passing", 
         type = "l", lty = 1, col = 1, lwd = 2, ylim = c(0, 1))
    if (!is.null(k)) {
      lines(mu, PPQpp, lty = 1, col = 2, lwd = 2)
    }
    bl.mu <- mu[order(abs(rlpp - 0.9))[1:2]]
    hp.mu <- mu[order(abs(rlpp - 0.95))[1:2]]
    abline(h = 0.95, lty = 2, col = 16)
    abline(h = 0.9, lty = 2, col = 16)
    abline(h = 0.2, lty = 2, col = 16)
    abline(h = 0.05, lty = 2, col = 16)
    abline(v = bl.mu, lty = 2, col = 16)
    abline(v = hp.mu, lty = 2, col = 16)
    points(bl.mu[1], 0.2, pch = "*", col = 4, cex = 3)
    points(bl.mu[2], 0.2, pch = "*", col = 4, cex = 3)
    points(hp.mu[1], 0.05, pch = "*", col = 3, cex = 3)
    points(hp.mu[2], 0.05, pch = "*", col = 3, cex = 3)
    if (add.reference == TRUE) {
      bl.k <- optimize(f = function(x) {
        abs(PPQ.pp(mu = bl.mu[1], sigma = sigma, n = n, Llim = Llim, 
                   Ulim = Ulim, n.batch = n.batch, k = x) - 0.2)
      }, interval = c(0, 10), tol = 1e-06)$minimum
      bl.PPQpp <- sapply(X = mu, FUN = PPQ.pp, sigma = sigma, n = n, 
                         Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = bl.k)
      lines(mu, bl.PPQpp, lty = 2, col = 4, lwd = 2)
      hp.k <- optimize(f = function(x) {
        abs(PPQ.pp(mu = hp.mu[1], sigma = sigma, n = n, Llim = Llim, 
                   Ulim = Ulim, n.batch = n.batch, k = x) - 0.05)
      }, interval = c(0, 10), tol = 1e-06)$minimum
      hp.PPQpp <- sapply(X = mu, FUN = PPQ.pp, sigma = sigma, n = n, 
                         Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = hp.k)
      lines(mu, hp.PPQpp, lty = 2, col = 3, lwd = 2)
      if (is.null(k)) {
        legend("center", legend = c(paste0("Spec: ", Llim, "-", 
                                           Ulim, attr.unit), paste0("Interval with k = ", round(bl.k, 
                                                                                                3), " Baseline"), paste0("Interval with k = ", round(hp.k, 
                                                                                                                                                     3), " High Perform")), bg = "white", col = c(1, 4, 3), 
               lty = c(1, 2, 2), cex = 0.6)
        
      } else {
        legend("center", legend = c(paste0("Spec: ", Llim, "-", 
                                           Ulim, attr.unit), paste0("Interval with k = ", k), paste0("Interval with k = ", 
                                                                                                     round(bl.k, 3), " Baseline"), paste0("Interval with k = ", 
                                                                                                                                          round(hp.k, 3), " High Perform")), bg = "white", col = c(1, 
                                                                                                                                                                                                   2, 4, 3), lty = c(1, 1, 2, 2), cex = 0.6)
      }
    } else {
      legend("center", legend = c(paste0("Spec: ", Llim, "-", Ulim, 
                                         attr.unit), paste0("Interval with k = ", k)), bg = "white", 
             col = c("1", "2"), lty = c(1, 1), cex = 0.6)
    }
    title(paste0(attr.name, " OC curves for ", n.batch, " future result \n Assume process standard deviation = ", 
                 sigma, attr.unit), cex = 0.9)
  } else {
    stop("mu and sigma should be one single value and one vector!")
  }
}

#' Heatmap/Contour Plot for Assessing Power of the CQA PPQ Plan Using General Multiplier.
#'
#' The function for plotting the heatmap to evaluate the PPQ plan based on the specification test, given lower and upper specification limits.
#'
#' @usage PPQ.ctplot(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch, k, test.point)
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
#' @return
#' Heatmap (or Countour Plot) for PPQ Assessment.
#' @seealso \code{PPQ.pp} and \code{PPQ.occurve}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' mu <- seq(1.6,3.4,0.05)
#' sigma <- seq(0.05,0.8,0.01)
#' PPQ.ctplot(attr.name = "Total Protein", attr.unit = "mg/mL", Llim=1.5, Ulim=3.5,
#' mu = mu, sigma = sigma, k=2.373)
#'
#' ## Example verifying simulation resutls in the textbook page 249
#' mu <- seq(95, 105, 0.1)
#' sigma <- seq(0.2, 5, 0.1)
#' PPQ.ctplot(attr.name = "Composite Assay", attr.unit = "%LC", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373)
#' mu <- seq(90, 110, 0.5)
#' PPQ.ctplot(attr.name = "Composite Assay", attr.unit = "%LC", Llim=90, Ulim=110,
#' mu = mu, sigma = sigma, k=2.373)
#'
#' mu <- seq(95,105,0.1)
#' sigma <- seq(0.1,2.5,0.1)
#' PPQ.ctplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373)
#' test <- data.frame(mean=c(97,98.3,102.5), sd=c(0.55, 1.5, 1.2))
#' PPQ.ctplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373, test.point=test)
#' }
#' @import graphics
#' @export

PPQ.ctplot <- function(attr.name="", attr.unit="", Llim, Ulim, mu, sigma, n=10, n.batch=1, k, test.point=c()){

  para <- expand.grid(mu,sigma)
  ct.df <- data.frame(para, sapply(1:nrow(para), function(i) PPQ.pp(mu = para[i,1], sigma=para[i,2], n = n, n.batch = n.batch, Llim = Llim, Ulim = Ulim, k = k)))
  colnames(ct.df) <- c("Mean", "Std Dev", "Passing Probability")
  ct.mat<-matrix(ct.df$`Passing Probability`, nrow=length(mu), ncol=length(sigma), byrow=FALSE)

  filled.contour(mu,sigma,ct.mat,levels=c(0, 0.8, 0.9, 0.95, 0.99,1), # ylim=c(0.05,0.27),
                 col=c("red", "orange", "yellow", "light green", "green"), ylab=paste0("Standard Deviation",  " (", attr.unit, ")"),
                 xlab=paste0("Mean for ", attr.name, " (", attr.unit, ")"),
                 main=paste0("Heatmap for ", attr.name, "\nLSL = ", Llim, attr.unit, ", USL = ", Ulim, attr.unit,  ", k = ", k),
                 plot.axes={axis(1); axis(2); points(x=test.point[,1], y=test.point[,2], pch=8)},
                 key.axes=axis(4, at=c(0.80, 0.90, 0.95, 0.99)))
}



#' Heatmap/Contour GGPlot for Dynamically Assessing Power of the CQA PPQ Plan Using General Multiplier.
#'
#' The function for dynamically plotting (ggplot) the heatmap to evaluate the PPQ plan based on the specification test, given lower and upper specification limits.
#'
#' @usage PPQ.ggplot(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch, k,
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
#' @seealso \code{PPQ.pp} and \code{PPQ.occurve}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' mu <- seq(95, 105, 0.1)
#' sigma <- seq(0.1,1.7,0.1)
#' PPQ.ggplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373, dynamic = FALSE)
#' test <- data.frame(mu=c(97,98.3,102.5), sd=c(0.55, 1.5, 0.2))
#' PPQ.ggplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373, test.point = test)
#' }
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @export
PPQ.ggplot <- function(attr.name="", attr.unit="", Llim, Ulim, mu, sigma, n=10, n.batch=1, k, test.point =c(), dynamic=TRUE){

  para <- expand.grid(mu,sigma)
  ct.df <- data.frame(para, sapply(1:nrow(para), function(i) PPQ.pp(mu = para[i,1], sigma=para[i,2], n = n, Llim = Llim, Ulim = Ulim, k=k)))
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

#' Probability of Passing General Upper and/or Lower Specification Limit
#'
#' The function for calculating the probability of passing a general upper and/or lower boundary.
#'
#' @usage pp(Llim, Ulim, mu, sigma, n)
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) 
#' @return
#' A numeric value of the passing/acceptance probability
#' @seealso \code{rl.pp} and \code{PPQ.pp}.
#' @author Yalin Zhu
#' @export
pp <- function (Llim, Ulim, mu, sigma, n=1) {
  (pnorm(Ulim, mean = mu, sd = sigma) - pnorm(Llim, mean = mu, sd = sigma))^n
}


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

