#' Probability of Passing PPQ Test Using General Multiplier
#'
#' The function for calculating the probability of passing critical quality attributes (CQA) PPQ test .
#'
#' @usage PPQ_pp(Llim, Ulim, mu, sigma, n, n.batch, k)
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param k general multiplier for constructing the specific interval
#' @return
#' A numeric value of the passing/acceptance probability
#' @seealso \code{rl_pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' PPQ_pp(Llim = 90, Ulim = 110, mu=105, sigma=1.5, n=10, k=3.1034)
#'
#' # One-sided tolerance interval with k=0.753 (95/67.5 one-sided tolerance interval LTL)
#' PPQ_pp(sigma=0.03, mu=1.025, n=40, Llim=1, Ulim=Inf, k=0.753)
#'
#' sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN = PPQ_pp, mu=97, n=10, Llim=95, Ulim=105, k=2.373)
#' sapply(X=seq(0.1,10,0.1), FUN = PPQ_pp, mu=97, n=10, Llim=95, Ulim=105, k=2.373)
#'
#' sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN =  PPQ_pp, mu=100, n=10, Llim=95, Ulim=105, k=2.373)
#'
#' sigma <- seq(0.1, 4, 0.1)
#' pp1 <- sapply(X=sigma, FUN =  PPQ_pp, mu=97, n=10, Llim=95, Ulim=105, k=2.373)
#' pp2 <- sapply(X=sigma, FUN =  PPQ_pp, mu=98, n=10, Llim=95, Ulim=105, k=2.373)
#' pp3 <- sapply(X=sigma, FUN =  PPQ_pp, mu=99, n=10, Llim=95, Ulim=105, k=2.373)
#' pp4 <- sapply(X=sigma, FUN =  PPQ_pp, mu=100, n=10, Llim=95, Ulim=105, k=2.373)
#' plot(sigma, pp1, xlab="Standard Deviation", main="LSL=95, USL=105, k=2.373, n=10",
#' ylab="Probability of Passing", type="o", pch=1, col=1, lwd=1, ylim=c(0,1))
#' lines(sigma, pp2, type="o", pch=2, col=2)
#' lines(sigma, pp3, type="o", pch=3, col=3)
#' lines(sigma, pp4, type="o", pch=4, col=4)
#' legend("topright", legend=paste0(rep("mu=",4),c(97,98,99,100)), bg="white",
#' col=c(1,2,3,4), pch=c(1,2,3,4), lty=1, cex=0.8)
#'
#' mu <- seq(95, 105, 0.1)
#' pp5 <- sapply(X=mu, FUN =  PPQ_pp, sigma=0.5, n=10, Llim=95, Ulim=105, k=2.373)
#' pp6 <- sapply(X=mu, FUN =  PPQ_pp, sigma=1, n=10, Llim=95, Ulim=105, k=2.373)
#' pp7 <- sapply(X=mu, FUN =  PPQ_pp, sigma=1.5, n=10, Llim=95, Ulim=105, k=2.373)
#' pp8 <- sapply(X=mu, FUN =  PPQ_pp, sigma=2, n=10, Llim=95, Ulim=105, k=2.373)
#' pp9 <- sapply(X=mu, FUN =  PPQ_pp, sigma=2.5, n=10, Llim=95, Ulim=105, k=2.373)
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
PPQ_pp <- function(Llim, Ulim, mu, sigma, n=10, n.batch=1, k=2.373){
  Func <- function(V){
    (pnorm(q = Ulim-k*sqrt(V), mean = mu, sd = sigma/sqrt(n))-pnorm(q = Llim+k*sqrt(V), mean = mu, sd = sigma/sqrt(n)))*dchisq(x = (n-1)*V/sigma^2, df = n-1)
  }
  (min(1, integrate(Func, lower = 0, upper = ((Ulim-Llim)/(2*k))^2, rel.tol = 1e-10)$value*(n-1)/sigma^2))^n.batch
}
