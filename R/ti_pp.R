#' Probability of Passing PPQ Test using Tolerance Interval
#'
#' The function for calculating the probability of passing critical quality attributes (CQA) PPQ test .
#'
#' @usage ti_pp(Llim, Ulim, mu, sigma, n, n.batch, alpha, coverprob, side)
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param alpha significant level for constructing the tolerance interval
#' @param coverprob coverage probability for constructing the tolerance interval
#' @param side whether a 1-sided or 2-sided tolerance interval is required (determined by side = 1 or side = 2, respectively).
#' @return
#' A numeric value of the passing/acceptance probability
#' @seealso \code{rl_pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' ti_pp(sigma=0.5, mu=2.5, n=10, n.batch=1, Llim=1.5, Ulim=3.5, alpha=0.05)
#'
#' sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN = ti_pp, mu=97, n=10, Llim=95, Ulim=105,
#' n.batch=1, alpha=0.05)
#' sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN = ti_pp, mu=100, n=10, Llim=95, Ulim=105,
#' n.batch=1, alpha=0.05)
#' @import stats
#' @export
ti_pp <- function(Llim, Ulim, mu, sigma, n=10, n.batch=1, alpha=0.05, coverprob = 0.675, side = 2){
  k <- k_factor(n = n, alpha = alpha, P = coverprob, side = side)
  Func <- function(V){
    (pnorm(q = Ulim-k*sqrt(V), mean = mu, sd = sigma/sqrt(n))-pnorm(q = Llim+k*sqrt(V), mean = mu, sd = sigma/sqrt(n)))*dchisq(x = (n-1)*V/sigma^2, df = n-1)
  }
  (min(1, integrate(Func, lower = 0, upper = ((Ulim-Llim)/(2*k))^2, rel.tol = 1e-10)$value*(n-1)/sigma^2))^n.batch
}
