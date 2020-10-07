#' Probability of Passing Specification Test for a Release Batch
#'
#' The function for calculating the probability of passing critical quality attributes (CQA) specification test .
#'
#' @usage rl_pp(Llim, Ulim, mu, sigma, NV)
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param NV nominal volume for the specification test.
#' @return
#' A numeric value of the passing/acceptance probability
#' @seealso \code{PPQ_pp}, \code{pi_pp} and \code{ti_pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' rl_pp(Llim=1.5, Ulim=3.5, mu=2.5, sigma=0.8)
#' @export
rl_pp <- function(Llim, Ulim, mu, sigma, NV=10){
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
