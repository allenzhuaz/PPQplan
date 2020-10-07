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
#' @seealso \code{rl_pp} and \code{PPQ_pp}.
#' @author Yalin Zhu
#' @export
pp <- function (Llim, Ulim, mu, sigma, n=1) {
  (pnorm(Ulim, mean = mu, sd = sigma) - pnorm(Llim, mean = mu, sd = sigma))^n
}
