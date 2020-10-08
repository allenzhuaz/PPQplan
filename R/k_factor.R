#' Estimating K-factors for Tolerance Intervals Based on Howe's Method
#'
#' Estimates k-factors for tolerance intervals based on Howe's method with normality assumption.
#'
#' @param n Sample size
#' @param alpha The level chosen such that (1-alpha) is the confidence level.
#' @param P The proportion of the population to be covered by the tolerance interval.
#' @param side Whether a 1-sided or 2-sided tolerance interval is required (determined by \code{side = 1} or \code{side = 2}, respectively).
#'
#' @return The estimated k-factor for tolerance intervals assuming normality.
#'
#' @note This function is a simplified version of \code{tolerance::K.factor()}, only considering Howe's method.
#' @seealso \code{ti_pp}
#'
#' @examples
#' k_factor(10, P = 0.95, side = 2)
#'
#' @export
k_factor <- function (n, alpha = 0.05, P = 0.99, side = 1){

  if (side != 1 && side != 2) {
    stop(paste("Must specify a one-sided or two-sided procedure!",
               "\n"))
  }
  if (side == 1) {
    z.p <- qnorm(P)
    ncp <- sqrt(n) * z.p
    t.a <- suppressWarnings(qt(1 - alpha, df = n-1, ncp = ncp))
    K <- t.a/sqrt(n)
  }
  else {
    K.temp <- function(n, alpha, P) {
      f <- n-1
      chi.a <- qchisq(alpha, f)
      k2 <- sqrt(f * qchisq(P, 1, 1/n)/chi.a)
        TEMP4 <- function(n, P, alpha) {
          f <- n-1
          chi.a <- qchisq(alpha, f)
          z.p <- qnorm((1 + P)/2)
          z.a <- qnorm((2 - alpha)/2)
          df.cut <- n^2 * (1 + 1/z.a^2)
          V <- 1 + z.a^2/n + ((3 - z.p^2) * z.a^4)/(6 *
                                                      n^2)
          K.1 <- suppressWarnings(z.p * sqrt(V * (1 +
                                                    (n * V/(2 * f)) * (1 + 1/z.a^2))))
          G <- (f - 2 - chi.a)/(2 * (n + 1)^2)
          K.2 <- suppressWarnings(z.p * sqrt(((f * (1 +
                                                      1/n))/(chi.a)) * (1 + G)))
          if (f > df.cut) {
            K <- K.1
          }
          else {
            K <- K.2
            if (is.na(K))
              K <- 0
          }
          K
        }
        TEMP5 = Vectorize(TEMP4)
        K <- TEMP5(n, P, alpha)

    }
    TEMP <- Vectorize(K.temp)
    K <- TEMP(n = n, alpha = alpha, P = P)
  }
  K
}
