#' Operating Characteristic (OC) Curves for the CQA PPQ Plan Using General Multiplier.
#'
#' The function for  plotting the OC curve to show the PPQ plan, given lower and upper specification limits.
#'
#' @usage PPQ_occurve(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch, k, add.reference)
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
#' @seealso \code{PPQ_pp} and \code{rl_pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' PPQ_occurve(attr.name = "Sterile Concentration Assay", attr.unit="%", Llim=95, Ulim=105,
#' mu=97, sigma=seq(0.1, 10, 0.1), n=10, k=2.373, add.reference=TRUE)
#' PPQ_occurve(attr.name = "Sterile Concentration Assay", attr.unit="%", Llim=95, Ulim=105,
#' mu=100, sigma=seq(0.1, 10, 0.1), n=10, k=2.373, add.reference=TRUE)
#' PPQ_occurve(attr.name = "Sterile Concentration Assay", attr.unit="%", Llim=95, Ulim=105,
#' mu=seq(95,105,0.1), sigma=1, n=10, k=2.373)
#' PPQ_occurve(attr.name = "Sterile Concentration Assay", attr.unit="%", Llim=95, Ulim=105,
#' mu=seq(95,105,0.1), sigma=1, n=10, k=2.373, add.reference=TRUE)
#'
#' PPQ_occurve(attr.name = "Protein Concentration", attr.unit="%", Llim=90, Ulim=110,
#' mu=seq(90, 110, 0.1), sigma=1.25, k=2.373)
#'
#' ## Only display referece curves, leave k as NULL by default
#' PPQ_occurve(attr.name = "Sterile Concentration Assay", attr.unit="%LC", Llim=95, Ulim=105,
#' mu=98, sigma=seq(0.1, 10, 0.1), n=10, add.reference=TRUE)
#' }
#' @author Yalin Zhu
#' @export
PPQ_occurve <- function(attr.name = "", attr.unit = "", Llim = 1.5, Ulim = 3.5, mu = 2.5,
                        sigma = 0.1, n = 10, n.batch = 1, k = NULL, add.reference = FALSE) {
  rlpp <- rl_pp(Llim = Llim, Ulim = Ulim, mu = mu, sigma = sigma)
  if (length(mu) == 1 & length(sigma) > 1) {
    if (!is.null(k)) {
      PPQpp <- sapply(X = sigma, FUN = PPQ_pp, mu = mu, n = n, Llim = Llim,
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
        abs(PPQ_pp(sigma = bl.sigma, mu = mu, n = n, Llim = Llim,
                   Ulim = Ulim, n.batch = n.batch, k = x) - 0.2)
      }, interval = c(0, 10), tol = 1e-06)$minimum
      bl_ppQpp <- sapply(X = sigma, FUN = PPQ_pp, mu = mu, n = n,
                         Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = bl.k)
      lines(sigma, bl_ppQpp, lty = 2, col = 4, lwd = 2)
      hp.k <- optimize(f = function(x) {
        abs(PPQ_pp(sigma = hp.sigma, mu = mu, n = n, Llim = Llim,
                   Ulim = Ulim, n.batch = n.batch, k = x) - 0.05)
      }, interval = c(0, 10), tol = 1e-06)$minimum
      hp_ppQpp <- sapply(X = sigma, FUN = PPQ_pp, mu = mu, n = n,
                         Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = hp.k)
      lines(sigma, hp_ppQpp, lty = 2, col = 3, lwd = 2)
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
      PPQpp <- sapply(X = mu, FUN = PPQ_pp, sigma = sigma, n = n,
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
        abs(PPQ_pp(mu = bl.mu[1], sigma = sigma, n = n, Llim = Llim,
                   Ulim = Ulim, n.batch = n.batch, k = x) - 0.2)
      }, interval = c(0, 10), tol = 1e-06)$minimum
      bl_ppQpp <- sapply(X = mu, FUN = PPQ_pp, sigma = sigma, n = n,
                         Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = bl.k)
      lines(mu, bl_ppQpp, lty = 2, col = 4, lwd = 2)
      hp.k <- optimize(f = function(x) {
        abs(PPQ_pp(mu = hp.mu[1], sigma = sigma, n = n, Llim = Llim,
                   Ulim = Ulim, n.batch = n.batch, k = x) - 0.05)
      }, interval = c(0, 10), tol = 1e-06)$minimum
      hp_ppQpp <- sapply(X = mu, FUN = PPQ_pp, sigma = sigma, n = n,
                         Llim = Llim, Ulim = Ulim, n.batch = n.batch, k = hp.k)
      lines(mu, hp_ppQpp, lty = 2, col = 3, lwd = 2)
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
