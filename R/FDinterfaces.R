#'@import stats
#'@import graphics
#'@import utils
#'@import grDevices

# ---------------------------- FD.density ----------------------------

#' @title The Flexible Dirichlet Density Function
#'
#' @description Density function on the simplex for the Flexible Dirichlet distribution with parameters \code{a}, \code{p} and \code{t}.
#'
#' @param x vector of a point on the simplex. It must sum to one.
#' @param a vector of the non-negative alpha parameters.
#' @param p vector of the clusters' probabilities. It must sum to one.
#' @param t non-negative scalar tau parameter.
#'
#' @details Vectors \code{x}, \code{a} and \code{p} must be of the same length.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' x <- c(0.1,0.25,0.65)
#' alpha <- c(12,7,15)
#' prob <- c(0.3,0.4,0.3)
#' tau <- 8
#' FD.density(x,alpha,prob,tau)
#'
#' @seealso \code{\link{FD.theorcontours}}, \code{\link{FD.generate}}
#'
#' @export

FD.density <- function(x, a, p, t) {
  if (!is.matrix(x)) 
    if (is.data.frame(x)) 
      x <- as.matrix(x) else x <- t(x)
      if (any(round(rowSums(x), 7) != 1)) 
        stop("Values of 'x' must sum to 1")
      if (!is.matrix(a)) 
        a <- matrix(a, ncol = length(a), nrow = nrow(x), byrow = TRUE)
      if (any(dim(x) != dim(a))) 
        stop("Mismatch between dimensions of 'x' and 'alpha'")
      if (!is.matrix(p)) 
        p <- matrix(p, ncol = length(p), nrow = nrow(x), byrow = TRUE)
      if (any(dim(x) != dim(p))) 
        stop("Mismatch between dimensions of 'x' and 'p'")
      if (!is.matrix(t)) 
        t <- matrix(t, ncol = length(t), nrow = nrow(x), byrow = TRUE)
      pd <- vector(length = nrow(x))
      for (i in 1:nrow(x)) pd[i] <- dFDir(x[i, ], a[i, ], p[i, ], t[i, ])
      pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
      pd[apply(x, 1, function(z) all.equal(round(sum(z), 7), 1) != TRUE)] <- 0
      pd
}

# ---------------------------- FD.generate ----------------------------

#' @title The Flexible Dirichlet Random Generation
#'
#' @description Random generation from the Flexible Dirichlet distribution with parameters \code{a}, \code{p} and \code{t}.
#'
#' @param n number of points on the simplex to be generated.
#' @param a vector of the non-negative alpha parameters.
#' @param p vector of the clusters' probabilities. It must sum to one.
#' @param t non-negative scalar tau parameter.
#'
#' @details Vectors \code{a} and \code{p} must be of the same length.
#' The Flexible Dirichlet distribution derives from the normalization of a basis of positive dependent random variables obtained by starting from a basis of independent equally scaled gamma random variables, and randomly allocating to the \code{i}-th element a further independent gamma random variable.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, 1--21.
#' }
#'
#' @examples
#' n <- 100
#' alpha <- c(12,7,15)
#' prob <- c(0.3,0.4,0.3)
#' tau <- 8
#' data <- FD.generate(n,alpha,prob,tau)
#' data
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.density}}, \code{\link{FD.theorcontours}}, \code{\link{FD.subcomposition}}, \code{\link{FD.amalgamation}}
#'
#' @export

FD.generate <- function(n, a, p, t) {
  if (length(a) != length(p)) 
    stop("Error: Alpha and P do not have the same length") else if (length(t) != 1) 
      stop("Error: Tau must be a single number") else if (round(sum(p), 5) != 1) 
        stop("Error: P must sum to one") else if (any(a <= 0)) 
          stop("Error: Alpha parameters must be non-negative") else if (t <= 0) 
            stop("Error: Tau parameter must be non-negative") else rFDir(n = n, alpha = a, p = p, tau = t)
}

# ---------------------------- FD.estimation ----------------------------

#' @title Flexible Dirichlet Estimation
#'
#' @description Estimates the vector of parameters of a Flexible Dirichlet distribution through an EM-based maximum likelihood approach.
#'
#' @param data a matrix or a dataframe containing only the variables in the model. Rows must sum to one, or \code{normalize} must be set \code{TRUE}.
#' @param normalize if \code{TRUE}, each row of \code{data} will be divided by its own total to become a point of the simplex. Values in \code{data} must be positive.
#' @param iter.initial.SEM number of iterations for the initial SEM step. Default to 50.
#' @param iter.final.EM number of iterations for the final EM step. Default to 100.
#' @param verbose if \code{TRUE}, the progression of the elaboration and the results will be printed on screen.
#'
#' @return an object of class FDfitted. It's a list composed by:
#' \describe{
#'   \item{\code{alpha}}{Estimated values of the parameter vector Alpha}
#'   \item{\code{p}}{Estimated values of the parameter vector P}
#'   \item{\code{tau}}{Estimated value of the parameter Tau}
#'   \item{\code{logL}}{LogLikelihood}
#'   \item{\code{data}}{Normalized dataset}
#' }
#'
#' @details The procedure is made up of four stages:
#' \enumerate{
#'   \item Clustering: The algorithm applies many different clustering rules to the dataset, in order to exploit the specific cluster patterns that the parameter structure of the model involves.
#'   \item Labelling: Once the initial partitions are obtained, group labeling needs to be established because any clustering algorithm assigns the group labels randomly, but the FD cluster structure entails a precise labelling scheme.
#'   \item Initial SEM: A Stochastic E-M algorithm is applied at every initial partition and every possible label permutation identified.
#'   \item Final E-M: The previous step must be seen as a multiple initialization strategy. At this point only the best one is selected and a final E-M algorithm is used to find the point that maximizes the likelihood of the parameter vector.
#' }
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' summary(results)
#'
#' @seealso \code{\link{FD.generate}}, \code{\link{FD.stddev}}, \code{\link{FD.aicbic}}, \code{\link{FD.barycenters}}, \code{\link{FD.ternaryplot}}, \code{\link{FD.rightplot}}, \code{\link{FD.marginalplot}}
#'
#' @export

FD.estimation <- function(data, normalize = F, iter.initial.SEM = 50, iter.final.EM = 100, 
                          verbose = T) {
  dataset <- as.matrix(data)
  D <- dim(dataset)[2]
  # n <- dim(dataset)[1] proviamo
  colnames(dataset) <- c(paste("VAR", 1:D, sep = ""))
  if (!is.null(colnames(data))) 
    for (i in 1:D) if (colnames(data)[i] != "") 
      colnames(dataset)[i] <- colnames(data)[i]
  if (D < 2) 
    stop("Error: At least 2 variables are needed.\n") else if (!normalize && (sum(round(rowSums(dataset), 7) == 1) < (dim(dataset)[1]))) {
      print(head(dataset[(round(rowSums(dataset), 7) != 1), ]))
      stop("Error: Printed rows must sum to one. Set normalize=TRUE.\n")
    } else {
      if (normalize) 
        dataset <- FD.normalization(dataset)
      numD <- dim(dataset)[2]
      stima <- fd.estimation.short(dataset, SEM.iter = iter.initial.SEM, FIN.iter = iter.final.EM, 
                                   verbose)
      output <- list(alpha = stima[(numD + 2):(2 * numD + 1)], p = stima[(2):(numD + 
                                                                                1)], tau = stima[2 * numD + 2], logL = stima[1], data = dataset)
      if (verbose) {
        cat(c("Alpha:\n"))
        cat(c(round(output$alpha, 2), "\n"))
        cat(c("P:\n"))
        cat(c(round(output$p, 4), "\n"))
        cat(c("Tau:\n"))
        cat(c(round(output$tau, 2), "\n"))
        cat(c("LogLikelihood:\n"))
        cat(c(round(output$logL, 7), "\n"))
      }
      class(output) <- "FDfitted"
      invisible(output)
    }
}

# ---------------------------- print.FDfitted ----------------------------
#' @title Print Method for FDfitted Objects
#'
#' @description This method shows the results of \code{\link{FD.estimation}}.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#' @param ... additional arguments
#'
#' @export

print.FDfitted <- function(x, ...) {
  cat(c("Alpha:\n"))
  cat(c(round(x$alpha, 2), "\n"))
  cat(c("P:\n"))
  cat(c(round(x$p, 4), "\n"))
  cat(c("Tau:\n"))
  cat(c(round(x$tau, 2), "\n"))
  cat(c("LogLikelihood:\n"))
  cat(c(round(x$logL, 7), "\n"))
  invisible(x)
}

# ---------------------------- FD.stddev ----------------------------

#' @title Standard Deviation of the ML estimators of a Flexible Dirichlet
#'
#' @description Conditional Bootstrap evaluation of the standard errors of the maximum likelihood parameter estimates of a Flexible Dirichlet distribution.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#' @param iter.bootstrap number of iterations of the Bootstrap.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' FD.stddev(results)
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.aicbic}}, \code{\link{FD.barycenters}}
#'
#' @export

FD.stddev <- function(x, iter.bootstrap = 500) {
  data <- x$data
  alpha <- x$alpha
  p <- x$p
  tau <- x$tau
  D <- length(alpha)
  
  ris <- boot.FD(alpha = alpha, p = p, tau = tau, data = data, B = iter.bootstrap)
  colnames(ris) <- c(paste("alpha", 1:D, sep = ""), paste("p", 1:D, sep = ""), 
                     "tau")
  rownames(ris) <- "Std Dev"
  ris
}

# ---------------------------- FD.aicbic ----------------------------
#' @title Information Criterions of a Flexible Dirichlet Model
#'
#' @description
#' Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC) of a fitted Flexible Dirichlet model.
#' An Information Criterion for one fitted model object for which a log-likelihood value can be obtained is defined as
#' \eqn{ -2*log-likelihood + k*npar}, where \eqn{ npar} represents the number of parameters in the fitted model, and \eqn{ k = 2} for AIC, or \eqn{ k = log(n)} for BIC (\eqn{ n} being the number of observations).
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' FD.aicbic(results)
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.stddev}}, \code{\link{FD.barycenters}}
#'
#' @export

FD.aicbic <- function(x) {
  output <- AIC.Fdir(x)
  output
}

# ---------------------------- FD.ternaryplot ----------------------------

#' @title Ternary Plot of a Flexible Dirichlet
#'
#' @description Ternary plot and contour lines of the density function of a fitted Flexible Dirichlet distribution.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#' @param zoomed if \code{TRUE}, shows only the area where most of the density is concentrated. If \code{FALSE}, shows the whole area of the ternary diagram.
#' @param showgrid if \code{TRUE}, shows the axis and the labels. If \code{FALSE}, only the graph is printed.
#' @param showdata if \code{TRUE}, prints the data points. If \code{FALSE}, shows only the contour lines.
#' @param nlevels approximate number of contour lines to be drawn.
#'
#' @details The number of variables in the fitted model must be 3 to draw a ternary plot.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' FD.ternaryplot(results)
#' FD.ternaryplot(results, zoomed=FALSE, showgrid=TRUE, showdata=FALSE, nlevels=3)
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.rightplot}}, \code{\link{FD.marginalplot}}
#'
#' @export

FD.ternaryplot <- function(x, zoomed = T, showgrid = T, showdata = T, nlevels = 10) {
  if (length(x$alpha) != 3) 
    stop("Number of variables must be 3 to draw a ternary plot.")
  if (zoomed) 
    zoom_ternary_plot_FD(x, showgrid = showgrid, showdata = showdata, nlevels = nlevels) else ternary_plot_FD(x, showgrid = showgrid, showdata = showdata, nlevels = nlevels)
}


# ---------------------------- FD.rightplot ----------------------------
#' @title Right Triangle Plot of a Flexible Dirichlet
#'
#' @description Right triangle plot and contour lines of the density function of a fitted Flexible Dirichlet distribution.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#' @param var numeric vector containing the two variables to be plotted on the axis.
#' @param zoomed if \code{TRUE}, shows only the area where most of the density is concentrated. If \code{FALSE}, shows the whole area of the right triangle.
#' @param showgrid if \code{TRUE}, shows the axis and the labels. If \code{FALSE}, only the graph is printed.
#' @param showdata if \code{TRUE}, prints the data points. If \code{FALSE}, shows only the contour lines.
#' @param nlevels approximate number of contour lines to be drawn.
#'
#' @details The number of variables in the fitted model must be 3 to draw a plot on the right triangle.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' FD.rightplot(results)
#' FD.rightplot(results, var=c(3,2), zoomed=FALSE, showgrid=TRUE, showdata=FALSE, nlevels=3)
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.ternaryplot}}, \code{\link{FD.marginalplot}}
#'
#' @export

FD.rightplot <- function(x, var = c(1, 2), zoomed = T, showgrid = T, showdata = T, 
                         nlevels = 10) {
  if (length(x$alpha) != 3) 
    stop("Number of variables must be 3 to draw a right triangle plot.")
  if (!(length(var) == 2 && var[1] %in% (1:3) && var[2] %in% (1:3))) 
    stop("Variable var must contain a vector of 2 numbers out of 1, 2 or 3.")
  if (zoomed) 
    zoom_right_triangle_plot_FD(x, var, showgrid = showgrid, showdata = showdata, 
                                nlevels = nlevels) else right_triangle_plot_FD(x, var, showgrid = showgrid, showdata = showdata, 
                                                                               nlevels = nlevels)
}

# ---------------------------- FD.marginalplot ----------------------------

#' @title Marginal Plot of a Flexible Dirichlet
#'
#' @description Histogram of the observed marginal variable and estimated density function of the marginal variable of a fitted Flexible Dirichlet distribution.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#' @param var position of the variable to be plotted.
#' @param zoomed if \code{TRUE}, shows only the area where most of the density is concentrated. If \code{FALSE}, shows the whole range \code{[0;1]}.
#' @param showgrid if \code{TRUE}, shows the axis and the labels. If \code{FALSE}, only the graph is printed.
#' @param showdata if \code{TRUE}, prints the histogram of the data. If \code{FALSE}, shows only the density function.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' FD.marginalplot(results, var=2)
#' FD.marginalplot(results, var=2, zoomed=FALSE, showgrid=TRUE, showdata=FALSE)
#'
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.ternaryplot}}, \code{\link{FD.rightplot}}
#'
#' @export

FD.marginalplot <- function(x, var, zoomed = T, showgrid = T, showdata = T) {
  D <- length(x$alpha)
  if (!(length(var) == 1) || !(var %in% 1:D)) 
    stop("Variable var is not a valid number") else marginal_plot_FD(x, var, zoomed = zoomed, showdata = showdata, showgrid = showgrid)
}

# ---------------------------- FD.FDfitted ----------------------------

#' @title Plot Method for FDfitted Objects
#'
#' @description This method plots the results of \code{\link{FD.estimation}}, using the functions \code{\link{FD.ternaryplot}} or \code{\link{FD.rightplot}}.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#' @param type string containing \code{'ternary'} or \code{'right'}.
#' @param var numeric vector containing the two variables to be plotted on the axis. Used only if \code{type='right'}.
#' @param zoomed if \code{TRUE}, shows only the area where most of the density is concentrated. If \code{FALSE}, shows the whole area.
#' @param showgrid if \code{TRUE}, shows the axis and the labels. If \code{FALSE}, only the graph is printed.
#' @param showdata if \code{TRUE}, prints the data points. If \code{FALSE}, shows only the contour lines.
#' @param nlevels approximate number of contour lines to be drawn.
#' @param ... additional arguments
#'
#' @details The number of variables in the fitted model must be 3 to draw a plot.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' plot(results)
#' plot(results, type='right', var=c(3,2), zoomed=FALSE, showgrid=TRUE, showdata=FALSE, nlevels=3)
#'
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.ternaryplot}}, \code{\link{FD.rightplot}}, \code{\link{FD.marginalplot}}
#'
#' @export

plot.FDfitted <- function(x, type = "ternary", var = c(1, 2), zoomed = T, showgrid = T, 
                          showdata = T, nlevels = 10, ...) {
  D <- dim(x$data)[2]
  if (D != 3) 
    stop("Number of variables in the model must be 3 to draw a plot.")
  if (type != "ternary" && type != "right") 
    stop("Argument 'type' must be 'ternary' or 'right'.")
  if (type == "ternary") 
    FD.ternaryplot(x, zoomed = zoomed, showgrid = showgrid, showdata = showdata, 
                   nlevels = nlevels) else FD.rightplot(x, var, zoomed = zoomed, showgrid = showgrid, showdata = showdata, 
                                                        nlevels = nlevels)
}


# ---------------------------- FD.theorcontours ----------------------------

#' @title Contour Lines of a Flexible Dirichlet
#'
#' @description Contour lines of a Flexible Dirichlet with given parameters on the ternary diagram or on the right triangle.
#'
#' @param a vector of the non-negative alpha parameters.
#' @param p vector of the clusters' probabilities. It must sum to one.
#' @param t non-negative scalar tau parameter.
#' @param type string indicating whether to plot the contour lines on a ternary diagram \code{'ternary'}, or on a right triangle plot \code{'right'}.
#' @param var numeric vector containing the two variables to be plotted on the axis. Used only if \code{type='right'}.
#' @param zoomed if \code{TRUE}, shows only the area where most of the density is concentrated. If \code{FALSE}, shows the whole area.
#' @param showgrid if \code{TRUE}, shows the axis and the labels. If \code{FALSE}, only the graph is printed.
#' @param nlevels approximate number of contour lines to be drawn.
#'
#' @details The number of variables in the Flexible Dirichlet must be 3 to draw a plot. Vectors \code{a} and \code{p} must be of the same length.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' alpha <- c(12,7,15)
#' prob <- c(0.3,0.4,0.3)
#' tau <- 8
#' FD.theorcontours(alpha,prob,tau)
#' FD.theorcontours(alpha,prob,tau, type='right', var=c(3,2), zoomed=FALSE, showgrid=TRUE, nlevels=3)
#'
#' @seealso \code{\link{FD.generate}}, \code{\link{FD.density}}
#'
#' @export

FD.theorcontours <- function(a, p, t, type = "ternary", var = c(1, 2), zoomed = T, 
                             showgrid = T, nlevels = 10) {
  if (length(a) != length(p)) 
    stop("Error: Alpha and P do not have the same length")
  if (length(t) != 1) 
    stop("Error: Tau must be a single number")
  if (round(sum(p), 5) != 1) 
    stop("Error: P must sum to one")
  D <- length(a)
  if (D != 3) 
    stop("Number of variables in the model must be 3 to draw contour lines.")
  if (type != "ternary" && type != "right") 
    stop("Argument 'type' must be 'ternary' or 'right'.")
  
  datasim <- rFDir(n = 1000, alpha = a, p = p, tau = t)
  colnames(datasim) <- c(paste("VAR", 1:D, sep = ""))
  x <- list(alpha = a, p = p, tau = t, logL = NA, data = datasim)
  if (type == "ternary") 
    FD.ternaryplot(x, zoomed = zoomed, showgrid = showgrid, showdata = F, 
                   nlevels = nlevels) else FD.rightplot(x, var, zoomed = zoomed, showgrid = showgrid, showdata = F, 
                                                        nlevels = nlevels)
}


# ---------------------------- FD.barycenters ----------------------------

#' @title Cluster Barycenters of a Flexible Dirichlet model
#'
#' @description Cluster barycenters of a fitted Flexible Dirichlet distribution.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' FD.barycenters(results)
#'
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.clusterdistances}}, \code{\link{FD.moments}}
#'
#' @export

FD.barycenters <- function(x) {
  alpha <- x$alpha
  tau <- x$tau
  D <- length(alpha)
  bar <- mu.FD(alpha, tau)
  labs <- colnames(x$data)
  colnames(bar) <- labs
  rownames(bar) <- paste("cluster", 1:D, sep = "")
  bar
}

# ---------------------------- FD.clusterdistances ----------------------------

#' @title Flexible Dirichlet Cluster Distances
#'
#' @description Returns a measure of symmetrized Kullback-Leibler distance between mixture component densities of a fitted Flexible Dirichlet distribution.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' FD.clusterdistances(results)
#' 
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.barycenters}}, \code{\link{FD.moments}}
#'
#' @export

FD.clusterdistances <- function(x) {
  alpha <- x$alpha
  tau <- x$tau
  bar <- t(Kappa(tau, alpha))
  # labs <- colnames(x$data)
  rownames(bar) <- "Distances"
  bar
}

# ---------------------------- FD.moments ----------------------------

#' @title Flexible Dirichlet Moments
#'
#' @description Moments of a fitted Flexible Dirichlet distribution. The function returns the mean and variance vectors and the covariance and correlation matrices.
#'
#' @param x an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' FD.moments(results)
#' 
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.barycenters}}, \code{\link{FD.clusterdistances}}
#'
#' @export

FD.moments <- function(x) {
  alpha <- x$alpha
  p <- x$p
  tau <- x$tau
  labs <- colnames(x$data)
  
  mean <- t(E.Fdir(alpha, p, tau))
  colnames(mean) <- labs
  rownames(mean) <- ""
  var <- t(V.Fdir(alpha, p, tau))
  colnames(var) <- labs
  rownames(var) <- ""
  cov <- Cov.Fdir(alpha, p, tau)
  colnames(cov) <- labs
  rownames(cov) <- labs
  cor <- Cor.Fdir(alpha, p, tau)
  colnames(cor) <- labs
  rownames(cor) <- labs
  output <- list(mean = mean, var = var, cov = cov, cor = cor)
  output
}

# ---------------------------- summary.FDfitted ----------------------------

#' @title Summary Method for FDfitted Objects
#'
#' @description This method summarizes the results of \code{\link{FD.estimation}}, adding also information from the functions \code{\link{FD.stddev}} and \code{\link{FD.aicbic}}.
#'
#' @param object an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.
#' @param ... additional arguments
#'
#' @return A list composed by:
#' \describe{
#'   \item{\code{par}}{Estimated parameter vector}
#'   \item{\code{sd}}{Vector of the estimated standard deviations}
#'   \item{\code{goodness}}{Vector containing LogLikelihood, AIC and BIC}
#' }
#'
#' @references {
#'  Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
#' \cr
#'  Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
#' }
#'
#' @examples
#' data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
#' data
#' results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
#' results
#' summary(results)
#' 
#'
#' @seealso \code{\link{FD.estimation}}, \code{\link{FD.stddev}}, \code{\link{FD.aicbic}}
#'
#' @export

summary.FDfitted <- function(object, ...) {
  x <- object
  D <- length(x$alpha)
  par <- t(c(x$alpha, x$p, x$tau))
  colnames(par) <- c(paste("alpha", 1:D, sep = ""), paste("p", 1:D, sep = ""), 
                     "tau")
  rownames(par) <- "values"
  sd <- FD.stddev(x)
  aicbic <- FD.aicbic(x)
  goodness <- t(c(x$logL, aicbic$AIC, aicbic$BIC))
  colnames(goodness) <- c("logLik", "AIC", "BIC")
  rownames(goodness) <- ""
  # m <- FD.moments(x) barycenters <- FD.barycenters(x) clusterdist <-
  # FD.clusterdistances(x) output <- list(par=par,sd=sd,goodness=goodness,
  # mean=m$mean,var=m$var,cov=m$cov,cor=m$cor,
  # barycenters=barycenters,clusterdist=clusterdist)
  output <- list(par = par, sd = sd, goodness = goodness)
  output
}

# ---------------------------- FD.subcomposition ----------------------------

#' @title Subcomposition
#'
#' @description Given a matrix or a numeric dataframe, this function returns a subcomposition made up of the specified columns.
#'
#' @param data a matrix or a dataframe containing only variables in the model.
#' @param columns numeric vector containing the position of the columns to keep in the new composition.
#'
#' @details Values must be positive. In case one row-entry (or more) is NA, the whole row will be returned as NA.

#' @examples
#' data(oliveoil)
#' dataoil <- oliveoil
#' head(dataoil)
#' data <- FD.normalization(dataoil[,3:10])
#' head(data)
#' data.sub <- FD.subcomposition(data,c(1,3,4,5))
#' head(data.sub)
#' data.amalg <- FD.amalgamation(data,c(2,6,7,8),name='others')
#' head(data.amalg)
#'
#' @seealso \code{\link{FD.generate}}, \code{\link{FD.amalgamation}}, \code{\link{FD.normalization}}
#'
#' @export

FD.subcomposition <- function(data, columns) {
  if (any(data < 0) )
    stop("Error: Values must be positive") else {
  temp <- data[, columns]
  temp/rowSums(temp)
    }
}

# ---------------------------- FD.amalgamation ----------------------------

#' @title Amalgamation
#'
#' @description Given a matrix or a numeric dataframe, this function returns a composition where a set of specified columns is amalgamated together.
#' The compositional operation of amalgamation  provides sums of composition elements aimed at grouping homogeneous parts of the whole.
#'
#' @param data a matrix or a dataframe containing only variables to be transformed into compositional variables, after amalgamation.
#' @param columns numeric vector containing the position of the columns to be amalgamated together.
#' @param name string containing the name of the new column resulted from the amalgamation.
#'
#' @details Values must be positive. In case one row-entry (or more) is NA, the whole row will be returned as NA.

#'
#' @examples
#' data(oliveoil)
#' dataoil <- oliveoil
#' head(dataoil)
#' data <- FD.normalization(dataoil[,3:10])
#' head(data)
#' data.sub <- FD.subcomposition(data,c(1,3,4,5))
#' head(data.sub)
#' data.amalg <- FD.amalgamation(data,c(2,6,7,8),name='others')
#' head(data.amalg)
#'
#' @seealso \code{\link{FD.generate}}, \code{\link{FD.subcomposition}}, \code{\link{FD.normalization}}
#'
#' @export

FD.amalgamation <- function(data, columns, name = NULL) {
  amalg <- rowSums(data[, columns])
  output <- cbind(data[, -columns], amalg)
  D <- dim(output)[2]
  if (any(data < 0) )
    stop("Error: Values must be positive") else {
  if (!is.null(name)) 
    colnames(output)[D] <- name
  output/rowSums(output)
    }
}

# ---------------------------- FD.normalization ----------------------------

#' @title Normalization
#'
#' @description Given a matrix or a numeric dataframe, this function returns a composition (i.e. data summing up to 1).
#'
#' @details Values must be positive. In case one row-entry (or more) is NA, the whole row will be returned as NA.
#'
#' @param data a matrix or a dataframe containing only variables to be transformed into compositional variables.
#'
#' @examples
#' data(oliveoil)
#' dataoil <- oliveoil
#' head(dataoil)
#' data <- FD.normalization(dataoil[,3:10])
#' head(data)
#' data.sub <- FD.subcomposition(data,c(1,3,4,5))
#' head(data.sub)
#' data.amalg <- FD.amalgamation(data,c(2,6,7,8),name='others')
#' head(data.amalg)
#'
#' @seealso \code{\link{FD.generate}}, \code{\link{FD.subcomposition}}, \code{\link{FD.amalgamation}}
#'
#' @export

FD.normalization <- function(data) {
  if (any(data < 0) )
    stop("Error: Values must be positive") else {
  data/rowSums(data)
    }
}