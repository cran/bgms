% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgm.em.R
\name{bgm.em}
\alias{bgm.em}
\title{EM variable selection for a Markov Random Field model for ordinal variables.}
\usage{
bgm.em(
  x,
  precision = 0.975,
  convergence_criterion = sqrt(.Machine$double.eps),
  theta = 0.5,
  hierarchical = FALSE,
  indicator_alpha = 1,
  indicator_beta = 1,
  maximum_iterations = 1000,
  threshold_alpha = 1,
  threshold_beta = 1
)
}
\arguments{
\item{x}{A matrix with \code{n} rows and \code{p} columns, containing binary
and ordinal variables for \code{n} independent observations and \code{p}
variables in the network. Variables are recoded as non-negative integers
\code{(0, 1, ..., m)} if not done already. Unobserved categories are
collapsed into other categories after recoding. See \code{reformat_data} for
details.}

\item{precision}{A value between 0 and 1 representing the desired precision
for edge selection, equal to one minus the desired type-1 error rate. Default
is 0.975.}

\item{convergence_criterion}{The criterion for the pseudoposterior values'
convergence in the EM algorithm. Default is \code{sqrt(.Machine$double.eps)}.}

\item{theta}{The prior inclusion probability. A value of \code{0.5}, combined
with \code{hierarchical = FALSE}, specifies a uniform prior on the network
structures' space.}

\item{hierarchical}{If TRUE, a beta prior distribution with hyperparameters
\code{indicator_alpha, indicator_beta} is imposed on the prior inclusion
probability \code{theta}. A uniform prior on inclusion probability, using a
beta with \code{indicator_alpha = indicator_beta = 1}, specifies a uniform
prior on network structure complexity.}

\item{indicator_alpha, indicator_beta}{Hyperparameters for the beta prior
distribution on the prior inclusion probability theta when
\code{hierarchical = TRUE}. Default is 1.}

\item{maximum_iterations}{Maximum number of EM iterations used. Default is
1e3. A warning appears if the procedure hasn't converged within the maximum
number of iterations.}

\item{threshold_alpha, threshold_beta}{Shape parameters for the Beta-prime
prior on thresholds. Default is 1.}
}
\value{
A list containing:
\itemize{
\item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
containing the pairwise association estimates in the off-diagonal elements.
\item \code{gamma}: A matrix with \code{p} rows and \code{p} columns,
containing the expected values of edge inclusion variables (local posterior
probabilities of edge inclusion).
\item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
columns, containing the category thresholds for each node.
\item \code{theta} (if \code{hierarchical == TRUE}): A numeric value
representing the modal estimate of the prior inclusion probability.
}
}
\description{
The function \code{bgm.em} selects promising edges for the ordinal
MRF using the joint pseudolikelihood and a continuous spike and slab prior
distribution stipulated on the MRF's interaction or association parameters.
}
\examples{
\donttest{
 #Store user par() settings
 op <- par(no.readonly = TRUE)

 ##Analyse the Wenchuan dataset
 fit = bgm.em(x = Wenchuan)


 #------------------------------------------------------------------------------|
 # INCLUSION - EDGE WEIGHT PLOT
 #------------------------------------------------------------------------------|

 par(mar = c(6, 5, 1, 1))
 plot(x = fit$interactions[lower.tri(fit$interactions)],
      y = fit$gamma[lower.tri(fit$gamma)], ylim = c(0, 1),
      xlab = "", ylab = "", axes = FALSE, pch = 21, bg = "#bfbfbf", cex = 1.3)
 abline(h = 0, lty = 2, col = "#bfbfbf")
 abline(h = 1, lty = 2, col = "#bfbfbf")
 abline(h = .5, lty = 2, col = "#bfbfbf")
 mtext("Posterior Inclusion Probability", side = 1, line = 3, cex = 1.7)
 mtext("Posterior Mode Edge Weight", side = 2, line = 3, cex = 1.7)
 axis(1)
 axis(2, las = 1)


 #------------------------------------------------------------------------------|
 # THE LOCAL MEDIAN PROBABILITY NETWORK
 #------------------------------------------------------------------------------|

 library(qgraph) #For plotting the estimated network

 posterior.inclusion = fit$gamma[lower.tri(fit$gamma)]
 tmp = fit$interactions[lower.tri(fit$interactions)]
 tmp[posterior.inclusion < 0.5] = 0

 median.prob.model = matrix(0, nrow = ncol(Wenchuan), ncol = ncol(Wenchuan))
 median.prob.model[lower.tri(median.prob.model)] = tmp
 median.prob.model = median.prob.model + t(median.prob.model)

 rownames(median.prob.model) = colnames(Wenchuan)
 colnames(median.prob.model) = colnames(Wenchuan)

 qgraph(median.prob.model,
        theme = "TeamFortress",
        maximum = .5,
        fade = FALSE,
        color = c("#f0ae0e"), vsize = 10, repulsion = .9,
        label.cex = 1.1, label.scale = "FALSE",
        labels = colnames(Wenchuan))

 #Restore user par() settings
 par(op)
}
}
