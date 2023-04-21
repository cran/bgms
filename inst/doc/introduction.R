## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = identical(Sys.getenv("BUILD_VIGNETTE"), "true")
)

## ----setup--------------------------------------------------------------------
#  library(bgms)

## ----usage, eval=FALSE--------------------------------------------------------
#  bgm.em(
#    x,
#    precision = 0.975,
#    convergence_criterion = sqrt(.Machine$double.eps),
#    theta = 0.5,
#    hierarchical = FALSE,
#    indicator_alpha = 1,
#    indicator_beta = 1,
#    maximum_iterations = 1000,
#    threshold_alpha = 1,
#    threshold_beta = 1
#  )

## ----message=FALSE, warning=FALSE, cache=TRUE---------------------------------
#  fit <- bgm.em(x = Wenchuan)

## ----message=FALSE, warning=FALSE, fig.width= 7, fig.height= 7----------------
#  par(mar = c(6, 5, 1, 1))
#  plot(x = fit$interactions[lower.tri(fit$interactions)],
#      y = fit$gamma[lower.tri(fit$gamma)], ylim = c(0, 1),
#      xlab = "", ylab = "", axes = FALSE, pch = 21, bg = "#bfbfbf", cex = 1.3)
#  abline(h = 0, lty = 2, col = "#bfbfbf")
#  abline(h = 1, lty = 2, col = "#bfbfbf")
#  abline(h = .5, lty = 2, col = "#bfbfbf")
#  mtext("Posterior mode edge weight", side = 1, line = 3, cex = 1.7)
#  mtext("(Local) posterior inclusion probability", side = 2, line = 3, cex = 1.7)
#  axis(1)
#  axis(2, las = 1)

## ----message=FALSE, warning=FALSE, fig.width= 7, fig.height= 7----------------
#  library(qgraph) #For plotting the estimated network
#  
#  posterior.inclusion <- fit$gamma[lower.tri(fit$gamma)]
#  tmp <- fit$interactions[lower.tri(fit$interactions)]
#  tmp[posterior.inclusion < 0.5] = 0
#  
#  median.prob.model <- matrix(0, nrow = ncol(Wenchuan), ncol = ncol(Wenchuan))
#  median.prob.model[lower.tri(median.prob.model)] <- tmp
#  median.prob.model <- median.prob.model + t(median.prob.model)
#  
#  rownames(median.prob.model) <- colnames(Wenchuan)
#  colnames(median.prob.model) <- colnames(Wenchuan)
#  
#  qgraph(median.prob.model,
#        theme = "TeamFortress",
#        maximum = .5,
#        fade = FALSE,
#        color = c("#f0ae0e"), vsize = 10, repulsion = .9,
#        label.cex = 1.1, label.scale = "FALSE",
#        labels = colnames(Wenchuan))

## ----usage2, eval=FALSE-------------------------------------------------------
#  bgm(
#    x,
#    iter = 10000,
#    burnin = 1000,
#    interaction_prior = c("UnitInfo", "Cauchy"),
#    cauchy_scale = 2.5,
#    threshold_alpha = 1,
#    threshold_beta = 1,
#    save = FALSE,
#    display_progress = TRUE
#  )

## ----message=FALSE, warning=FALSE, cache=TRUE---------------------------------
#  fit <-  bgm(x = Wenchuan)

## ----message=FALSE, warning=FALSE, fig.width= 7, fig.height= 7----------------
#  par(mar = c(6, 5, 1, 1))
#  plot(x = fit$interactions[lower.tri(fit$interactions)],
#       y = fit$gamma[lower.tri(fit$gamma)], ylim = c(0, 1),
#       xlab = "", ylab = "", axes = FALSE, pch = 21, bg = "gray", cex = 1.3)
#  abline(h = 0, lty = 2, col = "gray")
#  abline(h = 1, lty = 2, col = "gray")
#  abline(h = .5, lty = 2, col = "gray")
#  mtext("Posterior mean edge weight", side = 1, line = 3, cex = 1.7)
#  mtext("Posterior inclusion probability", side = 2, line = 3, cex = 1.7)
#  axis(1)
#  axis(2, las = 1)

## ----message=FALSE, warning=FALSE, fig.width= 7, fig.height= 7----------------
#  library(qgraph) #For plotting the estimated network
#  
#  posterior.inclusion <- fit$gamma[lower.tri(fit$gamma)]
#  tmp <- fit$interactions[lower.tri(fit$interactions)]
#  tmp[posterior.inclusion < 0.5] = 0
#  
#  median.prob.model <- matrix(0, nrow = ncol(Wenchuan), ncol = ncol(Wenchuan))
#  median.prob.model[lower.tri(median.prob.model)] <- tmp
#  median.prob.model <- median.prob.model + t(median.prob.model)
#  
#  rownames(median.prob.model) <- colnames(Wenchuan)
#  colnames(median.prob.model) <- colnames(Wenchuan)
#  
#  qgraph(median.prob.model,
#         theme = "TeamFortress",
#         maximum = .5,
#         fade = FALSE,
#         color = c("#f0ae0e"), vsize = 10, repulsion = .9,
#         label.cex = 1.1, label.scale = "FALSE",
#         labels = colnames(Wenchuan))

## ----message=FALSE, warning=FALSE---------------------------------------------
#  # Calculate the inclusion BFs
#  prior.odds = 1
#  posterior.inclusion = fit$gamma[lower.tri(fit$gamma)]
#  posterior.odds = posterior.inclusion / (1 - posterior.inclusion)
#  log.bayesfactor = log(posterior.odds / prior.odds)
#  #The next line is used to truncate the extreme values of the Bayes factor in the plot
#  log.bayesfactor[log.bayesfactor > 5] = 5

## ----message=FALSE, warning=FALSE, fig.width= 7, fig.height= 7----------------
#  par(mar = c(5, 5, 1, 1) + 0.1)
#  plot(fit$interactions[lower.tri(fit$interactions)], log.bayesfactor, pch = 21, bg = "#bfbfbf",
#      cex = 1.3, axes = FALSE, xlab = "", ylab = "", ylim = c(-5, 5.5),
#      xlim = c(-0.5, 1.5))
#  axis(1)
#  axis(2, las = 1)
#  abline(h = log(1/10), lwd = 2, col = "#bfbfbf")
#  abline(h = log(10), lwd = 2, col = "#bfbfbf")
#  
#  text(x = 1, y = log(1 / 10), labels = "Evidence for exclusion", pos = 1,
#      cex = 1.7)
#  text(x = 1, y = log(10), labels = "Evidence for inclusion", pos = 3, cex = 1.7)
#  text(x = 1, y = 0, labels = "Weak evidence", cex = 1.7)
#  mtext("Log-inclusion Bayes factor", side = 2, line = 3, cex = 1.5, las = 0)
#  mtext("Posterior mean edge weights ", side = 1, line = 3.7, cex = 1.5, las = 0)

## ----message=FALSE, warning=FALSE, cache=TRUE---------------------------------
#  fit <-  bgm(x = Wenchuan, save = TRUE)

## ----message=FALSE, warning=FALSE, fig.width= 7, fig.height= 7----------------
#  den = density(fit$interactions[,1], bw = "SJ")
#  i = which.min(abs(den$x - mean(fit$interactions[,1])))[1]
#  x = den$x[i]
#  f = den$y[i]
#  par(cex.main = 1.5, mar = c(5, 6, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
#  plot(den, axes = FALSE, xlab="", ylab="", main = "", frame.plot = FALSE)
#  axis(1)
#  axis(2)
#  par(las = 0)
#  mtext(text = "Edge weight", side = 1, line = 2.5, cex = 1.5)
#  mtext(text = "Posterior density", side = 2, line = 2.5, cex = 1.5)
#  # Add a point to indicate the posterior mean
#  points(x, f, pch = 21, bg = "grey", cex = 1.7)

## ----message=FALSE, warning=FALSE---------------------------------------------
#  G = 2 * fit$gamma - 1
#  S = unique(G)
#  nrow(S)

## ----message=FALSE, warning=FALSE---------------------------------------------
#  Ps = vector(length = nrow(S))
#  for(r in 1:nrow(S)) {
#    s = S[r, ]
#    tmp = G %*% s
#    Ps[r] = sum(tmp == ncol(G))
#  }
#  Ps = Ps / nrow(G) * 100
#  max(Ps)

