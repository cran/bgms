## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 3.4
)

## ----eval=FALSE---------------------------------------------------------------
# library(bgms)
# data = Wenchuan[, 1:9]
# fit = bgm(data, seed = 1234)

## ----include=FALSE------------------------------------------------------------
library(bgms)
data = Wenchuan[, 1:9]
fit = bgm(data,
  seed = 1234, chains = 2, iter = 1000, warmup = 1000,
  display_progress = "none", verbose = FALSE
)

## -----------------------------------------------------------------------------
edges = verdicts(fit)
edges

## -----------------------------------------------------------------------------
edges[edges$fragile, c("parameter", "log_bf", "verdict")]

## ----fig.width = 12, fig.height = 4.5, out.width = "100%"---------------------
plot(fit)

## ----fig.height = 3.6, fig.width = 5------------------------------------------
plot_edge_posterior(fit, "intrusion", "upset")

## ----fig.height = 3.6, fig.width = 6------------------------------------------
strength = extract_centrality(fit)
head(summary(strength))
plot(strength)

## -----------------------------------------------------------------------------
calibration = calibration_check(fit, nrep = 200, seed = 1)
calibration

## ----fig.height = 5.5, fig.width = 7------------------------------------------
plot(calibration)

## -----------------------------------------------------------------------------
observed = as.matrix(data[complete.cases(data), ])

replicates = simulate(fit,
  nsim = nrow(observed), method = "posterior-sample",
  ndraws = 200, seed = 1, display_progress = "none"
)

scores = 0:(ncol(observed) * max(observed))
distribution = function(x) {
  as.numeric(table(factor(rowSums(x), levels = scores))) / nrow(x)
}

observed_share = distribution(observed)
replicated_share = vapply(replicates, distribution, observed_share)
band = apply(replicated_share, 1, quantile, probs = c(0.025, 0.975))

## ----fig.height = 3.8, fig.width = 6------------------------------------------
plot(scores, observed_share,
  type = "n", las = 1, xlab = "Sum score", ylab = "Share of cases",
  ylim = range(c(band, observed_share))
)
polygon(c(scores, rev(scores)), c(band[1, ], rev(band[2, ])),
  col = adjustcolor("grey55", 0.25), border = NA
)
lines(scores, observed_share, col = "#0072B2", lwd = 1.8)

mean(observed_share < band[1, ] | observed_share > band[2, ])

