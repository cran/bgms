## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## ----eval=FALSE---------------------------------------------------------------
# library(bgms)
# data = Wenchuan[, 1:5]
# fit = bgm(data, seed = 1234)

## ----include=FALSE------------------------------------------------------------
library(bgms)
data = Wenchuan[, 1:5]
fit = bgm(data, seed = 1234, chains = 2, display_progress = "none", verbose = FALSE)

## -----------------------------------------------------------------------------
summary(fit)$pairwise

## -----------------------------------------------------------------------------
summary(fit)$indicator

## ----fig.width= 7, fig.height= 7----------------------------------------------
param_index = 1
chains = fit$raw_samples$pairwise
nchains = length(chains)
cols = c("firebrick", "steelblue", "darkgreen", "goldenrod")

plot(chains[[1]][, param_index],
  type = "l", col = cols[1],
  xlab = "Iteration", ylab = "Value",
  main = "Traceplot of pairwise[1]",
  ylim = range(sapply(chains, function(ch) range(ch[, param_index])))
)
if(nchains > 1) {
  for(c in 2:nchains) {
    lines(chains[[c]][, param_index], col = cols[c])
  }
}

## -----------------------------------------------------------------------------
coef(fit)$indicator

## -----------------------------------------------------------------------------
# Example for one edge
p = coef(fit)$indicator[1, 5]
BF_10 = p / (1 - p)
BF_10

## -----------------------------------------------------------------------------
1 / BF_10

## -----------------------------------------------------------------------------
fit$nuts_diag$summary

## -----------------------------------------------------------------------------
fit$nuts_diag$warmup_check

