## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----eval = FALSE-------------------------------------------------------------
# library(bgms)
# data = Wenchuan[, 1:9]
# fit = bgm(data, seed = 1234)

## ----include = FALSE----------------------------------------------------------
library(bgms)
data = Wenchuan[, 1:9]
fit = bgm(data,
  seed = 1234, chains = 2, cores = 2,
  display_progress = "none", verbose = FALSE
)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
verdicts(fit)

## ----fig.width = 12, fig.height = 4.5, out.width = "100%"---------------------
plot(fit)

## ----eval = FALSE-------------------------------------------------------------
# fit = bgm(data,
#   interaction_prior = cauchy_prior(scale = 2.5),
#   edge_prior = beta_bernoulli_prior(alpha = 1, beta = 1)
# )

## ----eval = FALSE-------------------------------------------------------------
# # 200 draws from a five-variable Gaussian graphical model whose precision
# # matrix is a chain, 1-2-3-4-5.
# set.seed(1234)
# K = diag(5)
# K[cbind(1:4, 2:5)] = -0.4
# K[cbind(2:5, 1:4)] = -0.4
# continuous_data = matrix(rnorm(200 * 5), 200, 5) %*% chol(solve(K))
# 
# fit_ggm = bgm(continuous_data, variable_type = "continuous", seed = 1234)
# verdicts(fit_ggm)

## ----eval = FALSE-------------------------------------------------------------
# fit_compare = bgmCompare(
#   x = ADHD[ADHD$group == 1, 2:6],
#   y = ADHD[ADHD$group == 0, 2:6],
#   seed = 1234
# )
# verdicts(fit_compare)
# plot(fit_compare)

