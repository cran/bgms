## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 6
)

## ----eval = FALSE-------------------------------------------------------------
# library(bgms)
# data = Wenchuan

## ----include = FALSE----------------------------------------------------------
library(bgms)
data = Wenchuan

## ----eval = FALSE-------------------------------------------------------------
# fit = bgm(data[, 1:6], seed = 1234, chains = 2, cores = 2,
#           display_progress = "none", verbose = FALSE)
# ps = prior_sensitivity_check(fit, seed = 1234, iter = 3000, warmup = 1000,
#                              cores = 2)
# ps

## ----echo = FALSE-------------------------------------------------------------
# The report below is the real object that call produced, run once with those
# exact arguments by data-raw/make-prior-sensitivity.R and stored
# beside this file. Nothing here is transcribed: the check is not re-run on
# every build, but what is printed is still `print()` on a genuine
# prior_sensitivity_check() result.
ps = readRDS("prior-sensitivity-ps.rds")
ps

## ----eval = FALSE-------------------------------------------------------------
# plot(ps)

## ----eval = FALSE-------------------------------------------------------------
# head(ps$edges)

