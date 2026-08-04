## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## -----------------------------------------------------------------------------
library(bgms)

data_adhd = ADHD[ADHD$group == 1, -1]
data_adhd = data_adhd[, 1:5]
data_no_adhd = ADHD[ADHD$group == 0, -1]
data_no_adhd = data_no_adhd[, 1:5]

## ----eval = FALSE-------------------------------------------------------------
# fit = bgmCompare(x = data_adhd, y = data_no_adhd, seed = 1234)

## ----include=FALSE------------------------------------------------------------
fit = bgmCompare(
  x = data_adhd, y = data_no_adhd, seed = 1234, chains = 2,
  display_progress = "none", verbose = FALSE
)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
coef(fit)

## ----fig.width = 10, fig.height = 5, out.width = "100%"-----------------------
plot(fit, type = "groups")

