# Extracted from test-mixed-gradient.R:178

# setup ------------------------------------------------------------------------
library(testthat)
test_env = simulate_test_env(package = "bgms", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
mixed_fd_gradient = function(params, x, y, num_cats, is_ord, base_cat,
                             edge_ind, pl_mode, scale, eps = 1e-5) {
  n_total = length(params)
  fd = numeric(n_total)
  for(k in seq_len(n_total)) {
    p_plus = params
    p_minus = params
    p_plus[k] = p_plus[k] + eps
    p_minus[k] = p_minus[k] - eps
    fp = mixed_test_logp_and_gradient(
      p_plus, x, y, num_cats, as.integer(is_ord),
      base_cat, edge_ind, pl_mode, scale
    )$value
    fm = mixed_test_logp_and_gradient(
      p_minus, x, y, num_cats, as.integer(is_ord),
      base_cat, edge_ind, pl_mode, scale
    )$value
    fd[k] = (fp - fm) / (2 * eps)
  }
  fd
}
mixed_check_gradient = function(params, x, y, num_cats, is_ord, base_cat,
                                edge_ind, pl_mode, scale, eps = 1e-5) {
  res = mixed_test_logp_and_gradient(
    params, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, pl_mode, scale
  )
  fd = mixed_fd_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, pl_mode, scale, eps
  )
  ag = res$gradient
  denom = pmax(abs(ag), abs(fd), 1)
  rel_err = abs(ag - fd) / denom
  max(rel_err)
}

# test -------------------------------------------------------------------------
set.seed(42)
n = 60
p = 2
q = 3
x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
y = matrix(rnorm(n * q), n, q)
num_cats = rep(2L, p)
is_ord = rep(1L, p)
base_cat = rep(0L, p)
total = p + q
edge_ind = matrix(1L, total, total)
diag(edge_ind) = 0L
n_main = sum(num_cats)
n_pw = p * (p - 1) / 2
n_chol = q * (q + 1) / 2
set.seed(1)
params = rnorm(n_main + n_pw + q + p * q + n_chol, sd = 0.3)
err = mixed_check_gradient(
  params, x, y, num_cats, is_ord, base_cat,
  edge_ind, 2.5
)
