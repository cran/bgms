# Extracted from test-methods.R:38

# setup ------------------------------------------------------------------------
library(testthat)
test_env = simulate_test_env(package = "bgms", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
specs = get_bgms_fixtures()
labels = vapply(specs, `[[`, character(1), "label")
required = c("binary", "ordinal", "ggm", "mixed-mrf", "blume-capel")
for(r in required) {
  expect_true(r %in% labels, info = sprintf("missing required label '%s'", r))
}
expect_equal(length(specs), 24L,
  info = "bgms fixture count changed — update this guard if intentional"
)
