# Extracted from test-zratio-surface-build.R:93

# test -------------------------------------------------------------------------
skip_on_cran()
withr::local_options(
  bgms.zratio_surface_cache = FALSE,
  bgms.correction_table_cache = FALSE
)
zc = bgms:::zratio_constants(0.5 * log(12), 3)
s1 = bgms:::build_surfaces_allmc(zc, max_size = 8L, cores = 1L)
withr::local_options(bgms.zratio_surface_psock = TRUE)
s2 = bgms:::build_surfaces_allmc(zc, max_size = 8L, cores = 2L)
expect_identical(s1, s2)
