test_that("method and format choices have labels and codes", {
  expect_setequal(unname(batch_methods()), c("model", "combat_seq", "none"))
  expect_setequal(unname(plot_formats()), c("pdf", "png", "tiff"))
  expect_true(all(nzchar(names(batch_methods()))))
  expect_true(all(nzchar(names(plot_formats()))))
})

test_that("enrichment_label names service, version, organism and background", {
  label <- enrichment_label(list(service = "STRING", version = "12.0", organism = "Synechocystis (1148)",
                                 background = "all 10 tested genes (9 recognised by STRING)"))
  expect_identical(label, paste0("STRING (12.0); organism: Synechocystis (1148); ",
                                 "background: all 10 tested genes (9 recognised by STRING)."))
})

test_that("write_tma_script writes the lines of tma_script", {
  settings <- list(de_method = "DESeq2", meta_method = "fisher", batch_method = "none",
                   batch_plots = FALSE, reference = NULL)
  file <- withr::local_tempfile(fileext = ".R")
  expect_identical(write_tma_script(file, settings, data = "demo"), file)
  expect_identical(readLines(file), tma_script(settings, data = "demo"))
})

test_that("tma_geo_cache creates the cache folder in the user cache directory", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  dir <- tma_geo_cache()
  expect_true(dir.exists(dir))
  expect_identical(basename(dir), "geo")
  expect_true(startsWith(normalizePath(dir), normalizePath(Sys.getenv("R_USER_CACHE_DIR"))))
})

test_that("tma_use_background_workers sets a parallel plan only when the plan is sequential", {
  old <- future::plan(future::sequential)
  withr::defer(future::plan(old))
  tma_use_background_workers(workers = 2)
  expect_false(inherits(future::plan(), "sequential"))
  expect_true(inherits(future::plan(), c("multicore", "multisession")))

  # a plan chosen by the user is kept
  future::plan(future::multisession, workers = 2)
  tma_use_background_workers(workers = 3)
  expect_true(inherits(future::plan(), "multisession"))
  expect_identical(future::nbrOfWorkers(), 2L)
})
