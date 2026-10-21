test_that("runTMA sets the stop-on-session-end option only while the app runs", {
  seen <- NULL
  local_mocked_bindings(
    runApp = function(...) seen <<- getOption("tma.stop_on_session_end"),
    .package = "shiny"
  )
  withr::local_options(tma.stop_on_session_end = NULL)

  runTMA()

  expect_true(seen)
  expect_null(getOption("tma.stop_on_session_end"))
})

