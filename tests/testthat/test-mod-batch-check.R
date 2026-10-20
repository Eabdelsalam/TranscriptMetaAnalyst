test_that("mod_batch_check_server runs the check on click and shows the advice", {
  sim <- batch_data_set()
  data <- shiny::reactive(list(counts = sim$counts, design = sim$design))
  shiny::testServer(mod_batch_check_server, args = list(data = data), {
    expect_error(check(), class = "shiny.silent.error")
    session$setInputs(run = 1)
    expect_s3_class(check(), "tma_batch_check")
    expect_match(output$table, "Batch explains as much variance as condition or more")
    expect_match(output$table, "Variance explained by batch")
    expect_s3_class(pca(), "ggplot")
    expect_s3_class(rle(), "ggplot")
  })
})

test_that("the RNA-seq tab checks the uploaded files for batch effects", {
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_counts = list(datapath = tma_example_file("counts")),
                      rnaseq_design = list(datapath = tma_example_file("design")),
                      `rnaseq_batch_check-run` = 1)
    out <- output$`rnaseq_batch_check-table`
    expect_match(out, "no batch column")
    expect_match(out, "No sign of a batch effect")
  })
})
