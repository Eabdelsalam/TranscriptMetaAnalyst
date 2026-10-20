# A small ExpressionSet standing in for a downloaded GEO series
fake_eset <- function(two_colour = FALSE, prefix = "GSM", n_genes = 3) {
  base <- matrix(c(5, 6, 7, 8, 9, 10,
                   7, 7, 8, 8, 9, 9,
                   1, 2, 1, 2, 1, 2), nrow = 3, byrow = TRUE)
  if (n_genes > 3) {
    set.seed(1)
    extra <- matrix(stats::rnorm((n_genes - 3) * 6, 8, 1), ncol = 6)
    extra[1:5, 3:4] <- extra[1:5, 3:4] + 3
    base <- rbind(base, extra)
  }
  expr <- base
  dimnames(expr) <- list(paste0("p", seq_len(nrow(expr))), paste0(prefix, 1:6))
  pheno <- data.frame(title = paste0("sample", 1:6),
                      source_name_ch1 = "culture",
                      characteristics_ch1 = rep(c("treatment: control", "treatment: stress", "treatment: other"), each = 2),
                      label_ch1 = if (two_colour) rep(c("Cy3", "Cy5"), 3) else "biotin",
                      channel_count = "1", row.names = colnames(expr))
  orf <- c("sll0001", "sll0002", "", if (n_genes > 3) sprintf("slr%04d", seq_len(n_genes - 3)))
  feature <- data.frame(ID = rownames(expr), ORF = orf,
                        Gene.symbol = c("abcA", "abcB", "", if (n_genes > 3) sprintf("gene%d", seq_len(n_genes - 3))),
                        row.names = rownames(expr))
  Biobase::ExpressionSet(expr, phenoData = Biobase::AnnotatedDataFrame(pheno),
                         featureData = Biobase::AnnotatedDataFrame(feature), annotation = "GPL123")
}
