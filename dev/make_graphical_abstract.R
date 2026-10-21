# Exports the graphical abstract shown in the README,
# manual/TMA_graphical_abstract.pptx (one A4-portrait slide, edited in
# PowerPoint), to manual/TMA_graphical_abstract.png at 300 dpi (2481 x 3508 px):
# LibreOffice converts the slide to PDF, and pdftoppm renders the PDF.
# Needs LibreOffice (soffice) and poppler (pdftoppm). Run from the package root
# after editing the slide:
#
#   Rscript dev/make_graphical_abstract.R

source_pptx <- normalizePath(file.path("manual", "TMA_graphical_abstract.pptx"))
output_png <- file.path(normalizePath("manual"), "TMA_graphical_abstract.png")
for (tool in c("soffice", "pdftoppm")) {
  if (!nzchar(Sys.which(tool))) stop(tool, " is needed to export the abstract.", call. = FALSE)
}

work <- tempfile("abstract")
dir.create(work)
file.copy(source_pptx, work)
pptx <- file.path(work, basename(source_pptx))
# R's LD_LIBRARY_PATH stops LibreOffice from finding its own libraries
status <- system2("env", c("-u", "LD_LIBRARY_PATH", "soffice", "--headless", "--convert-to", "pdf",
                           "--outdir", shQuote(work), shQuote(pptx)),
                  stdout = FALSE, stderr = FALSE)
pdf <- sub("\\.pptx$", ".pdf", pptx)
if (status != 0 || !file.exists(pdf)) stop("LibreOffice could not convert the slide.", call. = FALSE)

status <- system2("pdftoppm", c("-png", "-r", "300", "-singlefile", shQuote(pdf),
                                shQuote(file.path(work, "abstract"))))
if (status != 0) stop("pdftoppm could not render the slide.", call. = FALSE)
file.copy(file.path(work, "abstract.png"), output_png, overwrite = TRUE)
unlink(work, recursive = TRUE)
message("Graphical abstract written to ", output_png)
