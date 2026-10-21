# Builds manual/TMA_Manual.pdf from manual/TMA_Manual.md: pandoc writes a
# self-contained HTML page (figures embedded), and headless Chrome prints it to
# PDF. Needs pandoc and Chrome or Chromium. Run from the package root after
# changing the manual or its figures:
#
#   Rscript dev/make_manual_pdf.R

manual_dir <- normalizePath("manual")
source_md <- file.path(manual_dir, "TMA_Manual.md")
output_pdf <- file.path(manual_dir, "TMA_Manual.pdf")

chrome <- Sys.getenv("CHROMOTE_CHROME", unset = "")
if (!nzchar(chrome)) chrome <- Sys.which(c("google-chrome", "chromium", "chromium-browser"))
chrome <- chrome[nzchar(chrome)][1]
if (is.na(chrome)) stop("Chrome or Chromium is needed to print the PDF.", call. = FALSE)
if (!nzchar(Sys.which("pandoc"))) stop("pandoc is needed to convert the manual.", call. = FALSE)

# print layout: A4 pages; headings stay with the text that follows them, and
# figures and table rows are not split across pages
css <- "
@page { size: A4; margin: 18mm 16mm; }
body { font-family: 'Source Sans Pro', 'Helvetica Neue', Arial, sans-serif; font-size: 10.5pt;
       line-height: 1.45; color: #222; max-width: none; margin: 0; }
h1 { font-size: 22pt; margin-top: 0; }
h2 { font-size: 16pt; border-bottom: 1px solid #ccc; padding-bottom: 2px; margin-top: 1.4em; }
h3 { font-size: 13pt; margin-top: 1.2em; }
h4 { font-size: 11pt; }
h2, h3, h4 { break-after: avoid; }
img { max-width: 100%; max-height: 210mm; display: block; margin: 0.6em auto; }
p:has(img), table:has(img) { break-inside: avoid; }
img[alt='Figure 4'] { max-width: 45%; }
table { border-collapse: collapse; margin: 0.6em 0; font-size: 9.5pt; break-inside: auto; }
th, td { border: 1px solid #ccc; padding: 3px 6px; vertical-align: top; }
tr { break-inside: avoid; }
th { background: #f2f2f2; }
code { font-size: 9pt; background: #f4f4f4; padding: 0 2px; }
pre { background: #f4f4f4; padding: 6px 8px; font-size: 9pt; white-space: pre-wrap; break-inside: avoid; }
blockquote { margin: 0.6em 0; padding: 0.2em 1em; border-left: 3px solid #4a90c2; color: #444; }
a { color: #1f5f99; text-decoration: none; }
"

work <- tempfile("manual-pdf")
dir.create(work)
css_file <- file.path(work, "manual.css")
html_file <- file.path(work, "TMA_Manual.html")
writeLines(css, css_file)

status <- system2("pandoc", c(
  shQuote(source_md), "--from", "gfm", "--to", "html5", "--standalone", "--embed-resources",
  "--resource-path", shQuote(manual_dir), "--css", shQuote(css_file),
  "--metadata", shQuote("pagetitle=TranscriptMetaAnalyst (TMA) User Manual"),
  "--output", shQuote(html_file)
))
if (status != 0) stop("pandoc failed.", call. = FALSE)

status <- system2(chrome, c(
  "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
  paste0("--user-data-dir=", shQuote(file.path(work, "chrome"))),
  paste0("--print-to-pdf=", shQuote(output_pdf)),
  shQuote(paste0("file://", html_file))
))
if (status != 0 || !file.exists(output_pdf)) stop("Chrome could not print the PDF.", call. = FALSE)

unlink(work, recursive = TRUE)
message("Manual written to ", output_pdf)
