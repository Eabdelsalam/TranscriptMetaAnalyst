#' Plot download formats
#'
#' @return A named character vector of file formats (and extensions), named by
#'   their labels.
#' @export
#' @examples
#' plot_formats()
plot_formats <- function() {
  c("PDF" = "pdf", "PNG (300 dpi)" = "png", "TIFF (300 dpi)" = "tiff")
}

#' Save a plot in a chosen format and size
#'
#' Writes a ggplot, a ComplexHeatmap object, or a function that draws a plot
#' (for base graphics and pheatmap) to a PDF, PNG or TIFF file.
#'
#' @param file Output file.
#' @param plot A ggplot, a ComplexHeatmap object (e.g. an UpSet plot), or a
#'   function without arguments that draws the plot.
#' @param format One of [plot_formats()].
#' @param width,height Size in cm. `NULL` or `NA` uses `default_width` and
#'   `default_height`.
#' @param default_width,default_height Size in cm used when `width` or
#'   `height` is not given.
#' @param dpi Resolution of PNG and TIFF files.
#' @return `file`, invisibly.
#' @export
#' @examples
#' p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
#'   ggplot2::geom_point()
#' f <- tempfile(fileext = ".png")
#' save_plot(f, p, format = "png", width = 10, height = 8)
#' save_plot(tempfile(fileext = ".pdf"), function() plot(1:10))
save_plot <- function(file, plot, format = "pdf", width = NULL, height = NULL,
                      default_width = 20, default_height = 15, dpi = 300) {
  format <- match.arg(format, plot_formats())
  size <- function(x, default) if (length(x) == 1 && !is.na(x) && x > 0) x else default
  width <- size(width, default_width)
  height <- size(height, default_height)
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(file, plot, device = format, width = width, height = height,
                    units = "cm", dpi = dpi)
    return(invisible(file))
  }
  inches <- function(cm) cm / 2.54
  switch(format,
         pdf = grDevices::pdf(file, width = inches(width), height = inches(height)),
         png = grDevices::png(file, width = width, height = height, units = "cm", res = dpi),
         tiff = grDevices::tiff(file, width = width, height = height, units = "cm", res = dpi,
                                compression = "lzw"))
  on.exit(grDevices::dev.off())
  if (is.function(plot)) plot() else ComplexHeatmap::draw(plot)
  invisible(file)
}
