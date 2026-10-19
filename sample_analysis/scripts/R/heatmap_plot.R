# Explicit figure coordinates preserve the top margin for rotated labels.
draw_selected_heatmap <- function(z, column_labels, row_labels, colors, breaks, family) {
  par(fig = c(0.17, 1, 0, 1), mar = c(0.8, 0.25, 4.5, 5.5),
      family = family, cex = 1, xpd = FALSE)
  image(
    x = seq_len(ncol(z)), y = seq_len(nrow(z)),
    z = t(z[nrow(z):1, , drop = FALSE]), col = colors, breaks = breaks,
    axes = FALSE, xlab = "", ylab = "", main = "",
    xlim = c(0.5, ncol(z) + 0.5), ylim = c(0.5, nrow(z) + 0.5),
    xaxs = "i", yaxs = "i"
  )
  axis(3, at = seq_len(ncol(z)), labels = column_labels, las = 2,
       cex.axis = 0.62, col.axis = "#000000", line = 0.1, tick = FALSE)
  axis(4, at = seq_len(nrow(z)), labels = row_labels, las = 1,
       cex.axis = 0.52, col.axis = "#000000", line = 0.1, tick = FALSE)
  box(lwd = 0.65)

  par(fig = c(0, 0.17, 0.18, 0.70), mar = c(0, 0, 0, 0),
      new = TRUE, family = family, cex = 1, xpd = FALSE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-3.2, 3.2), xaxs = "i", yaxs = "i")
  legend_y <- seq(min(breaks), max(breaks), length.out = length(colors) + 1)
  for (i in seq_along(colors)) {
    rect(0.55, legend_y[i], 0.80, legend_y[i + 1], col = colors[i], border = NA)
  }
  rect(0.55, min(breaks), 0.80, max(breaks), border = "#1F252B", lwd = 0.6)
  for (tick in c(-2.5, 0, 2.5)) text(0.48, tick, labels = format(tick), adj = c(1, 0.5), cex = 0.52)
  text(0.65, 2.95, labels = "Row Z", cex = 0.58)
}
