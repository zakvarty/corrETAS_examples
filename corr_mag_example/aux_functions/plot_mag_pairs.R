
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


plot_mag_pairs <- function(m ,b, ...){
  triggered_indices <- which(b>0)
  parent_indices  <- b[triggered_indices]
  parent_distribution <- b[parent_indices] > 0

  triggered_mags <- m[triggered_indices]
  parent_mags <- m[parent_indices]

  plot(x = parent_mags, y = triggered_mags, pch = 16, col = parent_distribution + 1, ...)
  # Add legend to top right, outside plot region
  # add_legend("topright", legend=c("parent background", "parent triggered"), pch=16,
  #            col=c("black", "red"),
  #            horiz=TRUE, bty='n', cex=1)
}

