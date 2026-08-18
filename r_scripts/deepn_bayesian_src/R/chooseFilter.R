chooseFilter <- function(Data, minRPM=1, maxRPM=100, plot=FALSE, ...) {
  RPM <- rpm(Data)
  thresh <- exp(seq(log(minRPM), log(maxRPM), length=19))
  g <- o <- numeric(length(thresh))
  for (j in 1:length(thresh)) {
    Data <- applyFilter(Data, thresh=thresh[j])
    g[j] <- nrow(Data$Vector)
    o[j] <- overdisp(Data)[2]
  }
  if (plot) {
    plot(g, o, type='l', bty='n', las=1, lwd=3, xlab="Genes that pass", ylab="Overdispersion")
    f <- approxfun(g, thresh)
    at <- seq(min(g), max(g), length=6)
    lab <- formatC(f(at), format='f', digits=1)
    axis(3, at=at, labels=lab)
  }
  #print(cbind(thresh, g2, o2, (g2 - min(g2)) / (max(g2) - min(g2)), (max(o2) - o2) / (max(o2) - min(o2))))
  #score <- (max(o2) - o2) / (max(o2) - min(o2)) + (g2 - min(g2)) / (max(g2) - min(g2))
  #thresh[which.max(score)]
  invisible(data.frame(Threshold=thresh, Genes=g, Overdispersion=o))
}
