mcvPlot <- function(Data, type=c("line", "point"), lwd=2, pch=16, col="#4339824D",
                    bty="n", xlab="Mean (Unselected)", ylab="SD (Selected) / Mean (Unselected)", ...) {
  type <- match.arg(type)
  V <- rpm(Data)$Vector
  m1 <- apply(V[,1,], 1, mean)
  v2 <- apply(V[,2,], 1, sd)
  s2 <- sqrt(v2)

  if (type=="line") {
    plot(sqrt(m1), sqrt(s2/m1), lwd=lwd, col=col, type='h', bty=bty,
         xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", ...)
    sqrtAxis()
  } else if (type=="point") {
    plot(sqrt(m1), sqrt(s2/m1), col=col, pch=pch, bty='n',
         xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", ...)
    sqrtAxis()
  }
}
