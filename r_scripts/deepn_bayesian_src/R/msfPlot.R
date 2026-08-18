msfPlot <- function(Data, type=c("line", "point"), lwd=2, pch=16, col="#4339824D",
                    bty="n", xlab="Mean (Unselected)", ylab="SD (Fold change)", ...) {
  type <- match.arg(type)
  V <- rpm(Data)$Vector
  m1 <- apply(V[,1,], 1, mean)
  FC <- V[,2,]/V[,1,]
  vF <- apply(FC, 1, var)
  sF <- sqrt(vF)

  if (type=="line") {
    plot(sqrt(m1), sqrt(sF), lwd=lwd, col=col, type='h', bty=bty,
         xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", ...)
    sqrtAxis()
  } else if (type=="point") {
    plot(sqrt(m1), sqrt(sF), col=col, pch=pch, bty='n',
         xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", ...)
    sqrtAxis()
  }
}
