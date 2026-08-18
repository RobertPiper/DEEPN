overdisp <- function(Data) {
  V <- Data$Vector
  omega <- c(edgeR::estimateCommonDisp(edgeR::DGEList(V[,1,]))$common.dispersion,
             edgeR::estimateCommonDisp(edgeR::DGEList(V[,2,]))$common.dispersion)
  names(omega) <- dimnames(V)[[2]]
  if (!is.null(Data$Bait)) {
    B <- if(Data$multiBait) Data$Bait[,1,] else Data$Bait[,1,drop=FALSE]
    S <- cbind(V[,1,], B)
    colnames(S) <- NULL
    omega[3] <- edgeR::estimateCommonDisp(edgeR::DGEList(S))$common.dispersion
    names(omega)[3] <- "baitEffect"
  }
  omega
}
