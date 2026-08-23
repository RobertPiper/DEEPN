overdisp <- function(Data) {
  # FIX: "Baseline" (non-selected) and "Selected" (Selection) dispersion are
  # supposed to characterize Vector on its own - but Data$Vector has already
  # been through applyFilter()'s JOINT vector+bait threshold filter by the
  # time overdisp() used to run, so the gene set feeding these two estimates
  # silently depended on whichever bait happened to be paired with Vector in
  # that run. Confirmed empirically: the same Vector12/Vector13 pair gave
  # "Selection" OD ranging 0.46-0.56 across 6 different bait pairings, purely
  # from this filtering artifact, not from any real difference in Vector's
  # own data. Fixed by computing these two from Data$Counts$Vector - the
  # original, full, unfiltered Vector counts that applyFilter() preserves
  # alongside its bait-filtered view - using a Vector-only threshold filter
  # (mirrors applyFilter()'s own filter logic, minus any bait columns).
  # "baitEffect" (omega[3]) is deliberately left using the bait-joint-filtered
  # view below - it exists specifically to compare Vector against a bait, so
  # bait-dependence there is correct, not a bug.
  Vraw <- Data$Counts$Vector
  thresh <- if (!is.null(Data$threshold)) Data$threshold else 0
  Vrpm <- sweep(Vraw, 2:3, Data$vtr, "/") * 1e6

  bPass_v <- apply(Vrpm[,1,] > thresh, 1, all)
  sPass_v <- apply(Vrpm[,2,] > thresh, 1, any)
  pass_v <- which(bPass_v & sPass_v)

  Vclean <- Vraw[pass_v,,]
  omega <- c(edgeR::estimateCommonDisp(edgeR::DGEList(Vclean[,1,]))$common.dispersion,
             edgeR::estimateCommonDisp(edgeR::DGEList(Vclean[,2,]))$common.dispersion)
  names(omega) <- dimnames(Data$Vector)[[2]]

  if (!is.null(Data$Bait)) {
    V <- Data$Vector
    B <- if(Data$multiBait) Data$Bait[,1,] else Data$Bait[,1,drop=FALSE]
    S <- cbind(V[,1,], B)
    colnames(S) <- NULL
    omega[3] <- edgeR::estimateCommonDisp(edgeR::DGEList(S))$common.dispersion
    names(omega)[3] <- "baitEffect"
  }
  omega
}
