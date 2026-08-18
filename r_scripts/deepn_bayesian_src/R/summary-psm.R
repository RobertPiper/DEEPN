summary.psm.deepn <- function(object, sort=FALSE, allGenes=FALSE, outfile, ...) {
  multi <- object$Data$multiBait
  sData <- summary(object$Data)

  # Summarize psm
  pm <- object$pm
  ps <- object$ps
  Z <- pm/ps
  if (multi) {
    sPSM <- data.frame(
      AdjEnr1 = pm[,1]/log(2),
      AdjEnr2 = pm[,2]/log(2),
      pBait1_Vec = pnorm(Z[,1]),
      pBait2_Vec = pnorm(Z[,2]),
      pBait1_Bait2 = pnorm(-Z[,3]),
      pBait2_Bait1 = pnorm(Z[,3]),
      pBait1 = pnorm(Z[,1])*pnorm(-Z[,3]),
      pBait2 = pnorm(Z[,2])*pnorm(Z[,3])
    )
  } else {
    sPSM <- data.frame(
      AdjEnr = pm/log(2),
      p = pnorm(Z)
    )
  }

  # Merge
  pass <- rownames(sData) %in% dimnames(object$Data$Vector)[[1]]
  if (allGenes) {
    sPSM_NULL <- sPSM[1,]
    sPSM_NULL[] <- NA
    Tab <- rbind(cbind(sData[pass,], sPSM),
                 cbind(sData[!pass,], sPSM_NULL))
  } else {
    Tab <- cbind(sData[pass,], sPSM)
  }

  # Sort
  m <- ncol(Tab)
  if (multi) {
    if (sort==1) {
      Tab <- Tab[order(-Tab[,m-1]),]
    } else if (sort==2) {
      Tab <- Tab[order(-Tab[,m]),]
    }
  } else {
    if (sort) Tab <- Tab[order(-Tab[,m]),]
  }

  # Output
  if (!missing(outfile)) {
    Tab <- cbind(Gene=rownames(Tab), Tab)
    write.csv(Tab, file=outfile, quote=FALSE, row.names=FALSE)
  } else {
    return(Tab)
  }
}
