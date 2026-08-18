summary.data.deepn <- function(object, sort=FALSE, outfile, ...) {
  Vector <- rpm(object)$Vector
  Bait <- rpm(object)$Bait
  multi <- object$multiBait
  Gene <- dimnames(Vector)[[1]]

  # Raw PPMs
  B <- cbind(Vector[,1,], if (multi) Bait[,1,] else Bait[,1])
  Base <- apply(B, 1, mean)
  Vec <- apply(Vector[,2,], 1, mean)
  if (multi) {
    Enr1 <- log2((Bait[,2,1]+0.05)/(Vec+0.05))
    Enr2 <- log2((Bait[,2,2]+0.05)/(Vec+0.05))
    Tab <- data.frame(
      Gene = Gene,
      Base = Base,
      Vec = Vec,
      Bait1 = Bait[,2,1],
      Bait2 = Bait[,2,2],
      Enr1 = Enr1,
      Enr2 = Enr2
    )
  } else {
    Enr <- log2((Bait[,2]+0.05)/(Vec+0.05))
    Tab <- data.frame(
      Gene = Gene,
      Base = Base,
      Vec = Vec,
      Bait = Bait[,2],
      Enr = Enr
    )
  }
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
  if (!missing(outfile)) {
    write.csv(Tab, file=outfile, quote=FALSE, row.names=FALSE)
  } else {
    return(Tab[,-1])
  }
}
