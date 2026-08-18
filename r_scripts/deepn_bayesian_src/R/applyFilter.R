# Need to make this work for when Data is just vector
applyFilter <- function(Data, thresh, base=TRUE) {
  RPM <- rpm(Data)

  # Baseline filter
  if (base) {
    if (Data$multiBait) {
      N <- cbind(apply(RPM$Vector[,1,], 1, mean), RPM$Bait[,1,])
    } else {
      N <- cbind(apply(RPM$Vector[,1,], 1, mean), RPM$Bait[,1])
    }
    bPass <- apply(N > thresh, 1, all)
  } else {
    bPass <- rep(TRUE, nrow(RPM$Vector))
  }

  # Selection filter
  if (Data$multiBait) {
    S <- cbind(RPM$Vector[,2,], RPM$Bait[,2,])
  } else {
    S <- cbind(RPM$Vector[,2,], RPM$Bait[,2])
  }
  sPass <- apply(S > thresh, 1, any)

  # Subset and return
  pass <- which(bPass & sPass)
  if (length(pass) <= 1) warning("Fewer than 2 genes pass this filter")
  Data$Vector <- Data$Counts$Vector[pass,,]
  Data$Bait <- if (Data$multiBait) Data$Counts$Bait[pass,,] else Data$Counts$Bait[pass,]
  Data
}
rpm <- function(Data) {
  Vector <- sweep(Data$Counts$Vector, 2:3, Data$vtr, "/")*1e6
  if (!is.null(Data$Bait)) {
    if (Data$multiBait) {
      Bait <- sweep(Data$Counts$Bait, 2:3, Data$btr, "/")*1e6
    } else {
      Bait <- sweep(Data$Counts$Bait, 2, Data$btr, "/")*1e6
    }
  } else{
    Bait <- NULL
  }
  list(Vector=Vector, Bait=Bait)
}
