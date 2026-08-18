readDeepn <- function(nFiles, sFiles) {
  if (length(nFiles)!=length(sFiles)) stop("Number of replicates under selection must equal number of repliates in non-selected conditions")
  K <- length(nFiles)
  for (k in 1:K) {
    for (j in 1:2) {
      # Read in raw
      if (j==1) fn <- nFiles[k] else fn <- sFiles[k]
      rawData <- read.csv(fn, colClasses=c(V1="character", V2="character", V3="numeric"), skip=4, header=FALSE, flush=TRUE)
      if (k==1 & j==1) {
        Gene <- gsub(' ', '', rawData$V2)
      } else {
        Gene.i <- gsub(' ', '', rawData$V2)
        if (any(Gene.i != Gene)) stop("PROBLEM: Genes do not match!")
      }

      # Convert to counts
      x <- rawData$V3
      y <- x / min(x[x!=0])
      if (max(abs(y - round(y))) > 0.01) stop("PROBLEM: Count rounding doesn't line up!")
      yy <- tapply(y, Gene, sum) # Combine across genes

      # Store in array
      if (k==1 & j==1) {
        Data <- array(NA, dim=c(length(yy), 2, K),
                      dimnames=list(names(yy), c("Baseline", "Selected"), LETTERS[1:K]))
      }
      Data[,j,k] <- yy
    }
  }
  storage.mode(Data) <- 'integer'
  drop(Data)
}
