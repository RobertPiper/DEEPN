importFromDeepn <- function(file) {
  raw <- readLines(file)
  config <- matrix(trimws(unlist(strsplit(raw, "="))), ncol=2, byrow=TRUE)
  key <- config[,1]
  val <- config[,2]
  ind <- order(key)
  key <- key[ind]
  val <- val[ind]
  vn <- val[grep("Vector_Non", key)]
  vs <- val[grep("Vector_Sel", key)]
  bn <- val[grep("Bait._Non", key)]
  bs <- val[grep("Bait._Sel", key)]
  Data <- import(vn, vs, bn, bs)
  applyFilter(Data, as.numeric(val[key=="Threshold"]))
}
