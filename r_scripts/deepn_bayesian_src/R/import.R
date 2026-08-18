import <- function(vn, vs, bn, bs) {
  Vector <- readDeepn(vn, vs)
  if (!missing(bn)) {
    Bait <- readDeepn(bn, bs)
    multiBait <- length(dim(Bait))==3
    if (!all.equal(dimnames(Vector)[[1]], dimnames(Bait)[[1]])) stop("Genes don't match between vector and bait")
    btr <- if (multiBait) apply(Bait, 2:3, sum) else apply(Bait, 2, sum)
  } else {
    Bait <- multiBait <- btr <- NULL
  }
  Counts <- list(Vector=Vector, Bait=Bait)
  structure(list(Counts=Counts,
                 Vector=Vector,
                 Bait=Bait,
                 multiBait = multiBait,
                 vtr = apply(Vector, 2:3, sum),
                 btr = btr),
            class="data.deepn")
}
