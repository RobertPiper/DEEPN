sqrtAxis <- function(ax=1:2) {
  for (i in ax) {
    axis(i, at=axTicks(i), labels=axTicks(i)^2, las=1)
  }
}
