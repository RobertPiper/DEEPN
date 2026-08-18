# displayGene <- function(ind) {
#   if (is.character(ind)) ind <- which(dimnames(Data)[[1]]==ind)
#   cat('\n')
#   cat(dimnames(Data)[[1]][ind], '\n')
#   cat('     Effect size: ', formatC(m[ind], digits=1, format='f'), '\n')
#   cat('     FDR:         ', formatC(q[ind], digits=2, format='g'), '\n')
#   print(Data[ind,,])
#   cat('\n')
# }
