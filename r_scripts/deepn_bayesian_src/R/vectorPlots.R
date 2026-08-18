# vectorPlots <- function(Vector, name) {
#   fname <- paste0('web/img/', name, '-', 1:6, '.png')
#
#   PPM <- apply(Vector, 2:3, function(x) x/sum(x)*1e6)
#   m1 <- apply(PPM[,1,], 1, mean)
#   m2 <- apply(PPM[,2,], 1, mean)
#   v1 <- apply(PPM[,1,], 1, var)
#   v2 <- apply(PPM[,2,], 1, sd)
#   s1 <- sqrt(v1)
#   s2 <- sqrt(v2)
#   F <- PPM[,2,]/PPM[,1,]
#   mF <- apply(F, 1, mean)
#   vF <- apply(F, 1, var)
#   sF <- sqrt(vF)
#
#   #plot(m1, s2, pch=16, col=alphaCol('slateblue', .3), las=1, bty='n')
#   #plot(m1, sqrt(s2), pch=16, col=alphaCol('slateblue', .3), las=1, bty='n')
#   #plot(m1, s2, lwd=2, col=alphaCol('slateblue', .3), type='h', las=1, bty='n')
#   #plot(sqrt(m1), sqrt(s2), lwd=2, col=alphaCol('slateblue', .3), type='h', las=1, bty='n')
#   #plot(sqrt(m1), sqrt(v2/m1), lwd=2, col=alphaCol('slateblue', .3), type='h', las=1, bty='n')
#   #plot(sqrt(m1), sqrt(sF), lwd=2, col=alphaCol('slateblue', .3), type='h', las=1, bty='n', )
#
#   at.m <- seq(0, 40, 10)
#   lab.m <- at.m^2
#
#   png(fname[1], 6, 6, units='in', res=200)
#   plot(sqrt(m1), sqrt(s2), lwd=2, col=alphaCol('slateblue', .3), type='h', bty='n',
#        xlab='Mean (Unselected)', ylab='SD (Selected)', xaxt='n', yaxt='n')
#   sqrtAxis()
#   dev.off()
#
#   png(fname[2], 6, 6, units='in', res=200)
#   plot(sqrt(m1), sqrt(s2), lwd=2, col=alphaCol('slateblue', .3), pch=16, bty='n',
#        xlab='Mean (Unselected)', ylab='SD (Selected)', xaxt='n', yaxt='n')
#   sqrtAxis()
#   dev.off()
#
#   png(fname[3], 6, 6, units='in', res=200)
#   plot(sqrt(m1), sqrt(s2/m1), lwd=2, col=alphaCol('slateblue', .3), type='h', bty='n',
#        xlab='Mean (Unselected)', ylab='SD (Selected) / Mean (Unselected)', xaxt='n', yaxt='n')
#   sqrtAxis()
#   dev.off()
#
#   png(fname[4], 6, 6, units='in', res=200)
#   plot(sqrt(m1), sqrt(s2/m1), lwd=2, col=alphaCol('slateblue', .3), pch=16, bty='n',
#        xlab='Mean (Unselected)', ylab='SD (Selected) / Mean (Unselected)', xaxt='n', yaxt='n')
#   sqrtAxis()
#   dev.off()
#
#   png(fname[5], 6, 6, units='in', res=200)
#   plot(sqrt(m1), sqrt(sF), lwd=2, col=alphaCol('slateblue', .3), type='h', bty='n',
#        xlab='Mean (Unselected)', ylab='SD (Fold change)', xaxt='n', yaxt='n')
#   sqrtAxis()
#   dev.off()
#
#   png(fname[6], 6, 6, units='in', res=200)
#   plot(sqrt(m1), sqrt(sF), lwd=2, col=alphaCol('slateblue', .3), pch=16, bty='n',
#        xlab='Mean (Unselected)', ylab='SD (Fold change)', xaxt='n', yaxt='n')
#   sqrtAxis()
#   dev.off()
# }
# sqrtAxis <- function(ax=1:2) {
#   for (i in ax) {
#     axis(i, at=axTicks(i), lab=axTicks(i)^2, las=1)
#   }
# }
