# Manual ------------------------------------------------------------------

vn <- c("csv/ssVectorA-n.csv", "csv/ssVectorB-n.csv")
vs <- c("csv/ssVectorA-s.csv", "csv/ssVectorB-s.csv")
Vector <- readDeepn(vn, vs)
bn <- c("csv/diUb1-n.csv")
bs <- c("csv/diUb1-s.csv")
Data <- import(vn, vs, bn, bs)
Data <- applyFilter(Data, thresh=3)
Data <- applyFilter(Data, thresh=30)
dim(Data$Vector)
Data <- import(vn, vs, bn, bs)
Data <- applyFilter(Data, thresh=30)
dim(Data$Vector)


# Using config ------------------------------------------------------------

Data <- importFromDeepn("deepn1.config")
Data$omega <- overdisp(Data)
Data <- applyFilter(Data, thresh=200)
runMCMC(Data)

object <- psm()
head(summary(object))
head(summary(object, sort=TRUE))
summary(object, outfile="out.txt")

analyzeDeepn("deepn1.config", debug=TRUE)

# Sandbox -----------------------------------------------------------------

