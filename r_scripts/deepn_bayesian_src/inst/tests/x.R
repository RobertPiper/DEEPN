# Manual ------------------------------------------------------------------

vn <- c("csv/ssVectorA-n.csv", "csv/ssVectorB-n.csv")
vs <- c("csv/ssVectorA-s.csv", "csv/ssVectorB-s.csv")
Vector <- readDeepn(vn, vs)
bn <- c("csv/Rab5GDP-n.csv", "csv/Rab5GTP-n.csv")
bs <- c("csv/Rab5GDP-s.csv", "csv/Rab5GTP-s.csv")
Data <- import(vn, vs, bn, bs)
Data <- applyFilter(Data, thresh=3)
Data <- applyFilter(Data, thresh=30)
dim(Data$Vector)
Data <- import(vn, vs, bn, bs)
Data <- applyFilter(Data, thresh=30)
dim(Data$Vector)


# Using config ------------------------------------------------------------

Data <- importFromDeepn("deepn2.config")
Data$omega <- overdisp(Data)
Data <- applyFilter(Data, thresh=200)
runMCMC(Data)

object <- psm()
head(summary(object))
head(summary(object, sort=TRUE))
head(summary(object, sort=2))
summary(object, outfile="out.txt")

analyzeDeepn("deepn2.config", debug=TRUE)

analyzeDeepn("deepn2.config", debug=TRUE, port=1995)

# Sandbox -----------------------------------------------------------------

