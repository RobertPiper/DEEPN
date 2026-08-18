runMCMC <- function(Data, n.chains=4) {
  # Setup
  ng <- dim(Data$Vector)[1]
  nvr <- dim(Data$Vector)[3]
  monitor <- c('gamma')
  if (Data$multiBait) {
    nb <- dim(Data$Bait)[3]
    jData <- with(Data, list(v=Vector, y=Bait, ng=ng, nb=nb, nvr=nvr, ytr=btr, vtr=vtr, omega=omega[3:2]))
    model <- system.file("model/sModel2.jag", package="deepn")
    nsig <- 2 + nb
  } else {
    jData <- with(Data, list(v=Vector, y=Bait, ng=ng, nvr=nvr, ytr=btr, vtr=vtr, omega=omega[3:2]))
    model <- system.file("model/sModel1.jag", package="deepn")
    nsig <- 3
  }
  FUN <- function(chain) {list(log.sigma=runif(nsig))}
  inits <- mapply(FUN, 1:n.chains, SIMPLIFY=FALSE)

  runjags.options(rng.warning=FALSE)
  jagsfit <- run.jags(model, monitor=monitor, data=jData, inits=inits, n.chains=n.chains,
                      method="rjparallel", adapt=1500, burnin=1500, sample=1500, summarise=FALSE)
  list(jagsfit=jagsfit, Data=Data)
}
