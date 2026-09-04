#cells and detectors should be set up so that detectors (X) are at grid centroids of the cell they are in
#this example is set up that way, be careful changing it.
library(coda)
library(viridisLite) #for plot colors
library(nimble) #data simulator uses nimble
library(truncnorm) #required for data simulator
source("sim.SCR.SignSearch multiK.R")
source("init.data multiK.R")
source("NimbleModel SCR Sign Search multiK.R")
source("Nimble Functions SCR Sign Search multiK.R") #nimble functions used in data simulator
source("sSampler Dcov RSF.R")

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

#state space. Must start at (0,0)
xlim <- c(0,75)
ylim <- c(0,75)
res <- 1.5 #resolution, cell width/height.
if(xlim[1]!=0|ylim[1]!=0)stop("xlim and ylim must start at 0.")
if((diff(range(xlim))/res)%%1!=0)stop("The range of xlim must be divisible by 'res'")
if((diff(range(ylim))/res)%%1!=0)stop("The range of ylim must be divisible by 'res'")

#make discrete state space objects
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res)
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res)
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#simulate a D.cov
library(fields)
set.seed(103543)
grid <- list(x=x.vals,y=y.vals) 
obj <- Exp.image.cov(grid=grid,aRange=20,setup=TRUE)
D.cov <- sim.rf(obj)
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

#simulate an rsf cov, lower cov.pars for finer scale cov
set.seed(24674345)
obj <- Exp.image.cov(grid=grid,aRange=5,setup=TRUE)
rsf.cov <- sim.rf(obj)
rsf.cov <- as.numeric(scale(rsf.cov)) #scale
image(x.vals,y.vals,matrix(rsf.cov,n.cells.x,n.cells.y),main="rsf.cov",xlab="X",ylab="Y",col=cols1)

#make state space mask - just making a circle
center <- colMeans(dSS)
dists <- sqrt((dSS[,1]-center[1])^2+(dSS[,2]-center[2])^2)
rad <- min(diff(xlim),diff(ylim))/2
InSS <- 1*(dists<rad)
tmp <- D.cov
tmp[InSS==0] <- -Inf
image(x.vals,y.vals,matrix(tmp,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

#make some detectors. Really, we should put them in transect formation
#but skipping that for now and spacing equally
sigma.target <- 3 #expected sigma for spacing
#X is only used for plotting and identifying cell numbers
X <- as.matrix(expand.grid(seq(9.75,65.25,2*sigma.target),seq(9.75,65.25,2*sigma.target)))
#for plots to look good, snap X to dSS
for(j in 1:nrow(X)){
  x.diffs <- abs(X[j,1]-x.vals)
  pick <- which(x.diffs==min(x.diffs))
  X[j,1] <- x.vals[pick]
  y.diffs <- abs(X[j,2]-y.vals)
  pick <- which(y.diffs==min(y.diffs))
  X[j,2] <- y.vals[pick]
}
#remove any outside state space
dists2 <- sqrt((X[,1]-center[1])^2+(X[,2]-center[2])^2)
rad2 <- rad - 2*sigma.target
X <- X[dists2<rad2,]
points(X,pch=4,lwd=2)
J <- nrow(X)
K <- 2 #number of detection occasions; this model version requires K>1

D.beta0 <- -4.25 #baseline D
D.beta1 <- 1.0 #density coefficient 
rsf.beta <- 0.5 #rsf coefficient
# beta0.lam <- 3.5 #detection rate intercept
beta0.lam <- rep(3.5,K) #occasion-specific detection rate intercept by occasion
beta1.lam <- 1 #effort coefficient
sigma <- 3 #spatial scale of availability distribution
n.tel.inds <- 5 #number of telemetry individuals
K.tel <- 10 #number of telemetry locations per individual

set.seed(32444) #change this for new data set

survey <- matrix(1,K,J) #1 if detector j operated on occasion k, 0 otherwise
E <- matrix(log(runif(K*J,0,1)),K,J) #occasion x detector log effort
lambda.detect <- survey*exp(beta0.lam + beta1.lam*E)
plot(lambda.detect[1,]~exp(E[1,])) #plot first occasion; function of linear effort

data <- sim.SCR.SignSearch(D.beta0=D.beta0,D.beta1=D.beta1,rsf.beta=rsf.beta,
               beta0.lam=beta0.lam,beta1.lam=beta1.lam,E=E,K=K,sigma=sigma,X=X,
               xlim=xlim,ylim=ylim,res=res,InSS=InSS,D.cov=D.cov,rsf.cov=rsf.cov,
               n.tel.inds=n.tel.inds,K.tel=K.tel,survey=survey)

data$truth$lambda.N #expected abundance from D cov inputs
data$truth$N #simulated realized abundance
data$truth$n #number of inds detected
# table(rowSums(data$capture$y>0)) #number of inds captures X times (events not counts)
table(apply(data$capture$y>0,1,sum)) #number of individual detector x occasion detection events

#cells and locations of observatons
str(data$capture$u.cell) #n x max(y)
str(data$capture$u) #n x max(y) x 2
#detector numbers of each location (for reference, using u.cell as data)
str(data$truth$this.j) #n x max(y)
#used in model instead of X
data$constants$cell.to.detector #maps cells to detectors

#crappy plot
tmp <- rep(0,n.cells)
image(x.vals,y.vals,matrix(tmp,length(x.vals),length(y.vals)),
      main="",xlab="X",ylab="Y",col="white")
grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80")
points(X,pch=4,lwd=2)
for(i in 1:data$capture$n){
  if(data$capture$n.u.ind[i]==1){
    s <- data$capture$u[i,1,]
  }else{
    s <- colMeans(data$capture$u[i,1:data$capture$n.u.ind[i],])
  }
  points(s[1],s[2],pch=16,col="darkred",cex=1)
  for(l in 1:data$capture$n.u.ind[i]){
    lines(x=c(s[1],data$capture$u[i,l,1]),
          y=c(s[2],data$capture$u[i,l,2]))
    points(data$capture$u[i,l,1],data$capture$u[i,l,2],pch=16,
           col="goldenrod",cex=0.75)
  }
}

##Fit Model
nimbleOptions(determinePredictiveNodesInModel = FALSE) #must run this line
M <- 250 #data augmentation limit. Must be larger than simulated N. If N posterior hits M, need to raise M and try again.
if(M<=data$truth$N)stop("Raise M to be larger than simulated N.")

inits <- list(sigma=4) #needs to be set somewhere in the ballpark of truth
nimbuild <- init.data(data=data,inits=inits,M=M)

Niminits <- list(z=nimbuild$z,N=nimbuild$N, #must init N to be sum(z.init)
                 s=nimbuild$s,
                 D0=sum(nimbuild$z)/(sum(InSS)*res^2),D.beta1=0,
                 sigma=inits$sigma,
                 # beta0.lam=0,beta1.lam=0,
                 beta0.lam=rep(0,K),beta1.lam=0,
                 rsf.beta=2,
                 s.tel=nimbuild$s.tel)

#constants for Nimble
#Precompute grid edges once; getAvail1D reuses them for every individual/update.
x.vals.edges <- c(data$constants$x.vals-data$constants$res/2,
                  data$constants$x.vals[data$constants$n.cells.x]+data$constants$res/2)
y.vals.edges <- c(data$constants$y.vals-data$constants$res/2,
                  data$constants$y.vals[data$constants$n.cells.y]+data$constants$res/2)
#standard-normal cutoff for trimming negligible movement availability outside +/- avail.z SD
avail.z <- qnorm(1-1e-10)
#Precompute detector x/y cell indices so dObs does no cell-index conversion during MCMC.
detector.cell.x <- data$constants$detector.to.cell%%data$constants$n.cells.x
detector.cell.y <- floor(data$constants$detector.to.cell/data$constants$n.cells.x)+1
detector.cell.y[detector.cell.x==0] <- detector.cell.y[detector.cell.x==0]-1
detector.cell.x[detector.cell.x==0] <- data$constants$n.cells.x
detector.cell.x <- as.integer(detector.cell.x)
detector.cell.y <- as.integer(detector.cell.y)
#using this to keep telemetry individauls inside the state space, using same one as detected individuals here,
#but could be different.
pi.cell.tel <- InSS/sum(InSS)
constants <- list(M=M,J=J,K=K,n.cells=data$constants$n.cells,n.cells.x=data$constants$n.cells.x,
                  n.cells.y=data$constants$n.cells.y,res=data$constants$res,
                  x.vals.edges=x.vals.edges,y.vals.edges=y.vals.edges,avail.z=avail.z,
                  detector.cell.x=detector.cell.x,detector.cell.y=detector.cell.y,
                  n.locs.ind=data$constants$n.locs.ind,n.tel.inds=data$constants$n.tel.inds,
                  D.cov=D.cov,cellArea=data$constants$res^2,
                  n=data$capture$n,n.u.ind=data$capture$n.u.ind) #indexing for detection locations inside cells

#supply data to nimble
Nimdata <- list(y=nimbuild$y,E=data$constants$E,survey=data$constants$survey,#log transformed if log transformed above
                u.tel=data$telemetry$u.tel,u.cell.tel=data$telemetry$u.cell.tel,
                InSS=data$constants$InSS,
                z=nimbuild$z,rsf.cov=rsf.cov,pi.cell.tel=pi.cell.tel,
                u.cell=data$capture$u.cell,u=data$capture$u)

# set parameters to monitor
parameters <- c('beta0.lam','beta1.lam','rsf.beta','D.beta1','sigma','N','D0','lambda.N')

#can also monitor a different set of parameters with a different thinning rate
nt <- 2 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#tell nimble which nodes to configure so we don't waste time for samplers we will replace below
config.nodes <- c('beta0.lam','beta1.lam','sigma','rsf.beta') #will add block updates for density parameters below
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,nodes=config.nodes) 

###*required* N/z sampler
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
#nodes used for update, calcNodes + z nodes
y.nodes <- Rmodel$expandNodeNames("y")
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
s.nodes <- Rmodel$expandNodeNames("s")
calcNodes <- c(N.node,y.nodes,s.nodes)
ind.detected <- 1*(apply(nimbuild$y,1,sum)>0)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,K=K,ind.detected=ind.detected,
                                                 y.nodes=y.nodes,N.node=N.node,z.nodes=z.nodes,s.nodes=s.nodes,
                                                 calcNodes=calcNodes,
                                                 res=data$constants$res,x.vals.edges=x.vals.edges,
                                                 y.vals.edges=y.vals.edges,avail.z=avail.z,
                                                 n.cells=data$constants$n.cells,n.cells.x=data$constants$n.cells.x,
                                                 n.cells.y=data$constants$n.cells.y),silent = TRUE)

#add sSampler
for(i in 1:M){
  calcNodes <- Rmodel$getDependencies(paste("s[",i,",1:2]"))
  conf$addSampler(target = paste("s[",i,",1:2]", sep=""),
                  type = 'sSamplerDcovRSF',control = list(i=i,xlim=data$constants$xlim,ylim=data$constants$ylim,
                                                          calcNodes=calcNodes), silent = TRUE)
}

#add efficient s sampler for telemetry inds
for(i in 1:data$constants$n.tel.inds){
  calcNodes.s.tel <- Rmodel$getDependencies(paste("s.tel[",i,",1:2]"))
  conf$addSampler(target = paste("s.tel[",i,",1:2]", sep=""),
                  type = 'sSamplerDcovRSF.tel',control = list(i=i,xlim=xlim,ylim=ylim,
                                                              calcNodes.s.tel=calcNodes.s.tel), silent = TRUE)
}

#AF slice is slower but tends to pay off for density covariates
conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)

#Add block RW update, if posteriors correlated. This was 2x as efficient as AF_slice in a few examples I looked at with 2 occasions
conf$addSampler(target = c("beta0.lam","beta1.lam"),
                type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2000,reset=FALSE) #can keep running this line to extend sampler
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[-c(1:burnin),])) #discarding some burnin here. Can't plot 1st sample which is all NA

data$truth$lambda.N #target expected abundance
data$truth$N #target realized abundance

#Check posterior correlation
rem.idx <- which(colnames(mvSamples)%in%c("N","lambda.N"))
tmp <- cor(mvSamples[-c(1:burnin),-rem.idx])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)
round(tmp,2)

