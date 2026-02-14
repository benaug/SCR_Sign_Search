getCellR <- function(u,res,cells,xlim,ylim){
  inout <- 1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
  if(inout==1){
    this.cell.init <- cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
  }else{
    this.cell.init <- 0
  }
  return(this.cell.init)
}

init.data <- function(data=NA,M=NA,inits=NA){
  data <- c(data$constants,data$capture,data$telemetry) #restructure data list
  J <- data$J
  K <- data$K
  xlim <- data$xlim
  ylim <- data$ylim
  cells <- data$cells
  dSS <- data$dSS
  InSS <- data$InSS

  #get some inits, actually sigma is all we need
  sigma <- inits$sigma
  
  n <- nrow(data$y)
  y <- matrix(0,M,J)
  y[1:n,] <- data$y
  #Initialize z, just using observed z's
  z.init <- c(rep(1,n),rep(0,M-n))
  s.init <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
  for(i in 1:M){
    if(sum(y[i,])>0){#if captured
      trapcaps <- which(y[i,]>0)
      if(length(trapcaps)>1){
        s.init[i,] <- colMeans(X[trapcaps,])
      }else{
        s.init[i,] <- X[trapcaps,]
      }
    }
  }
  #move any initialized outside state space
  for(i in 1:M){
    s.cell.init <- getCellR(s.init[i,],res,cells,xlim,ylim)
    if(InSS[s.cell.init]==0){#not in SS, move to nearest cell
      dists <- sqrt((dSS[s.cell.init,1]-dSS[,1])^2+(dSS[s.cell.init,2]-dSS[,2])^2)
      dists[InSS==0] <- Inf
      pick <- which(dists==min(dists))[1] #if more than 1, just use first
      s.init[i,] <- dSS[pick,]
    }
  }
  s.tel.init <- apply(data$u.tel,c(1,3),mean,na.rm=TRUE)
  #move any initialized outside state space
  for(i in 1:data$n.tel.inds){
    s.cell.init <- getCellR(s.tel.init[i,],res,cells,xlim,ylim)
    if(InSS[s.cell.init]==0){#not in SS, move to nearest cell
      dists <- sqrt((dSS[s.cell.init,1]-dSS[,1])^2+(dSS[s.cell.init,2]-dSS[,2])^2)
      dists[InSS==0] <- Inf
      pick <- which(dists==min(dists))[1] #if more than 1, just use first
      s.tel.init[i,] <- dSS[pick,]
    }
  }
  
  return(list(y=y,z=z.init,s=s.init,s.tel=s.tel.init,N=sum(z.init)))
}