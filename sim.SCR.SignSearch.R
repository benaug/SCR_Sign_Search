
sim.SCR.SignSearch <-
  function(D.beta0=NA,D.beta1=NA,rsf.beta=NA,beta0.lam=NA,beta1.lam=NA,E=NA,
           sigma=NA,X=X,
           xlim=NA,ylim=NA,res=NA,InSS=NA,D.cov=NA,rsf.cov=NA,effort=NA,survey=NA,
           K.tel=0,n.tel.inds=0){
    if(xlim[1]!=0|ylim[1]!=0)stop("xlim and ylim must start at 0.")
    if((diff(range(xlim))/res)%%1!=0)stop("The range of xlim must be divisible by 'res'")
    if((diff(range(ylim))/res)%%1!=0)stop("The range of ylim must be divisible by 'res'")
    # make discrete state space
    x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res)
    y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res)
    dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
    cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    cellArea <- res^2
    
    #Density Model
    D.intercept <- exp(D.beta0)*cellArea
    lambda.cell <- InSS*exp(D.beta1*D.cov)
    pi.denom <- sum(lambda.cell)
    pi.cell <- lambda.cell/pi.denom
    lambda.N <- D.intercept*pi.denom
    lambda.cell.plot <- lambda.cell
    lambda.cell.plot[lambda.cell==0] <- -Inf
    image(x.vals,y.vals,matrix(lambda.cell.plot,length(x.vals),length(y.vals)),
          main="Spatially Explicit D and Realized ACs",xlab="X",ylab="Y",col=rev(viridisLite::mako(100)))
    N <- rpois(1,lambda.N)

    #Activity centers
    # s.cell <- rcat(N,pi.cell)
    s.cell <- sample(1:n.cells,N,replace=TRUE,prob=pi.cell)
    s <- matrix(NA,N,2)
    for(i in 1:N){
      s.xlim <- dSS[s.cell[i],1] + c(-res,res)/2
      s.ylim <- dSS[s.cell[i],2] + c(-res,res)/2
      s[i,1] <- runif(1,s.xlim[1],s.xlim[2])
      s[i,2] <- runif(1,s.ylim[1],s.ylim[2])
    }
    points(s[,1],s[,2],pch=16)
    
    #RSF Model
    rsf <- exp(rsf.beta*rsf.cov)
    #get BVN availability and use distributions
    avail.dist <- use.dist <- matrix(NA,N,n.cells)
    for(i in 1:N){
      avail.dist[i,] <- getAvail(s=s[i,1:2],sigma=sigma,res=res,x.vals=x.vals,
                                 y.vals=y.vals,n.cells.x=n.cells.x,n.cells.y=n.cells.y)
      use.dist[i,] <- rsf*avail.dist[i,]
      use.dist[i,] <- use.dist[i,]/sum(use.dist[i,])
    }
    
    #map traps to cells
    J <- nrow(X)
    trap.to.cell <- rep(NA,J)
    for(j in 1:J){
      trap.to.cell[j] <- getCell(X[j,],res=res,cells=cells,xlim=xlim,ylim=ylim)
    }
    
    #simulate cell-level detection
    lambda.detect <- exp(beta0.lam + beta1.lam*E)
    y <- matrix(NA,N,J)
    for(i in 1:N){
      for(j in 1:J){
        lam <- lambda.detect[j]*use.dist[i,trap.to.cell[j]] #use.dist evaluated at cell trap j is in
        y[i,j] <- rpois(1,lam)
      }
    }
    #simulate locations of detections
    n.u.ind <- rowSums(y)
    n.u <- sum(n.u.ind)
    max.n.u <- max(n.u.ind)
    u <- array(NA,dim=c(N,max.n.u,2))
    this.j <- u.cell <- matrix(NA,N,max.n.u)
    for(i in 1:N){
      if(n.u.ind[i]>0){
        this.j[i,1:n.u.ind[i]] <- rep(which(y[i,]>0),times=y[i,which(y[i,]>0)])
        for(l in 1:n.u.ind[i]){
          u.cell[i,l] <- trap.to.cell[this.j[i,l]] #cell this trap is in
          u.xlim <- dSS[u.cell[i,l],1] + c(-res,res)/2
          u.ylim <- dSS[u.cell[i,l],2] + c(-res,res)/2
          u[i,l,1] <- rtruncnorm(1,a=u.xlim[1],b=u.xlim[2],
                                 mean=s[i,1],sd=sigma)
          u[i,l,2] <- rtruncnorm(1,a=u.ylim[1],b=u.ylim[2],
                                 mean=s[i,2],sd=sigma)
        }
      }
    }
    
    image(x.vals,y.vals,matrix(colSums(use.dist),length(x.vals),length(y.vals)),
          main="Capture Plot with Intensity of Use Summed over Individuals",xlab="X",ylab="Y",col=rev(viridisLite::mako(100)))
    grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80")
    detected <- which(rowSums(y)>0)
    points(s[,1],s[,2],pch=16,col="darkred",cex=1)
    y2D <- apply(y,c(1,2),sum)
    for(i in detected){
      for(l in 1:n.u.ind[i]){
        lines(x=c(s[i,1],u[i,l,1]),
              y=c(s[i,2],u[i,l,2]),col="darkred")
        points(u[i,l,1],u[i,l,2],pch=16,col="goldenrod",cex=0.75)
      }
    }
    points(X,pch=4,lwd=2)
    
    
    ####Telemetry data#####
    #simulating from same D model here
    if(n.tel.inds>0){
      s.tel.cell <- rcat(n.tel.inds,pi.cell)
      s.tel <- matrix(NA,n.tel.inds,2)
      for(i in 1:n.tel.inds){
        s.xlim <- dSS[s.tel.cell[i],1] + c(-res,res)/2
        s.ylim <- dSS[s.tel.cell[i],2] + c(-res,res)/2
        s.tel[i,1] <- runif(1,s.xlim[1],s.xlim[2])
        s.tel[i,2] <- runif(1,s.ylim[1],s.ylim[2])
      }
      #get BVN availability and use distributions
      avail.dist.tel <- use.dist.tel <- matrix(NA,N,n.cells)
      for(i in 1:n.tel.inds){
        avail.dist.tel[i,] <- getAvail(s=s.tel[i,1:2],sigma=sigma,res=res,x.vals=x.vals,
                                       y.vals=y.vals,n.cells.x=n.cells.x,n.cells.y=n.cells.y)
        use.dist.tel[i,] <- rsf*avail.dist.tel[i,]
        use.dist.tel[i,] <- use.dist.tel[i,]/sum(use.dist.tel[i,])
      }
      #simulate u's
      u.tel <- array(NA,dim=c(n.tel.inds,K.tel,2))
      u.cell.tel <- matrix(NA,n.tel.inds,K.tel)
      u.xlim.tel <- u.ylim.tel <- array(NA,dim=c(n.tel.inds,K.tel,2))
      for(i in 1:n.tel.inds){
        for(k in 1:K.tel){
          # u.cell.tel[i,k] <- rcat(1,use.dist.tel[i,])
          u.cell.tel[i,k] <- sample(1:n.cells,1,replace=TRUE,prob=use.dist.tel[i,])
          u.xlim.tel[i,k,] <- dSS[u.cell.tel[i,k],1] + c(-res,res)/2
          u.ylim.tel[i,k,] <- dSS[u.cell.tel[i,k],2] + c(-res,res)/2
          u.tel[i,k,1] <- rtruncnorm(1,a=u.xlim.tel[i,k,1],b=u.xlim.tel[i,k,2],
                                     mean=s.tel[i,1],sd=sigma)
          u.tel[i,k,2] <- rtruncnorm(1,a=u.ylim.tel[i,k,1],b=u.ylim.tel[i,k,2],
                                     mean=s.tel[i,2],sd=sigma)
        }
      }
      image(x.vals,y.vals,matrix(rsf.cov,length(x.vals),length(y.vals)),
            main="Telemetry Individual ACs and Site Use with rsf.cov",xlab="X",ylab="Y",
            col=rev(viridisLite::mako(100)))
      grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80")
      for(i in 1:n.tel.inds){
        for(k in 1:K.tel){
          lines(x=c(s.tel[i,1],u.tel[i,k,1]),y=c(s.tel[i,2],u.tel[i,k,2]),col="darkred")
        }
      }
      points(u.tel[,,1],u.tel[,,2],pch=16,cex=1,col="goldenrod")
      points(s.tel[,1],s.tel[,2],pch=16,cex=1.25,col="darkred")
      
      n.locs.ind <- rep(K.tel,n.tel.inds)
      n.locs.max <- max(n.locs.ind)
    }else{
      s.tel.cell <- NA
      s.tel <- NA
      u.tel <- NA
      u.cell.tel <- NA
      n.locs.ind <- NA
      n.locs.max <- NA
    }

    #discard uncaptured inds and disaggregate data
    caught <- which(rowSums(y)>0)
    n <- length(caught)
    y <- y[caught,]
    s <- s[caught,]
    s.cell <- s.cell[caught]
    
    #discard individuals with no detections from location data
    u.obs <- u[caught,,]
    this.j.obs <- this.j[caught,]
    u.cell.obs <- u.cell[caught,]
    n.u.ind.obs <- n.u.ind[caught]
    
    constants <- list(X=X,J=J,K.tel=K.tel,xlim=xlim,ylim=ylim,dSS=dSS,res=res,cells=cells,x.vals=x.vals,y.vals=y.vals,
                      n.tel.inds=n.tel.inds,n.locs.ind=n.locs.ind,trap.to.cell=trap.to.cell,
                      InSS=InSS,n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,D.cov=D.cov,rsf.cov=rsf.cov,
                      E=E)
    truth <- list(lambda.N=lambda.N,lambda.cell=lambda.cell,rsf=rsf,N=N,s=s,s.cell=s.cell,
                  n=n,s.tel=s.tel,s.tel.cell=s.tel.cell,
                  use.dist=use.dist,avail.dist=avail.dist,
                  u=u,this.j=this.j,u.cell=u.cell,n.u.ind=n.u.ind)
    capture <- list(y=y,u=u.obs,this.j=this.j.obs,u.cell=u.cell.obs,n.u.ind=n.u.ind.obs,n=n)
                    
    if(n.tel.inds>0){
      telemetry <- list(u.tel=u.tel,u.cell.tel=u.cell.tel,u.xlim.tel=u.xlim.tel,
                        u.ylim.tel=u.ylim.tel) #observed telemetry data
    }else{
      telemetry <- NA
    }
    out <- list(constants=constants,truth=truth,capture=capture,telemetry=telemetry)
    return(out)
  }
