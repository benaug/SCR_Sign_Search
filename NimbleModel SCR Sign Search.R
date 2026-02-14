NimModel <- nimbleCode({
  #Priors
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  #RSF coefficients
  rsf.beta ~ dnorm(0,sd=10) 
  #availability distribution spatial scale
  sigma ~ dunif(0,20)
  #detection intensity within searched cell a function of effort
  beta0.lam ~ dnorm(0,sd=10)
  beta1.lam ~ dnorm(0,sd=10) #effort coefficient
  
  #Density model
  D.intercept <- D0*cellArea
  # D.intercept <- exp(D.beta0)*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells]) 
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N)
  #Resource selection function evaluated across all cells
  rsf[1:n.cells] <- exp(rsf.beta*rsf.cov[1:n.cells])
  #detection rate a function of effort (log-transform E)
  for(j in 1:J){
    log(lambda.detect[j]) <- beta0.lam + beta1.lam*E[j]
  }
  #ACs, avail and use dists, detection model
  for(i in 1:M){
    #continuous activity center likelihood inside cell
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1] #extract activity center cell
    #categorical activity center likelihood for this cell, equivalent to zero's trick
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]])
    #Individual availability distributions conditioned on cells, bivariate Normal centered on activity center
    avail.dist[i,1:n.cells] <- getAvail(s=s[i,1:2],sigma=sigma,res=res,x.vals=x.vals[1:n.cells.x],
                                        y.vals=y.vals[1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y)
    #Individual use distributions - multiply rsf and available distribution, normalize
    use.dist[i,1:n.cells] <- getUse(rsf=rsf[1:n.cells],avail=avail.dist[i,1:n.cells])
    #extract ind by trap use probability from cells with traps, multiply by detection intensity
    for(j in 1:J){
      # expected detections at trap is proportional to use of the cell containing trap
      lam[i,j] <- lambda.detect[j]*use.dist[i,trap.to.cell[j]] #trap.to.cell maps traps to cells
    }
    y[i,1:J] ~ dPoissonVector(lambda=lam[i,1:J],z=z[i])
  }
  #Model for observed locations|detection (conditioned on cell of detection)
  for(i in 1:n){ #only detected individuals have observed locations
    for(l in 1:n.u.ind[i]){ #locations for this indivivdual, possibly at different traps
      #u.cell[i,l] links to cell of detection
      u[i,l,1:2] ~ duInCell(s=s[i,1:2],u.cell=u.cell[i,l],sigma=sigma,n.cells.x=n.cells.x,res=res)
    }
  }
  
  #optional telemetry
  for(i in 1:n.tel.inds){
    #can use telemetry data to inform D cov estimation, assumes inds captured at random wrt to response to covs
    #not linking telemetry individuals to D covs here, just uniform
    s.tel[i,1] ~ dunif(xlim[1],xlim[2])
    s.tel[i,2] ~ dunif(ylim[1],ylim[2])
    s.cell.tel[i] <- cells[trunc(s.tel[i,1]/res)+1,trunc(s.tel[i,2]/res)+1] #extract activity center cell
    dummy.data.tel[i] ~ dCell(pi.cell[s.cell.tel[i]])
    #Individual available distributions - bivariate Normal centered on activity center.
    avail.dist.tel[i,1:n.cells] <- getAvail(s=s.tel[i,1:2],sigma=sigma,res=res,x.vals=x.vals[1:n.cells.x],
                                            y.vals=y.vals[1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y)
    #Individual use distributions - multiply rsf and available distribution, normalize
    use.dist.tel[i,1:n.cells] <- getUse(rsf=rsf[1:n.cells],avail=avail.dist.tel[i,1:n.cells])
    u.cell.tel[i,1:n.locs.ind[i]] ~ dCatVector(use.dist.tel[i,1:n.cells],n.locs.ind=n.locs.ind[i])
    u.tel[i,1:n.locs.ind[i],1:2] ~ dTruncNormVector(s=s.tel[i,1:2],sigma=sigma,n.locs.ind=n.locs.ind[i],
                                                    u.xlim=u.xlim.tel[i,1:n.locs.ind[i],1:2],
                                                    u.ylim=u.ylim.tel[i,1:n.locs.ind[i],1:2])
  }
})# end model
