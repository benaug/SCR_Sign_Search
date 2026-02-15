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
  log(lambda.detect[1:J]) <- beta0.lam + beta1.lam*E[1:J]
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
    #dynamic sparse matrix representation of cells with nonzero prob of availability to trim 
    #use dist and obsmod likelihood calculations below, only needs to be recomputed when s_i or sigma update
    #forces detection cells to be in pos.cells
    pos.cells[i,1:n.InSS.cells] <- getPosCells(avail.dist=avail.dist[i,1:n.cells],
                                               InSS.cells=InSS.cells[1:n.InSS.cells],
                                               det.cells=det.cells[i,1:max.det],n.det=n.det[i])
    n.pos.cells[i] <- getNPosCells(pos.cells=pos.cells[i,1:n.InSS.cells])
    #Individual use distributions - multiply rsf and available distribution, normalize. trimmed
    use.dist[i,1:n.cells] <- getUse(rsf=rsf[1:n.cells],avail=avail.dist[i,1:n.cells],
                                    pos.cells=pos.cells[i,1:n.InSS.cells],
                                    n.pos.cells=n.pos.cells[i],n.cells=n.cells)
    #observation model that is trimmed below
    # for(j in 1:J){
    #   # expected detections at detector is proportional to use of the cell containing detector
    #   lam[i,j] <- lambda.detect[j]*use.dist[i,detector.to.cell[j]] #detector.to.cell maps detectors to cells
    #   y[i,j] ~ dpois(lambda=lam[i,j]*z[i])
    # }
    #trimmed obsmod, only evaluate detectors in pos.cells
    y[i,1:J] ~ dObs(pos.cells=pos.cells[i,1:n.InSS.cells],n.pos.cells=n.pos.cells[i],
                    cell.to.detector=cell.to.detector[1:n.cells],
                    lambda.detect=lambda.detect[1:J],use.dist=use.dist[i,1:n.cells],z=z[i])
  }
  #Model for observed locations|detection (conditioned on cell of detection)
  for(i in 1:n){ #only detected individuals have observed locations
    for(l in 1:n.u.ind[i]){ #locations for this individual, possibly at different detectors
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
    use.dist.tel[i,1:n.cells] <- getUseTel(rsf=rsf[1:n.cells],avail=avail.dist.tel[i,1:n.cells])
    u.cell.tel[i,1:n.locs.ind[i]] ~ dCatVector(use.dist.tel[i,1:n.cells],n.locs.ind=n.locs.ind[i])
    u.tel[i,1:n.locs.ind[i],1:2] ~ dTruncNormVector(s=s.tel[i,1:2],sigma=sigma,n.locs.ind=n.locs.ind[i],
                                                    u.xlim=u.xlim.tel[i,1:n.locs.ind[i],1:2],
                                                    u.ylim=u.ylim.tel[i,1:n.locs.ind[i],1:2])
  }
})# end model
