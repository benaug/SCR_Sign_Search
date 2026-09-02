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
  for(k in 1:K){
    beta0.lam[k] ~ dnorm(0,sd=10)
  }
  beta1.lam ~ dnorm(0,sd=10) #effort coefficient
  
  #Density model
  D.intercept <- D0*cellArea
  # D.intercept <- exp(D.beta0)*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells]) 
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N)
  #Resource selection function evaluated across all cells in state space (InSS)
  rsf[1:n.cells] <- InSS[1:n.cells]*exp(rsf.beta*rsf.cov[1:n.cells])
  #detection rate a function of effort (log-transform E)
  for(k in 1:K){
    for(j in 1:J){
      lambda.detect[k,j] <- survey[k,j]*exp(beta0.lam[k] + beta1.lam*E[k,j])
    }
  }
  for(i in 1:M){
    # z-gated AC distribution. When z=0, s is fixed at c(0,0); when z=1, s is drawn
    # from pi.cell and then uniformly within the selected cell.
    s[i,1:2] ~ dAC(pi.cell=pi.cell[1:n.cells],res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=z[i])
    #Factored individual availability distribution. The BVN cell probability is
    #avail.x[cell.x]*avail.y[cell.y], no n.cells-length availability vector is stored.
    avail.x[i,1:n.cells.x] <- getAvail1D(s=s[i,1],sigma=sigma,res=res,vals.edges=x.vals.edges[1:(n.cells.x+1)],
                                         n.cells=n.cells.x,avail.z=avail.z,z=z[i])
    avail.y[i,1:n.cells.y] <- getAvail1D(s=s[i,2],sigma=sigma,res=res,vals.edges=y.vals.edges[1:(n.cells.y+1)],
                                         n.cells=n.cells.y,avail.z=avail.z,z=z[i])
    #normalizing constant for the RSF-weighted use distribution.
    #avail.x/avail.y can be reused when rsf.beta or detection parameters update.
    use.denom[i] <- getUseDenom(rsf=rsf[1:n.cells],avail.x=avail.x[i,1:n.cells.x],
                                avail.y=avail.y[i,1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=z[i])
    #Observation model
    for(k in 1:K){
      y[i,k,1:J] ~ dObs(detector.cell.x=detector.cell.x[1:J],detector.cell.y=detector.cell.y[1:J],
                        lambda.detect=lambda.detect[k,1:J],rsf=rsf[1:n.cells],
                        avail.x=avail.x[i,1:n.cells.x],avail.y=avail.y[i,1:n.cells.y],
                        use.denom=use.denom[i],n.cells.x=n.cells.x,z=z[i])
    }
  }
  #Model for observed locations|detection (conditioned on cell of detection)
  for(i in 1:n){ #only detected individuals have observed locations
    for(l in 1:n.u.ind[i]){ #locations for this individual, possibly at different detectors
      #u.cell[i,l] links to cell of detection
      #Reuse 1-D cell probabilities instead of recomputing normal CDF differences.
      u[i,l,1:2] ~ duInCell(s=s[i,1:2],u.cell=u.cell[i,l],sigma=sigma,n.cells.x=n.cells.x,res=res,
                            avail.x=avail.x[i,1:n.cells.x],avail.y=avail.y[i,1:n.cells.y])
    }
  }
  
  #optional telemetry
  for(i in 1:n.tel.inds){
    #use same density process as detected individuals
    # s.tel[i,1:2] ~ dAC(pi.cell=pi.cell[1:n.cells],res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)
    #uniform over InSS
    s.tel[i,1:2] ~ dAC(pi.cell=pi.cell.tel[1:n.cells],res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)
    #Factored telemetry availability and normalizing constants.
    avail.x.tel[i,1:n.cells.x] <- getAvail1D(s=s.tel[i,1],sigma=sigma,res=res,
                                             vals.edges=x.vals.edges[1:(n.cells.x+1)],
                                             n.cells=n.cells.x,avail.z=avail.z,z=1) #telemetry always active, so z=1
    avail.y.tel[i,1:n.cells.y] <- getAvail1D(s=s.tel[i,2],sigma=sigma,res=res,
                                             vals.edges=y.vals.edges[1:(n.cells.y+1)],
                                             n.cells=n.cells.y,avail.z=avail.z,z=1) #telemetry always active, so z=1
    use.denom.tel[i] <- getUseDenom(rsf=rsf[1:n.cells],avail.x=avail.x.tel[i,1:n.cells.x],
                                    avail.y=avail.y.tel[i,1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)
    #Cell-use probability and within-cell truncation cancel for telemetry locations.
    u.tel[i,1:n.locs.ind[i],1:2] ~ dTelemetryFactored(u.cell=u.cell.tel[i,1:n.locs.ind[i]],
                                                      s=s.tel[i,1:2],sigma=sigma,rsf=rsf[1:n.cells],
                                                      use.denom=use.denom.tel[i],n.cells.x=n.cells.x,
                                                      n.locs.ind=n.locs.ind[i])
  }
})# end model