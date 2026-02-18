getPosCells <- nimbleFunction(
  run = function(avail.dist=double(1),InSS.cells=double(1),det.cells=double(1),n.det=integer(0)){
    returnType(double(1))
    n.InSS.cells <- nimDim(InSS.cells)[1]
    pos.cells <- rep(0,n.InSS.cells)
    idx <- 1
    #1) baseline trimming
    for(c in 1:n.InSS.cells){
      #sets level of trimming used to get use distribution and calculate y marginal logprobs
      if(avail.dist[InSS.cells[c]]>1e-10){ #not effectively zero
        pos.cells[idx] <- InSS.cells[c]
        idx <- idx + 1
      }
    }
    #2) force detection cells to be included
    if(n.det>0){
      for(l in 1:n.det){
        this.det <- det.cells[l]
        already <- 0
        for(c in 1:(idx-1)){
          if(pos.cells[c] == this.det){
            already <- 1
          }
        }
        if(already == 0){ #if not already in, stick it in first empty slot
          pos.cells[idx] <- this.det
          idx <- idx + 1
        }
      }
    }
    return(pos.cells)
  }
)

getNPosCells <- nimbleFunction(
  run = function(pos.cells=double(1)) {
    returnType(integer(0))
    n.pos.cells <- sum(pos.cells>0)
    return(n.pos.cells)
  }
)

dObs <- nimbleFunction(
  run = function(x = double(1),pos.cells = double(1),n.pos.cells = integer(0),
                 cell.to.detector = double(1),lambda.detect = double(1),
                 use.dist = double(1), z = integer(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z==1){
      for(c in 1:n.pos.cells){#loop over cells with nonnegligible use
        this.cell <- pos.cells[c]
        this.detector <- cell.to.detector[this.cell]
        if(this.detector>0){ #cell surveyed (is a detector cell)
          lam <- lambda.detect[this.detector]*use.dist[this.cell] #detection rate proportional to cell use
          logProb <- logProb + dpois(x[this.detector],lam,log=TRUE)
        }
      }
    }else{
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual (never occurs with all known IDs)
        logProb <- -Inf
      }
    }
    return(logProb)
  })

#make dummy random vector generator to make nimble happy
rObs <- nimbleFunction(
  run = function(n = integer(0), pos.cells = double(1),n.pos.cells = integer(0),
                 cell.to.detector = double(1),lambda.detect = double(1),
                 use.dist = double(1), z = integer(0)) {
    returnType(double(1))
    J <- nimDim(lambda.detect)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual (never occurs with all known IDs)
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lambda, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0), lambda = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(lambda)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

getCell <- nimbleFunction(#cell 0 not allowed in this model, but leaving in as an error check
  run = function(u = double(1),res=double(0),cells=integer(2),xlim=double(1),ylim=double(1)) {
    returnType(double(0))
    inout <- 1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
    if(inout==1){
      u.cell <- cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
    }else{
      u.cell <- 0
    }
    return(u.cell)
  }
)

getAvail <- nimbleFunction(
  run = function(s = double(1),sigma=double(0),res=double(0),x.vals=double(1),y.vals=double(1),n.cells.x=integer(0),n.cells.y=integer(0)) {
    returnType(double(1))
    avail.dist.x <- rep(0,n.cells.x)
    avail.dist.y <- rep(0,n.cells.y)
    delta <- 1e-10 #this sets the degree of trimming used to get individual availability distributions
    x.limits <- qnorm(c(delta,1-delta),mean=s[1],sd=sigma)
    y.limits <- qnorm(c(delta,1-delta),mean=s[2],sd=sigma)
    #convert to grid edges instead of centroids (could precompute)
    x.vals.edges <- c(x.vals - res/2, x.vals[n.cells.x] + res/2)
    y.vals.edges <- c(y.vals - res/2, y.vals[n.cells.y] + res/2)
    #trim in x and y direction
    if(x.vals.edges[1]<x.limits[1]){
      x.start <- floor((x.limits[1] - x.vals.edges[1]) / res) + 1
    }else{
      x.start <- 1
    }
    if(x.vals.edges[n.cells.x]>x.limits[2]){
      x.stop <- ceiling((x.limits[2] - x.vals.edges[1]) / res)
    }else{
      x.stop <- n.cells.x
    }
    if(y.vals.edges[1]<y.limits[1]){
      y.start <- floor((y.limits[1] - y.vals.edges[1]) / res) + 1
    }else{
      y.start <- 1
    }
    if(y.vals.edges[n.cells.y]>y.limits[2]){
      y.stop <- ceiling((y.limits[2] - y.vals.edges[1]) / res)
    }else{
      y.stop <- n.cells.y
    }
    pnorm.x <- rep(0,n.cells.x+1)
    pnorm.y <- rep(0,n.cells.y+1)
    #get pnorms
    for(l in x.start:(x.stop+1)){
      pnorm.x[l] <- pnorm(x.vals.edges[l],mean=s[1],sd=sigma)
    }
    for(l in y.start:(y.stop+1)){
      pnorm.y[l] <- pnorm(y.vals.edges[l],mean=s[2],sd=sigma)
    }
    for(l in (x.start):(x.stop)){
      avail.dist.x[l] <- pnorm.x[l+1] - pnorm.x[l]
    }
    for(l in (y.start):(y.stop)){
      avail.dist.y[l] <- pnorm.y[l+1] - pnorm.y[l]
    }
    avail.dist.tmp <- matrix(0,n.cells.x,n.cells.y)
    sum.dist <- 0
    for(i in x.start:x.stop){
      for(j in y.start:y.stop){
        avail.dist.tmp[i,j] <- avail.dist.x[i]*avail.dist.y[j]
        sum.dist <- sum.dist + avail.dist.tmp[i,j]
      }
    }
    avail.dist <- c(avail.dist.tmp)
    #if any probability mass is outside state space, normalize
    if(sum.dist<1){
      avail.dist <- avail.dist/sum.dist
    }
    return(avail.dist)
  }
)


# getUse <- nimbleFunction(
#   run = function(rsf = double(1),avail.dist=double(1)) {
#     returnType(double(1))
#     use.dist <- rsf*avail.dist
#     use.dist <- use.dist/sum(use.dist)
#     return(use.dist)
#   }
# )

getUse <- nimbleFunction(
  run = function(rsf = double(1),avail.dist=double(1),pos.cells=double(1),n.pos.cells=double(0),n.cells=double(0)){
    returnType(double(1))
    use.dist <- rep(0,n.cells)
    sum.dist <- 0
    for(c in 1:n.pos.cells){
      this.cell <- pos.cells[c]
      use.dist[this.cell] <- rsf[this.cell]*avail.dist[this.cell]
      sum.dist <- sum.dist + use.dist[this.cell]
    }
    for(c in 1:n.pos.cells){
      this.cell <- pos.cells[c]
      use.dist[this.cell] <- use.dist[this.cell]/sum.dist
    }
    return(use.dist)
  }
)

getUseTel <- nimbleFunction(
  run = function(rsf = double(1),avail.dist=double(1)) {
    returnType(double(1))
    use.dist <- rsf*avail.dist
    use.dist <- use.dist/sum(use.dist)
    return(use.dist)
  }
)

dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

#GPS vector distributions
dTruncNormVector <- nimbleFunction(
  run = function(x = double(2), s = double(1), sigma = double(0), u.xlim = double(2), u.ylim = double(2), n.locs.ind = double(0), log = integer(0)) {
    returnType(double(0))
    if(n.locs.ind>0){
      logProb <- 0 
      eps <- 1e-100 #prevent underflow 
      for(i in 1:n.locs.ind){ 
        # log density at point 
        logpx <- dnorm(x[i,1], mean = s[1], sd = sigma, log = TRUE)
        logpy <- dnorm(x[i,2], mean = s[2], sd = sigma, log = TRUE) # truncation probabilities for the cell bounds 
        denomx <- pnorm(u.xlim[i,2], mean = s[1], sd = sigma) - pnorm(u.xlim[i,1], mean = s[1], sd = sigma)
        denomy <- pnorm(u.ylim[i,2], mean = s[2], sd = sigma) - pnorm(u.ylim[i,1], mean = s[2], sd = sigma) 
        # if bounds are numerically zero/negative, reject 
        if(denomx <= 0.0 | denomy <= 0.0){ 
          return(-Inf) 
        } 
        # if extremely tiny, treat as essentially impossible 
        if(denomx < eps | denomy < eps){
          return(-Inf)
        } 
        logProb <- logProb + logpx - log(denomx) + logpy - log(denomy)
      } 
    }else{ 
      logProb <- 0 
    } 
    return(logProb) 
  })

rTruncNormVector <- nimbleFunction(
  run = function(n = integer(0), s = double(1), sigma = double(0), u.xlim = double(2), u.ylim = double(2), n.locs.ind = double(0)) {
    returnType(double(2))
    return(matrix(0,n.locs.ind,2))
  }
)

dCatVector <- nimbleFunction(
  run = function(x = double(1), use.dist = double(1), n.locs.ind = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    for(i in 1:n.locs.ind){
      logProb <- logProb + log(use.dist[x[i]])
    }
    return(logProb)
  }
)
rCatVector <- nimbleFunction(
  run = function(n = integer(0),use.dist = double(1), n.locs.ind = double(0)) {
    returnType(double(1))
    return(rep(0,n.locs.ind))
  }
)

#need this to add back in continuous locations
duInCell <- nimbleFunction(
  run = function(x = double(1), s = double(1), u.cell = double(0),
                 sigma = double(0),n.cells.x = integer(0),res=double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(u.cell>0){
      u.cell.x <- u.cell%%n.cells.x
      u.cell.y <- floor(u.cell/n.cells.x)+1
      if(u.cell.x==0){
        u.cell.x <- n.cells.x
        u.cell.y <- u.cell.y-1
      }
      xlim.cell <- c(u.cell.x-1,u.cell.x)*res
      ylim.cell <- c(u.cell.y-1,u.cell.y)*res
      logProb <- log(dnorm(x[1],s[1],sigma,log=FALSE)/
                       (pnorm(xlim.cell[2],s[1],sd=sigma) - pnorm(xlim.cell[1],s[1],sd=sigma))) + #x continuous likelihood
        log(dnorm(x[2],s[2],sigma,log=FALSE)/
              (pnorm(ylim.cell[2],s[2],sd=sigma) - pnorm(ylim.cell[1],s[2],sd=sigma))) #y continuous likelihood
    }else{
      logProb <- 0
    }
    return(logProb)
  }
)

#dummy RNG to make nimble happy
ruInCell <- nimbleFunction(
  run = function(n = integer(0), s = double(1), u.cell = double(0),
                 sigma = double(0),n.cells.x = integer(0),res=double(0)) {
    returnType(double(1))
    return(c(0,0))
  })

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    ind.detected <- control$ind.detected
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(ind.detected[pick]==1){#is this an individual with captured?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)