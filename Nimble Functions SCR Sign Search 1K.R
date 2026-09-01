#Factored observation likelihood: evaluate use only at the J detector cells.
dObs <- nimbleFunction(
  run = function(x = double(1),detector.cell.x = double(1),detector.cell.y = double(1),
                 lambda.detect = double(1),rsf = double(1),avail.x = double(1),avail.y = double(1),
                 use.denom = double(0),n.cells.x = integer(0),z = integer(0),log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z==1){
      J <- nimDim(lambda.detect)[1]
      if(use.denom<=0){
        return(-Inf)
      }
      for(j in 1:J){
        cell <- detector.cell.x[j] + (detector.cell.y[j]-1)*n.cells.x
        use <- rsf[cell]*avail.x[detector.cell.x[j]]*avail.y[detector.cell.y[j]]/use.denom
        lam <- lambda.detect[j]*use
        logProb <- logProb + dpois(x[j],lam,log=TRUE)
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
  run = function(n = integer(0),detector.cell.x = double(1),detector.cell.y = double(1),
                 lambda.detect = double(1),rsf = double(1),avail.x = double(1),avail.y = double(1),
                 use.denom = double(0),n.cells.x = integer(0),z = integer(0)) {
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

#Factored replacement for getAvail in the fitted model. Cell availability is
#avail.x[cell.x]*avail.y[cell.y], so only the two 1-D vectors are stored.
getAvail1D <- nimbleFunction(
  run = function(s = double(0),sigma = double(0),res = double(0),vals.edges = double(1),
                 n.cells = integer(0),avail.z = double(0),z = integer(0)) {
    returnType(double(1))
    avail <- rep(0,n.cells)
    if(z==0){ #skip availability calculations for z=0 individuals
      return(avail)
    }
    lower <- s-sigma*avail.z
    upper <- s+sigma*avail.z
    if(vals.edges[1]<lower){
      idx.start <- floor((lower-vals.edges[1])/res)+1
    }else{
      idx.start <- 1
    }
    if(vals.edges[n.cells]>upper){
      idx.stop <- ceiling((upper-vals.edges[1])/res)
    }else{
      idx.stop <- n.cells
    }
    pnorm.vals <- rep(0,n.cells+1)
    for(l in idx.start:(idx.stop+1)){
      pnorm.vals[l] <- pnorm(vals.edges[l],mean=s,sd=sigma)
    }
    for(l in idx.start:idx.stop){
      avail[l] <- pnorm.vals[l+1]-pnorm.vals[l]
    }
    return(avail)
  }
)

#Scalar normalizing constant for RSF-weighted use. The scan identifies the
#nonzero x/y ranges so the nested loop does not use the full grid when sigma is small relative to grid extent.
getUseDenom <- nimbleFunction(
  run = function(rsf = double(1),avail.x = double(1),avail.y = double(1),
                 n.cells.x = integer(0),n.cells.y = integer(0),z = integer(0)) {
    returnType(double(0))
    if(z==0){ #skip use-denominator calculations for z=0 individuals
      return(0)
    }
    x.start <- n.cells.x
    x.stop <- 0
    y.start <- n.cells.y
    y.stop <- 0
    for(i in 1:n.cells.x){
      if(avail.x[i]>0){
        if(i<x.start){
          x.start <- i
        }
        if(i>x.stop){
          x.stop <- i
        }
      }
    }
    for(j in 1:n.cells.y){
      if(avail.y[j]>0){
        if(j<y.start){
          y.start <- j
        }
        if(j>y.stop){
          y.stop <- j
        }
      }
    }
    use.denom <- 0
    if(x.stop>0 & y.stop>0){
      for(j in y.start:y.stop){
        for(i in x.start:x.stop){
          cell <- i+(j-1)*n.cells.x
          use.denom <- use.denom + rsf[cell]*avail.x[i]*avail.y[j]
        }
      }
    }
    return(use.denom)
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

#Within-cell continuous likelihood using cached 1-D availability probabilities.
duInCell <- nimbleFunction(
  run = function(x = double(1),s = double(1),u.cell = double(0),sigma = double(0),
                 n.cells.x = integer(0),res = double(0),avail.x = double(1),avail.y = double(1),
                 log = integer(0)) {
    returnType(double(0))
    if(u.cell>0){
      u.cell.x <- u.cell%%n.cells.x
      u.cell.y <- floor(u.cell/n.cells.x)+1
      if(u.cell.x==0){
        u.cell.x <- n.cells.x
        u.cell.y <- u.cell.y-1
      }
      if(avail.x[u.cell.x]<=0 | avail.y[u.cell.y]<=0){
        return(-Inf)
      }
      logProb <- dnorm(x[1],s[1],sigma,log=TRUE)-log(avail.x[u.cell.x]) +
        dnorm(x[2],s[2],sigma,log=TRUE)-log(avail.y[u.cell.y])
    }else{
      logProb <- 0
    }
    return(logProb)
  }
)

#dummy RNG to make nimble happy
ruInCell <- nimbleFunction(
  run = function(n = integer(0),s = double(1),u.cell = double(0),sigma = double(0),
                 n.cells.x = integer(0),res = double(0),avail.x = double(1),avail.y = double(1)) {
    returnType(double(1))
    return(c(0,0))
  })

#Telemetry likelihood with the categorical cell probability and within-cell
#truncated Normal combined. The availability cell masses cancel exactly.
dTelemetryFactored <- nimbleFunction(
  run = function(x = double(2),u.cell = double(1),s = double(1),sigma = double(0),
                 rsf = double(1),use.denom = double(0),n.cells.x = integer(0),
                 n.locs.ind = double(0),log = integer(0)) {
    returnType(double(0))
    if(use.denom<=0){
      return(-Inf)
    }
    logProb <- 0
    for(l in 1:n.locs.ind){
      this.cell <- u.cell[l]
      if(rsf[this.cell]<=0){
        return(-Inf)
      }
      logProb <- logProb + log(rsf[this.cell])-log(use.denom) +
        dnorm(x[l,1],s[1],sigma,log=TRUE) + dnorm(x[l,2],s[2],sigma,log=TRUE)
    }
    return(logProb)
  }
)

rTelemetryFactored <- nimbleFunction(
  run = function(n = integer(0),u.cell = double(1),s = double(1),sigma = double(0),
                 rsf = double(1),use.denom = double(0),n.cells.x = integer(0),
                 n.locs.ind = double(0)) {
    returnType(double(2))
    return(matrix(0,n.locs.ind,2))
  }
)

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
    res <- control$res
    x.vals.edges <- control$x.vals.edges
    y.vals.edges <- control$y.vals.edges
    avail.z <- control$avail.z
    n.cells <- control$n.cells
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
  },
  run = function() {
    #Build undetected on/off lists once, then update after accepted proposals.
    z.on <- rep(0,M)
    z.off <- rep(0,M)
    non.curr <- 0
    noff.curr <- 0
    for(i in 1:M){
      if(ind.detected[i]==0){
        if(model$z[i]==1){
          non.curr <- non.curr+1
          z.on[non.curr] <- i
        }else{
          noff.curr <- noff.curr+1
          z.off[noff.curr] <- i
        }
      }
    }
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        non.init <- non.curr
        if(non.init>0){
          pick.pos <- rcat(1,rep(1/non.init,non.init))
          pick <- z.on[pick.pos]
          N.init <- model$N[1]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])

          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0

          #spatial process is gated by z; keep the off-state movement/use objects at zero
          model$avail.x[pick,1:n.cells.x] <<- rep(0,n.cells.x)
          model$avail.y[pick,1:n.cells.y] <<- rep(0,n.cells.y)
          model$use.denom[pick] <<- 0

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0

          #MH step
          #Restricting subtraction to undetected active individuals changes the cancellation
          #between P(z|N)=1/choose(M,N) and the individual-selection probabilities.
          #The remaining correction is non.init/N.init.
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) + log(non.init/N.init)
          accept <- decide(log_MH_ratio)

          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            mvSaved["avail.x",1][pick,1:n.cells.x] <<- model[["avail.x"]][pick,1:n.cells.x]
            mvSaved["avail.y",1][pick,1:n.cells.y] <<- model[["avail.y"]][pick,1:n.cells.y]
            mvSaved["use.denom",1][pick] <<- model[["use.denom"]][pick]
            #move guy from on list to off list
            z.on[pick.pos] <- z.on[non.curr]
            z.on[non.curr] <- 0
            non.curr <- non.curr-1
            noff.curr <- noff.curr+1
            z.off[noff.curr] <- pick
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model[["avail.x"]][pick,1:n.cells.x] <<- mvSaved["avail.x",1][pick,1:n.cells.x]
            model[["avail.y"]][pick,1:n.cells.y] <<- mvSaved["avail.y",1][pick,1:n.cells.y]
            model[["use.denom"]][pick] <<- mvSaved["use.denom",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          noff.init <- noff.curr
          if(noff.init>0){
            pick.pos <- rcat(1,rep(1/noff.init,noff.init))
            pick <- z.off[pick.pos]
            N.init <- model$N[1]
            #get initial logprobs for N and y
            lp.initial.N <- model$getLogProb(N.node)
            # lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
            lp.initial.y <- 0

            #propose new N/z
            model$N[1] <<-  model$N[1] + 1
            model$z[pick] <<- 1

            #construct the focal movement/use state directly when z is turned on, following JS z.super sampler
            model$avail.x[pick,1:n.cells.x] <<- getAvail1D(s=model$s[pick,1],sigma=model$sigma[1],res=res,
                                                           vals.edges=x.vals.edges,n.cells=n.cells.x,
                                                           avail.z=avail.z,z=1)
            model$avail.y[pick,1:n.cells.y] <<- getAvail1D(s=model$s[pick,2],sigma=model$sigma[1],res=res,
                                                           vals.edges=y.vals.edges,n.cells=n.cells.y,
                                                           avail.z=avail.z,z=1)
            model$use.denom[pick] <<- getUseDenom(rsf=model$rsf[1:n.cells],
                                                  avail.x=model$avail.x[pick,1:n.cells.x],
                                                  avail.y=model$avail.y[pick,1:n.cells.y],
                                                  n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)

            #get proposed logprobs for N and y
            lp.proposed.N <- model$calculate(N.node)
            lp.proposed.y <- model$calculate(y.nodes[pick])

            #MH step
            #Reverse subtraction selects only undetected active individuals. The remaining
            #combinatorial/proposal correction is (N.init+1)/(non.curr+1).
            log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) + log((N.init+1)/(non.curr+1))
            accept <- decide(log_MH_ratio)
            if(accept) {
              mvSaved["N",1][1] <<- model[["N"]]
              mvSaved["z",1][pick] <<- model[["z"]][pick]
              mvSaved["avail.x",1][pick,1:n.cells.x] <<- model[["avail.x"]][pick,1:n.cells.x]
              mvSaved["avail.y",1][pick,1:n.cells.y] <<- model[["avail.y"]][pick,1:n.cells.y]
              mvSaved["use.denom",1][pick] <<- model[["use.denom"]][pick]
              #move guy from off list to on list
              z.off[pick.pos] <- z.off[noff.curr]
              z.off[noff.curr] <- 0
              noff.curr <- noff.curr-1
              non.curr <- non.curr+1
              z.on[non.curr] <- pick
            }else{
              model[["N"]] <<- mvSaved["N",1][1]
              model[["z"]][pick] <<- mvSaved["z",1][pick]
              model[["avail.x"]][pick,1:n.cells.x] <<- mvSaved["avail.x",1][pick,1:n.cells.x]
              model[["avail.y"]][pick,1:n.cells.y] <<- mvSaved["avail.y",1][pick,1:n.cells.y]
              model[["use.denom"]][pick] <<- mvSaved["use.denom",1][pick]
              model$calculate(y.nodes[pick])
              model$calculate(N.node)
            }
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)