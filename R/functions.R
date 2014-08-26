################################################################################
### Model data structures

emptyParameters <- function(num.bases,num.clusters,num.covariates) {
  probs <- rep(NA,num.clusters)
  coefs <- matrix(NA,num.bases,num.clusters)
  loadings <- matrix(NA,2,num.covariates)
  list(coefs=coefs,probs=probs,loadings=loadings)
}

emptyHypers <- function(num.bases, num.covariates) {
  alpha <- 2
  mu.coefs <- rep(0,num.bases)
  sigma.coefs <- 1e5*diag(num.bases)
  mu.loadings <- rep(0,num.covariates)
  sigma.loadings <- 1e5*diag(num.covariates)
  sigma.offset <- 1e-5*diag(2)
  amp <- 1
  len <- 1
  sig <- 1
  list(alpha=alpha,mu.coefs=mu.coefs,sigma.coefs=sigma.coefs,
       mu.loadings=mu.loadings,sigma.loadings=sigma.loadings,
       sigma.offset=sigma.offset,amp=amp,len=len,sig=sig)
}

emptyModel <- function(t.range,num.bases,num.clusters,num.covariates) {
  parameters <- emptyParameters(num.bases,num.clusters,num.covariates)
  hypers <- emptyHypers(num.bases,num.covariates)
  list(t.range=t.range,num.bases=num.bases,num.clusters=num.clusters,
       num.covariates=num.covariates,parameters=parameters,hypers=hypers)
}

#' @import MCMCpack
initModel <- function(model,data,z) {
  if (missing(z)) {
    z <- sample(model$num.clusters,length(data),TRUE)
  }
  model$parameters$probs[] <- rdirichlet(1,rep(1,model$num.clusters))
  model$parameters$loadings[] <- rnorm(2*model$num.covariates,sd=0.01)
  model$parameters$coefs[] <- 0
  for (i in 1:model$num.clusters) {
    t <- do.call('c',lapply(data[z==i],times))
    y <- do.call('c',lapply(data[z==i],values))
    X <- bSplineDesign(t,model$t.range,model$num.bases)
    model$parameters$coefs[,i] <- coef(lm(y~X-1))
  }
  model
}

################################################################################
### Sequence data structures

sqnc <- function(id,t,y,x) {
  s <- list(id=id[1],t=t,y=y,x=x)
  class(s) <- c('sqnc', 'list')
  s
}

id <- function(s) s$id
times <- function(s) s$t
values <- function(s) s$y
covariates <- function(s) s$x

sequenceMean <- function(s,z,model) {
  t <- times(s)
  x <- covariates(s)
  beta <- model$parameters$coefs[,z]
  m1 <- bSplineDesign(t,model$t.range,model$num.bases) %*% beta
  m2 <- linearDesign(t,model$t.range) %*% (model$parameters$loadings %*% x)
  as.numeric(m1 + m2)
}

sequenceCovariance <- function(s,model) {
  t <- times(s)
  amp <- model$hypers$amp
  len <- model$hypers$len
  sig <- model$hypers$sig
  v1 <- sqExpCov(t,amp,len,sig)
  ld <- linearDesign(t,model$t.range)
  v2 <- ld %*% model$hypers$sigma.offset %*% t(ld)
  v1 + v2
}

sequencePredictions <- function(s, model, new.t) {
  if (missing(new.t)) {
    new.t <- times(s)
  }
  t <- times(s)
  y <- values(s)
  x <- covariates(s)

  beta <- model$parameters$coefs
  amp <- model$hypers$amp
  len <- model$hypers$len
  sig <- model$hypers$sig

  inf <- sequencePosterior(s, model)
  z <- which.max(inf$q)
  y1 <- bSplineDesign(t, model$t.range, model$num.bases) %*% beta[, z]
  new.y1 <- bSplineDesign(new.t, model$t.range, model$num.bases) %*% beta[, z]

  A <- linearDesign(t, model$t.range)
  y.sigma <- diag(sig^2, length(y))
  b.mu <- model$parameters$loadings %*% x
  b.sigma <- model$hypers$sigma.offset
  b.sigma.post <- solve(solve(b.sigma) + t(A) %*% solve(y.sigma, A))
  b.mu.post <- b.sigma.post %*% (t(A) %*% solve(y.sigma, y-y1) + solve(b.sigma, b.mu))
  y2 <- A %*% b.mu.post
  new.y2 <- linearDesign(new.t, model$t.range) %*% b.mu.post

  K <- sqExpGram(new.t,t,amp,len)
  S <- sqExpCov(t,amp,len,sig)
  new.y3 <- K %*% solve(S,y-y1-y2)

  cbind(new.y1, new.y2, new.y3, new.y1+new.y2+new.y3)
}

################################################################################
### Inference functions

#' @import mvtnorm
sequencePosterior <- function(s,model) {
  y <- values(s)
  lp <- log(model$parameters$probs)
  sigma <- sequenceCovariance(s,model)
  for(i in seq(along=lp)) {
    mu <- sequenceMean(s,i,model)
    ll <- dmvnorm(y,mu,sigma,log=TRUE)
    lp[i] <- lp[i] + ll
  }
  logl <- logsumexp(lp)
  q <- exp(lp - logl)
  list(q=q,logl=logl)
}

################################################################################
### Estimation functions

eStep <- function(data,model) {
  q <- matrix(NA,model$num.clusters,length(data))
  logl <- 0
  for (i in seq(along=data)) {
    inf <- sequencePosterior(data[[i]],model)
    q[,i] <- inf$q
    logl <- logl + inf$logl
  }
  list(q=q,logl=logl)
}

#' @import mvtnorm
mStep <- function(data,q,model,eps=1e-4,maxiter=1e2) {
  emptyCoefSuffStats <- function(model) {
    b <- model$num.bases
    g <- model$num.clusters
    eta1 <- array(0,c(b,b,g))
    eta2 <- matrix(0,b,g)
    for (i in 1:g) {
      eta1[,,i] <- solve(model$hypers$sigma.coefs)
      eta2[,i] <- solve(model$hypers$sigma.coefs,model$hypers$mu.coefs)
    }
    list(eta1=eta1,eta2=eta2)
  }

  emptyLoadingSuffStats <- function(model) {
    nc <- model$num.covariates
    eta1 <- array(0,c(nc,nc,2))
    eta2 <- matrix(0,nc,2)
    for (i in 1:2) {
      eta1[,,i] <- solve(model$hypers$sigma.loadings)
      eta2[,i] <- solve(model$hypers$sigma.loadings,model$hypers$mu.loadings)
    }
    list(eta1=eta1,eta2=eta2)
  }

  mObjective <- function(model) {
    coefs <- model$parameters$coefs
    loadings <- model$parameters$loadings
    mu.coefs <- model$hypers$mu.coefs
    sigma.coefs <- model$hypers$sigma.coefs
    mu.loadings <- model$hypers$mu.loadings
    sigma.loadings <- model$hypers$sigma.loadings

    obj <- 0
    obj <- obj + dmvnorm(loadings[1,],mu.loadings,sigma.loadings)
    obj <- obj + dmvnorm(loadings[2,],mu.loadings,sigma.loadings)
    for (i in 1:ncol(coefs)) {
      obj <- obj + dmvnorm(coefs[,i],mu.coefs,sigma.coefs)
    }
    for (i in seq(along=data)) {
      s <- data[[i]]
      y <- values(s)
      v <- sequenceCovariance(data[[i]],model)
      for (j in seq(len=model$num.clusters)) {
        m <- sequenceMean(data[[i]],j,model)
        obj <- obj + q[j,i]*dmvnorm(y,m,v,log=TRUE)
      }
    }
    obj
  }

  orig.model <- model
  orig.obj <- mObjective(orig.model)
  model$parameters$probs <- rowSums(q)/sum(q)
  obj <- orig.obj
  convergence <- Inf
  iter <- 0

  while (convergence > eps && (iter <- iter + 1) <= maxiter) {
    coef.ss <- emptyCoefSuffStats(model)
    for (i in seq(along=data)) {
      s <- data[[i]]
      t <- times(s)
      y <- values(s)
      x <- covariates(s)
      d1 <- bSplineDesign(t,model$t.range,model$num.bases)
      d2 <- linearDesign(t,model$t.range)
      cv <- sequenceCovariance(s,model)
      eta1 <- t(d1) %*% solve(cv,d1)
      r <- y - d2 %*% model$parameters$loadings %*% x
      eta2 <- t(d1) %*% solve(cv,r)
      for (j in seq(len=model$num.clusters)) {
        coef.ss$eta1[,,j] <- coef.ss$eta1[,,j] + q[j,i]*eta1
        coef.ss$eta2[,j] <- coef.ss$eta2[,j] + q[j,i]*eta2
      }
    }
    for (i in seq(len=model$num.clusters)) {
      eta1 <- coef.ss$eta1[,,i]
      eta2 <- coef.ss$eta2[,i]
      model$parameters$coefs[,i] <- solve(eta1,eta2)
    }

    loading.ss <- emptyLoadingSuffStats(model)
    for (i in seq(along=data)) {
      s <- data[[i]]
      t <- times(s)
      y <- values(s)
      x <- covariates(s)
      d1 <- bSplineDesign(t,model$t.range,model$num.bases)
      d2 <- linearDesign(t,model$t.range)
      cv <- sequenceCovariance(s,model)
      eta1 <- outer(x,d2[,1]) %*% solve(cv,outer(d2[,1],x))
      ym <- (d1 %*% model$parameters$coefs) %*% q[,i]
      r <- y - ym - d2[,2] * (model$parameters$loadings[2,] %*% x)
      eta2 <- outer(x,d2[,1]) %*% solve(cv,r)
      loading.ss$eta1[,,1] <- loading.ss$eta1[,,1] + eta1
      loading.ss$eta2[,1] <- loading.ss$eta2[,1] + eta2
    }
    eta1 <- loading.ss$eta1[,,1]
    eta2 <- loading.ss$eta2[,1]
    model$parameters$loadings[1,] <- solve(eta1,eta2)

    for (i in seq(along=data)) {
      s <- data[[i]]
      t <- times(s)
      y <- values(s)
      x <- covariates(s)
      d1 <- bSplineDesign(t,model$t.range,model$num.bases)
      d2 <- linearDesign(t,model$t.range)
      cv <- sequenceCovariance(s,model)
      eta1 <- outer(x,d2[,2]) %*% solve(cv,outer(d2[,2],x))
      ym <- (d1 %*% model$parameters$coefs) %*% q[,i]
      r <- y - ym - d2[,1] * (model$parameters$loadings[1,] %*% x)
      eta2 <- outer(x,d2[,2]) %*% solve(cv,r)
      loading.ss$eta1[,,2] <- loading.ss$eta1[,,2] + eta1
      loading.ss$eta2[,2] <- loading.ss$eta2[,2] + eta2
    }
    eta1 <- loading.ss$eta1[,,2]
    eta2 <- loading.ss$eta2[,2]
    model$parameters$loadings[2,] <- solve(eta1,eta2)

    obj.old <- obj
    obj <- mObjective(model)
    convergence <- (obj-obj.old)/abs(orig.obj)
  }
  model
}

runEm <- function(data,model,tol=1e-4,maxiter=1e2) {
  for (iter in 1:maxiter) {
    inf <- eStep(data,model)
    model <- mStep(data,inf$q,model)

    if (iter == 1) {
      convergence <- Inf
    } else {
      if (old.logl > inf$logl) break
      convergence <- abs(inf$logl - old.logl)/abs(old.logl)
    }

    message(sprintf('iter=%03d LL=%.02f Convergence=%.04f',
                    iter,inf$logl,convergence))
    if (convergence < tol) break
    old.logl <- inf$logl
  }

  model$clusters <- apply(inf$q,2,which.max)
  model$logl <- inf$logl
  model
}

################################################################################
### Visualization functions

clusterMeans <- function(model,n.points=50) {
  t <- seq(model$t.range[1],model$t.range[2],length.out=n.points)
  X <- bSplineDesign(t,model$t.range,model$num.bases)
  y <- X %*% model$parameters$coefs
  list(t=t,y=y)
}

sequencePoints <- function(data) {
  t <- do.call('c',lapply(data,'[[','t'))
  y <- do.call('c',lapply(data,'[[','y'))
  list(t=t,y=y)
}

#' @import ggplot2
plotResults <- function(model,data) {
  k <- model$num.clusters
  z <- model$clusters
  logl <- model$logl
  p <- ggplot() + labs(title=sprintf('LL=%.02f',logl))
  for (i in seq(k)) {
    sp <- sequencePoints(data[z==i])
    d <- cbind(as.data.frame(sp),z=i)
    p <- p + geom_point(aes(t,y),data=d,alpha=0.25)
  }
  traj <- clusterMeans(model,100)
  traj <- data.frame(t=traj$t,y=as.numeric(traj$y),z=rep(1:k,each=100))
  p <- p + geom_line(aes(t,y),data=traj,color='red',size=1.5)
  p <- p + facet_wrap(~z)
  p
}

################################################################################
### Utility functions

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

linearDesign <- function(t,t.range) {
  tbar <- mean(t.range)
  cbind(1,t-tbar)
}

#' @import splines
bSplineDesign <- function(t,t.range,num.bases) {
  order <- 4
  spacing <- seq(t.range[1],t.range[2],length=num.bases-order+2)
  bknots <- spacing[c(1,length(spacing))]
  iknots <- spacing[-c(1,length(spacing))]
  d <- bs(t,df=NULL,iknots,order-1,TRUE,bknots)
  unname(d)
}

sqExpGram <- function(x1,x2,amp,len) {
  d <- outer(x1,x2,'-')
  amp^2 * exp(-0.5 * d^2 / len^2)
}

sqExpCov <- function(x,amp,len,sig) {
  k <- sqExpGram(x,x,amp,len)
  k + diag(sig^2,length(x))
}
