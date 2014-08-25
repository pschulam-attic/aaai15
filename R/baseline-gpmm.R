newGpmm <- function(data,k,a,l) {
  model <- list()
  model$data <- data
  model$k <- k
  model$a <- a
  model$l <- l
  model$s <- 1
  model$probs <- rep(NA,k)
  model$bases <- vector('list',k)
  model$coefs <- vector('list',k)
  class(model) <- c('gpmm', 'list')
  model
}

fitGpmm <- function(data,k,a,l,z,tol=1e-2,maxiter=1e2) {
  initProbs <- function(i) {
    if (is.na(i)) {
      q <- runif(k)
      q <- q / sum(q)
    } else {
      q <- as.integer(1:k == i)
    }
    q
  }
  model <- newGpmm(data,k,a,l)
  q <- apply(as.matrix(z),1,initProbs)
  model <- gpmmMaxStep(data,q,model)
  for (iter in seq(maxiter)) {
    inf <- gpmmExpStep(data,model)
    if (iter == 1) {
      convergence <- Inf
    } else {
      convergence <- (inf$logl - old.logl) / abs(old.logl)
    }
    message(sprintf("LL=%.02f, conv=%.02f",inf$logl,convergence))
    if (convergence < tol) break
    model <- gpmmMaxStep(data,inf$q,model)
    old.logl <- inf$logl
  }
  model
}

gpmmExpStep <- function(data,model) {
  q <- matrix(log(model$probs),model$k,length(data))
  for (i in seq(model$k)) {
    beta <- model$coefs[[i]]
    t.bases <- model$bases[[i]]
    for (j in seq(along=data)) {
      s <- data[[j]]
      X <- sqExpGram(s$t,t.bases,model$a,model$l)
      y.hat <- X %*% beta
      q[i, j] <- sum(dnorm(s$y - y.hat,sd=model$s,log=TRUE))
    }
  }
  logl <- sum(apply(q,2,logsumexp))
  q <- apply(q,2,function(x) exp(x - logsumexp(x)))
  list(q=q,logl=logl)
}

gpmmMaxStep <- function(data,q,model) {
  model$probs <- rowSums(q)/sum(q)
  penalties <- numeric(model$k)
  residuals <- numeric(model$k)
  n.points <- 0
  z <- apply(q,2,which.max)
  for (i in seq(model$k)) {
    z.data <- data[z == i]
    z.t <- do.call('c',lapply(z.data,'[[','t'))
    z.y <- do.call('c',lapply(z.data,'[[','y'))
    n.points <- n.points + length(z.t)
    model$bases[[i]] <- z.t
    K <- sqExpGram(z.t,z.t,model$a,model$l)
    C <- K + diag(model$s^2,nrow(K))
    beta <- solve(C,z.y)
    penalties[i] <- -0.5 * z.y %*% beta - 0.5 * log(det(C))
    residuals[i] <- sum((z.y - K %*% beta)^2)
    model$coefs[[i]] <- beta
  }
  model$penalty <- sum(model$probs * penalties)
  model$s <- sqrt(sum(residuals)/n.points)
  model
}

plotGpmm <- function(model,t.range) {
  t.grid <- seq(t.range[1],t.range[2],0.1)
  p <- ggplot()
  for (i in seq(model$k)) {
    beta <- model$coefs[[i]]
    t.bases <- model$bases[[i]]
    X <- sqExpGram(t.grid,t.bases,model$a,model$l)
    y <- X %*% beta
    d <- data.frame(x=t.grid,y=y,z=i)
    p <- p + geom_line(aes(x,y),data=d)
  }
  p <- p + facet_wrap(~ z)
  p
}

demoGpmm <- function() {
  n <- 50
  beta <- matrix(rnorm(3*2),3,2)
  beta[1,] <- beta[1,] + rnorm(2,0,10)
  x <- runif(n,0,10)
  X <- cbind(1,x,x^2)
  y <- numeric(n)
  z <- sample(2,50,TRUE)
  y[z==1] <- X[z==1,] %*% beta[,1] + rnorm(sum(z==1),0,0.5)
  y[z==2] <- X[z==2,] %*% beta[,2] + rnorm(sum(z==2),0,0.5)
  model <- fitGpmm(x,y,2,1,5,0.5,tol=1e-5)

  plot(x,y)
  xgrid <- seq(0,10,0.1)
  ygrid <- sqExpGram(xgrid,model$x,model$a,model$l) %*% model$coefs
  matlines(xgrid,ygrid)
}
