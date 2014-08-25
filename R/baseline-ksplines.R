# --------------------------------------------------------------------------
# top-level functions

ksplines <- function(formula, data, ngroups, xrange, nbases,
                     intercept=TRUE, degree=3, tol=1e-4)
{
  mm <- match.call()
  mf <- model.frame(formula, data)

  y <- mf[[1]]
  x <- mf[[2]]
  by <- mf[[3]]

  curveset <- .make_curveset(x, y, by)
  basis <- .bspline_basis(xrange, nbases, intercept, degree)
  model <- .new_ksplines_model(ngroups, nbases, basis)

  model$train_info <- list()
  model$train_info$call <- mm
  model$train_info$mframe <- mf

  model <- .init_ksplines_model(curveset, model)
  emfit <- .ksplines_em(curveset, model, tol=tol)
  emfit
}

ksplines_infer <- function(model, data=NULL)
{
  if (is.null(data))
    mf <- model$train_info$mframe
  else
    mf <- model.frame(formula(model$train_info$mframe), data)

  y <- mf[[1]]
  x <- mf[[2]]
  by <- mf[[3]]

  curveset <- .make_curveset(x, y, by)
  m <- length(curveset$curves)
  curve_id <- vapply(curveset$curves, "[[", vector(class(by), 1), "by")

  group <- integer(m)
  posteriors <- matrix(NA, m, model$ngroups)
  colnames(posteriors) <- 1:model$ngroups

  for (i in 1:m)
  {
    curve <- curveset$curves[[i]]
    inf <- .ksplines_inference(curve, model)

    group[i] <- which.max(inf$z)
    posteriors[i, ] <- inf$z
  }

  posteriors <- cbind(curve_id, group, posteriors)
  colnames(posteriors)[1] <- colnames(mf)[3]

  posteriors
}

plot.ksplines_model <- function(x, ...)
{
  model <- x
  xrange <- eval(model$train_info$call$xrange)

  x <- seq(xrange[1], xrange[2], length=50)
  y <- predict(model, x)

  xname <- colnames(model$train_info$mframe)[2]
  yname <- colnames(model$train_info$mframe)[1]

  matplot(x, y, xlab=xname, ylab=yname, ...)
}

predict.ksplines_model <- function(object, x, ...)
{
  model <- object
  X <- model$basis(x)
  y <- X %*% model$beta

  colnames(y) <- 1:model$ngroups
  y
}

# --------------------------------------------------------------------------
# estimation functions

.ksplines_em <- function(curveset, model, tol=1e-8)
{
  iter <- 0
  convergence <- 1
  likelihood_old <- 0

  while (convergence > tol)
  {
    iter <- iter + 1

    likelihood <- 0
    ss <- .new_ksplines_suffstats(model, length(curveset$curves))

    for (i in seq(along=curveset$curves))
    {
      curve <- curveset$curves[[i]]
      estep <- .ksplines_estep(curve, i, model, ss)
      ss <- estep$ss
      likelihood <- likelihood + estep$likelihood
    }

    model <- .ksplines_mle(model, ss, curveset)

    convergence <- (likelihood_old - likelihood) / likelihood_old
    likelihood_old <- likelihood

    .msg(sprintf("iter=%04d, likelihood=%.2f, convergence=%.8f",
                 iter, likelihood, convergence))
  }

  list(model=model, likelihood=likelihood)
}

.ksplines_estep <- function(curve, curve.ix, model, ss)
{
  inf <- .ksplines_inference(curve, model)

  z <- inf$z
  x <- curve$x
  y <- curve$y
  X <- model$basis(x)
  eta1 <- crossprod(X)
  eta2 <- t(X) %*% y

  for (i in 1:model$ngroups)
  {
    ss$theta_suffstats[i] <- ss$theta_suffstats[i] + z[i]
    ss$beta_eta1[, , i] <- ss$beta_eta1[, , i] + z[i] * eta1
    ss$beta_eta2[, i] <- ss$beta_eta2[, i] + z[i] * eta2
  }
  ss$z[, curve.ix] <- z

  list(ss=ss, likelihood=inf$likelihood)
}

.ksplines_mle <- function(model, ss, curveset, alpha=1, fudge=1e-3)
{
  n <- vapply(curveset$curves, '[[', integer(1), 'n')
  x <- do.call('c', lapply(curveset$curves, '[[', 'x'))
  X <- model$basis(x)
  y <- do.call('c', lapply(curveset$curves, '[[', 'y'))
  z <- ss$z[, rep(seq(along=n), n)]

  theta <- ss$theta_suffstats + alpha
  model$theta <- theta / sum(theta)

  for (i in 1:model$ngroups)
  {
    A <- ss$beta_eta1[, , i] + diag(fudge, model$nbases)
    b <- ss$beta_eta2[, i]
    beta <- solve(A, b)
    y.hat <- X %*% beta
    y.res <- y - y.hat
    sigma.sq <- sum(z[i, ] * y.res^2) / sum(z[i, ])
    model$beta[, i] <- beta
    model$sigma[i] <- sqrt(sigma.sq)
  }

  model
}

# --------------------------------------------------------------------------
# inference functions


.ksplines_inference <- function(curve, model)
{
  X <- model$basis(curve$x)
  logpost <- numeric(model$ngroups)

  for (i in 1:model$ngroups)
  {
    lp <- log(model$theta[i])
    mean <- X %*% model$beta[, i]
    sigma <- diag(model$sigma[i]^2, nrow(X))
    logpost[i] <- lp + mvtnorm::dmvnorm(curve$y, mean, sigma, log=TRUE)
  }

  likelihood <- logsumexp(logpost)
  z <- exp(logpost - likelihood)

  list(z=z, likelihood=likelihood)
}

# --------------------------------------------------------------------------
# model functions

.new_ksplines_model <- function(G, P, basis)
{
  model <- structure(list(), class="ksplines_model")
  model$ngroups <- G
  model$nbases <- P
  model$basis <- basis
  model$theta <- rep(NA, G)
  model$beta <- matrix(NA, P, G)
  model$sigma <- rep(NA, G)
  model
}

.init_ksplines_model <- function(curveset, model)
{
  G <- model$ngroups
  P <- model$nbases

  theta <- runif(G)
  model$theta <- theta / sum(theta)
  model$sigma <- rep(1, G)

  init_groups <- sample(G, length(curveset$curves), TRUE)

  for (i in 1:model$ngroups)
  {
    curves <- curveset$curves[init_groups == i]
    x <- do.call("c", lapply(curves, "[[", "x"))
    y <- do.call("c", lapply(curves, "[[", "y"))
    X <- model$basis(x)
    A <- crossprod(X) + diag(1e-2, P)
    b <- t(X) %*% y
    model$beta[, i] <- solve(A, b)
  }

  model
}

.new_ksplines_suffstats <- function(model, ncurves)
{
  G <- model$ngroups
  P <- model$nbases

  ss <- structure(list(), class="ksplines_suffstats")
  ss$theta_suffstats <- numeric(G)
  ss$beta_eta1 <- array(0, c(P, P, G))
  ss$beta_eta2 <- matrix(0, P, G)
  ss$z <- matrix(NA, G, ncurves)
  ss
}

.bspline_basis <- function(xrange, nbases, intercept, degree=3)
{
  from <- xrange[1]
  to   <- xrange[2]

  if (intercept)
      nknots <- nbases - degree + 1
  else
      nknots <- nbases - degree + 2

  knots <- seq(from, to, length=nknots)
  boundary <- knots[ c(1, nknots)]
  interior <- knots[-c(1, nknots)]

  basis <- function(x)
  {
    X <- splines::bs(x, degree=degree, knots=interior,
                     Boundary.knots=boundary, intercept=intercept)

    unname(X)
  }

  basis
}

# --------------------------------------------------------------------------
# Data Functions

.make_curveset <- function(x, y, by)
{
  byf <- as.factor(by)
  curve_data <- data.frame(x=x, y=y, id=as.integer(byf), by=by)

  curveset <- structure(list(), class="curveset")
  curveset$curves <- lapply(split(curve_data, curve_data$id), .make_curve)
  curveset$xrange <- range(x)
  curveset$yrange <- range(y)
  curveset
}

.make_curve <- function(curve_data)
{
  curve <- structure(list(), class="curve")
  curve$id <- curve_data$id[1]
  curve$by <- curve_data$by[1]
  curve$n <- nrow(curve_data)
  curve$x <- curve_data$x
  curve$y <- curve_data$y
  curve
}

# --------------------------------------------------------------------------
# utility functions

.msg <- function(s)
{
  time <- format(Sys.time(), "%X")
  cat(sprintf("%s %s\n", time, s))
}

# --------------------------------------------------------------------------
# demo

## library(ggplot2)

## source("ksplines.R")

## load("patients.RData")

## # fit the model to the 'patients' data

## fit <- ksplines(y ~ x:ptid, data=patients, ngroups=8,
##                 xrange=c(0, 15), nbases=4)

## # extract the model and likelihood from the fit

## model <- fit$model
## likelihood <- fit$likelihood

## # plot the mean curves

## plot(model)

## # extract the posteriors and merge with the data

## posteriors <- ksplines_infer(model)
## patients <- merge(patients, posteriors[, c("ptid", "group")])

## # plot the patients by the most likely posterior group

## ggplot(patients) +
##   geom_point(aes(x, y), alpha=0.5) +
##   facet_wrap(~ group)
