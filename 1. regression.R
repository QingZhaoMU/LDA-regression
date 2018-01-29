rm(list=ls())
library(boot)

path <- 'c:/A. UFL/1. simulation/'

nsim <- 1

niter <- 100
nburn <- round(niter * .8)
adapt <- round(niter * .5)

#==============
# basic values
#==============
nc <- 3 # number of communities
nx <- 3 # number of covariates plus one
        # the first column for x is always 1, for intercept

nl <- 500 # number of sites
ns <- 200 # number of species
nm <- 1 # number of surveys per site
sigma <- 1 # standard deviation of the errors in logistic regression

for (k in 1:nsim) {
#=================
# data simulation
#=================
# beta's, the first row is the intercept
beta <- matrix(runif(nx*(nc-1), -1, 1), nrow=nx, ncol=nc-1)

# covariates, the first column is 1 for the intercept
x <- cbind(1, matrix(rnorm(nl*(nx-1), 0, 1), nrow=nl, ncol=nx-1))

# theta, proportion of communities for each site
mu <- x %*% beta
vmat <- matrix(exp(rnorm(nl*nc, mu, sigma)), nrow=nl, ncol=nc)
vmat[,nc] <- 1
theta <- vmat / rowSums(vmat)
#colMeans(theta)

# phi, proportion of species for each community
phi <- matrix(rbeta(nc*ns, .2, .8), nrow=nc, ncol=ns)

# pi, probability of occurrence for a species in a site
pi <- theta %*% phi

# number of detections during nmat surveys
y <- matrix(rbinom(nl*ns, nm, pi), nrow=nl, ncol=ns)

#==================
# define functions
#==================
tnorm <- function(n, lo, hi, mu, sig) {
  # generates truncated normal variates based on cumulative normal distribution normal truncated lo and hi

  if (length(lo) == 1 & length(mu) > 1) 
    lo <- rep(lo, length(mu))
  if (length(hi) == 1 & length(mu) > 1) 
    hi <- rep(hi, length(mu))
    
  q1 <- pnorm(lo, mu, sig)  #cumulative distribution
  q2 <- pnorm(hi, mu, sig)  #cumulative distribution
    
  z <- runif(n, q1, q2)
  z <- qnorm(z, mu, sig)
  z[z == -Inf] <- lo[z == -Inf]
  z[z == Inf] <- hi[z == Inf]
  z
} # tnorm

fix.MH <- function(lo, hi, old1, new1, jump) {
  jold <- pnorm(hi, mean = old1, sd = jump) - pnorm(lo, mean = old1, sd = jump)
  jnew <- pnorm(hi, mean = new1, sd = jump) - pnorm(lo, mean = new1, sd = jump)
  log(jold) - log(jnew)  #add this to pnew
} # fix.MH

get.logl <- function(theta, phi, y, nm) {
  prob <- theta %*% phi
  cond <- prob < 1e-05
  prob[cond] <- 1e-05
  cond <- prob > 0.99999
  prob[cond] <- 0.99999
  dbinom(y, size = nm, prob = prob, log = T)
} # get.logl

acceptMH <- function(p0, p1, x0, x1, BLOCK) {
  # accept for M, M-H if BLOCK, then accept as a block, otherwise, accept individually
    
  nz <- length(x0)  #no. to accept
  if (BLOCK) 
    nz <- 1
    
  a <- exp(p1 - p0)  #acceptance PR
  z <- runif(nz, 0, 1)
  keep <- which(z < a)
    
  if (BLOCK & length(keep) > 0) 
    x0 <- x1
  if (!BLOCK) 
    x0[keep] <- x1[keep]
  accept <- length(keep)
    
  list(x = x0, accept = accept)
} # acceptMH

update.theta <- function(param, jump, nl, nc, y, x, nm) {
  phi <- param$phi
  beta <- param$beta
  mu <- x %*% beta

  vmat.ori <- vmat.old <- param$vmat
  vmat.tmp <- exp(rnorm(nl*(nc-1), mean=log(vmat.old[,-nc]), sd=jump[,-nc]))
  proposed <- cbind(matrix(vmat.tmp, nl, nc-1), 1)

  for (j in 1:(nc-1)) {
    # last column has to be 1
    vmat.new <- vmat.old
    vmat.new[,j] <- proposed[,j]

    theta.old <- vmat.old / rowSums(vmat.old)
    theta.new <- vmat.new / rowSums(vmat.new)

    prob.old <- get.logl(theta=theta.old, phi=phi, y=y, nm=nm)
    prob.new <- get.logl(theta=theta.new, phi=phi, y=y, nm=nm)

    pold <- rowSums(prob.old) + dnorm(log(vmat.old[,j]), mu[,j], sigma, log = T)
    pnew <- rowSums(prob.new) + dnorm(log(vmat.new[,j]), mu[,j], sigma, log = T)

    k <- acceptMH(p0=pold, p1=pnew, x0=vmat.old[,j], x1=vmat.new[,j], BLOCK=F)
    vmat.old[,j] <- k$x
  }
  vmat <- vmat.old
  theta <- vmat / rowSums(vmat)
  list(theta=theta, vmat=vmat, accept=vmat.ori!=vmat.old)
} # update.theta

update.phi <- function(param, jump, nc, ns, y, nm, a.phi, b.phi) {
  theta <- param$theta

  phi.ori <- phi.old <- param$phi
  proposed <- matrix(tnorm(nc*ns, lo=0, hi=1, mu=phi.old, sig=jump), nc, ns)
  adj <- fix.MH(lo=0, hi=1, old1=phi.old, new1=proposed, jump=jump)

  for (j in 1:nc) {
    phi.new <- phi.old
    phi.new[j,] <- proposed[j,]

    prob.old <- get.logl(theta=theta, phi=phi.old, y=y, nm=nm)
    prob.new <- get.logl(theta=theta, phi=phi.new, y=y, nm=nm)

    pold <- colSums(prob.old) + dbeta(phi.old[j,], a.phi, b.phi, log=T)
    pnew <- colSums(prob.new) + dbeta(phi.new[j,], a.phi, b.phi, log=T)

    k <- acceptMH(p0=pold, p1=pnew + adj[j,], x0=phi.old[j,], x1=phi.new[j,], BLOCK=F)
    phi.old[j,] <- k$x
  }
  phi <- phi.old
  list(phi=phi, accept=phi.ori!=phi.old)
} # update.phi

update.beta <- function(param, jump, nx, nc, x) {
  vmat <- param$vmat

  beta.old <- param$beta
  beta.new <- matrix(rnorm(length(beta.old), beta.old, jump), nx, nc-1)
  mu.old <- x %*% beta.old
  mu.new <- x %*% beta.new

  prob.old <- dnorm(log(vmat[,-nc]), mu.old, sigma, log=T)
  prob.new <- dnorm(log(vmat[,-nc]), mu.new, sigma, log=T)

  pold <- sum(prob.old) + sum(dnorm(beta.old, 0, 5, log=T))
  pnew <- sum(prob.new) + sum(dnorm(beta.new, 0, 5, log=T))

  k <- acceptMH(p0=pold, p1=pnew, x0=beta.old, x1=beta.new, BLOCK=F)
  beta <- k$x
  list(beta=beta, accept=beta!=beta.old)
} # update.beta

jumpTune <- function(accept, jump, ni, adapt=2000, low=.3, high=.8) {
  nstart <- ifelse(ni>=100, ni-99, 1)
  if (ni > adapt) {
    jump.new <- jump
  } else if (ni %% 10 > 0) {
    jump.new <- jump
  } else {
    accept.rate <- apply(accept[,,nstart:ni], 1:2, mean)
    jump.new <- jump
    jump.new[which(accept.rate < low)] <- 
      jump[which(accept.rate < low)] * rnorm(length(which(accept.rate < low)),.5,.01)
    jump.new[which(accept.rate > high)] <- 
      jump[which(accept.rate > high)] * rnorm(length(which(accept.rate > high)),2,.01)
  }
  jump <- jump.new
  list(jump=jump)
} # jumpTune

#===========
# run model
#===========
theta.jump <- matrix(.05, nrow=nl, ncol=nc)
phi.jump <- matrix(.05, nrow=nc, ncol=ns)
beta.jump <- matrix(.05, nrow=nx, ncol=nc-1)

param <- list()

param$vmat <- vmat
param$theta <- theta
param$phi <- phi
param$beta <- beta

vmat.post <- theta.post <- theta.jump.post <- theta.accept <- array(, dim=c(nl, nc, niter))
vmat.post[,,1] <- param$vmat
theta.post[,,1] <- param$theta
theta.jump.post[,,1] <- theta.jump
theta.accept[,,1] <- FALSE

phi.post <- phi.jump.post <- phi.accept <- array(, dim=c(nc, ns, niter))
phi.post[,,1] <- param$phi
phi.jump.post[,,1] <- phi.jump
phi.accept[,,1] <- FALSE

beta.post <- beta.jump.post <- beta.accept <- array(, dim=c(nx, nc-1, niter))
beta.post[,,1] <- param$beta
beta.jump.post[,,1] <- beta.jump
beta.accept[,,1] <- FALSE

for (i in 2:niter) {
  theta.up <- update.theta(param, jump=theta.jump, nl, nc, y, x, nm)
  vmat.post[,,i] <- param$vmat <- theta.up$vmat
  theta.post[,,i] <- param$theta <- theta.up$theta
  theta.accept[,,i] <- theta.up$accept

  phi.up <- update.phi(param, jump=phi.jump, nc, ns, y, nm, a.phi=1, b.phi=1)
  phi.post[,,i] <- param$phi <- phi.up$phi
  phi.accept[,,i] <- phi.up$accept

  beta.up <- update.beta(param, jump=beta.jump, nx, nc, x)
  beta.post[,,i] <- param$beta <- beta.up$beta
  beta.accept[,,i] <- beta.up$accept

  theta.jump.up <- jumpTune(accept=theta.accept, jump=theta.jump, ni=i, adapt=adapt, low=.3, high=.8)
  theta.jump.post[,,i] <- theta.jump <- theta.jump.up$jump

  phi.jump.up <- jumpTune(accept=phi.accept, jump=phi.jump, ni=i, adapt=adapt, low=.3, high=.8)
  phi.jump.post[,,i] <- phi.jump <- phi.jump.up$jump

  beta.jump.up <- jumpTune(accept=beta.accept, jump=beta.jump, ni=i, adapt=adapt, low=.3, high=.8)
  beta.jump.post[,,i] <- beta.jump <- beta.jump.up$jump
}

#==============
# save results
#==============
theta.est <- apply(theta.post[,,(nburn+1):niter], 1:2, median)
phi.est <- apply(phi.post[,,(nburn+1):niter], 1:2, median)
beta.est <- apply(beta.post[,,(nburn+1):niter], 1:2, median)

#===============
# graph results
#===============
par(mfrow=c(4, nc))

for (i in 1:nc) {
  plot(theta.est[,i] ~ theta[,i], xlim=c(0,1), ylim=c(0,1))
  abline(0, 1, col=2)
}

for (i in 1:nc) {
  plot(phi.est[i,] ~ phi[i,], xlim=c(0,1), ylim=c(0,1))
  abline(0, 1, col=2)
}

for (i in 1:(nc-1)) {
  for (j in 1:nx) {
    plot(beta.post[j,i,], type='l')
    abline(h=beta[j,i], col=2)
  }
}

#==============
# save results
#==============
out <- list()
out$x <- x
out$y <- y
out$nc <- nc
out$nm <- nm
out$sigma <- sigma
out$vmat <- vmat
out$phi <- phi
out$beta <- beta
out$theta.est <- theta.est
out$phi.est <- phi.est
out$beta.post <- beta.post[,,(nburn+1):niter]
out$beta.est <- beta.est

save(out, file=paste(c(path, 'regression_', nc, '_', nx, '_', nm, '_', k, '.RData'), collapse=''))

} # k


