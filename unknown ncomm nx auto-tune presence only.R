rm(list=ls())
library(boot)

#==============
# basic values
#==============
nsim <- 3

nl <- 2000
nl.miss <- 2000

niter <- 10000
nburn <- round(niter * .5)

nc <- 5
nc.max <- 10
nx.true <- 3
nx <- 2

ns <- 25
nmat <- 1
jump.beta <- .03

gamma <- .6
beta <- matrix(c(-2.15, .15, 1.25, -1.5, 0, 
                 .7, -2.35, 1.85, -.5, 0, 
                 .85, -0.8, -.85, 1.55, 0), 3, 5, byrow=T)
beta <- beta[1:nx.true, 1:nc]

#==================
# define functions
#==================
vmat2theta <- function(vmat) {
  nl <- dim(vmat)[1]
  nc <- dim(vmat)[2]
  theta <- vmat
  for (i in 1:nc) {
    if (i==1) {
      prod <- rep(1,nl)
    } else {
      prod <- prod * (1-vmat[,i-1])
    }
    theta[,i] <- vmat[,i] * prod
  } # i
  return(theta)
} # vmat2theta

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

get.logl <- function(theta, phi, y, nmat) {
  prob <- theta %*% phi
  cond <- prob < 1e-05
  prob[cond] <- 1e-05
  cond <- prob > 0.99999
  prob[cond] <- 0.99999
  dbinom(y, size = nmat, prob = prob, log = T)
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

jumpTune <- function(param, adapt=2000, low=.3, high=.8) {
  accept <- param$accept
  jump.old <- param$jump
  ni <- max(which(!(is.na(accept))))

  if (ni > adapt) {
    jump <- jump.old
  } else if (ni %% 100 > 0) {
    jump <- jump.old
  } else {
    accept.rate <- mean(accept[(ni-99):ni])
    if (accept.rate < low) {
      jump <- jump.old / 2
    } else if (accept.rate > high) {
      jump <- jump.old * 2
    } else {
      jump <- jump.old
    }
  }

  return(jump)
} # jumpTune

update.theta <- function(param, jump, nl, nc, y, x, nmat, gamma) {
  eta <- matrix(, nl, nc)
  eta <- inv.logit(logit(param$gamma) + x %*% param$beta)

  v.old <- param$vmat
  v.tmp <- tnorm(nl * (nc-1), lo = 0, hi = 1, mu = v.old[,-nc], sig = jump)
  proposed <- cbind(matrix(v.tmp, nl, nc-1), 1)

  for (j in 1:(nc-1)) {
    # last column has to be 1
    v.new <- v.old
    v.new[,j] <- proposed[,j]

    theta.old <- vmat2theta(vmat = v.old)
    theta.new <- vmat2theta(vmat = v.new)

    adj <- fix.MH(lo = 0, hi = 1, old1 = v.old[,j], new1 = v.new[,j], jump = jump)

    prob.old <- get.logl(theta = theta.old, phi = param$phi, y = y, nmat = nmat)
    prob.new <- get.logl(theta = theta.new, phi = param$phi, y = y, nmat = nmat)

    pold <- rowSums(prob.old, na.rm=T) + dbeta(v.old[,j], eta[,j], 1-eta[,j], log = T)
    pnew <- rowSums(prob.new, na.rm=T) + dbeta(v.new[,j], eta[,j], 1-eta[,j], log = T)

    k <- acceptMH(p0 = pold, p1 = pnew + adj, x0=v.old[,j], x1=v.new[,j], BLOCK = F)
    v.old[,j] <- k$x
  }
  vmat <- v.old
  theta <- vmat2theta(vmat=vmat)
  list(theta=theta, vmat=vmat)
} # update.theta

update.phi <- function(param, jump, nc, ns, y, nmat, a.phi, b.phi) {
  phi.old <- param$phi
  proposed <- matrix(tnorm(nc*ns, lo = 0, hi = 1, mu = phi.old, sig = jump), nc, ns)
    
  for (i in 1:nc) {
    phi.new <- phi.old
    phi.new[i,] <- proposed[i,]

    adj <- fix.MH(lo = 0, hi = 1, old1 = phi.old[i,], new1 = phi.new[i,], jump = jump)

    prob.old <- get.logl(theta = param$theta, phi = phi.old, y = y, nmat = nmat)
    prob.new <- get.logl(theta = param$theta, phi = phi.new, y = y, nmat = nmat)

    pold <- colSums(prob.old, na.rm=T) + dbeta(phi.old[i,], a.phi, b.phi, log = T)
    pnew <- colSums(prob.new, na.rm=T) + dbeta(phi.new[i,], a.phi, b.phi, log = T)

    k <- acceptMH(p0 = pold, p1 = pnew + adj, x0 = phi.old[i,], x1 = phi.new[i,], BLOCK = F)
    phi.old[i,] <- k$x
  }
  phi <- phi.old
  return(phi)
} # update.phi

update.beta <- function(param, nl, nc, x) {
  jump <- param$jump
  vmat <- param$vmat

  gamma.old <- param$gamma
  gamma.new <- inv.logit(rnorm(1, logit(gamma.old), jump))

  beta.old <- param$beta
  beta.new <- matrix(rnorm(length(beta.old), beta.old, jump), nx, nc)

  eta.old <- eta.new <- matrix(, nl, nc)
  eta.old <- inv.logit(logit(gamma.old) + x %*% beta.old)
  eta.new <- inv.logit(logit(gamma.new) + x %*% beta.new)

  prob.old <- prob.new <- matrix(, nl, nc-1)
  for (i in 1:nl) {
    for (j in 1:(nc-1)) {
      prob.old[i,j] <- dbeta(vmat[i,j], eta.old[i,j], 1-eta.old[i,j], log=T)
      prob.new[i,j] <- dbeta(vmat[i,j], eta.new[i,j], 1-eta.old[i,j], log=T)
    }  #j
  } # i

  pold <- sum(prob.old) + sum(dnorm(beta.old[-(nc+1)], 0, 5, log=T)) + dunif(gamma.old, 0, 1, log=T)
  pnew <- sum(prob.new) + sum(dnorm(beta.new[-(nc+1)], 0, 5, log=T)) + dunif(gamma.new, 0, 1, log=T)

  if (exp(pnew - pold) > runif(1,0,1)) {
    gamma.old <- gamma.new
    beta.old <- beta.new
    accept <- 1
  } else {
    gamma.old <- gamma.old
    beta.old <- beta.old
    accept <- 0
  }
  gamma <- gamma.old
  beta <- beta.old
  list(gamma=gamma, beta=beta, accept=accept)
} # update.beta

for (k in 1:nsim) {
#===============
# simulate data
#===============
# create independent x-values 
x <- matrix(, nl, nx.true)
for (i in 1:nx.true) {
  x[,i] <- sample(seq(-2, 2, length.out=nl), nl, replace=F)
}

eta <- matrix(, nl, nc)
eta <- inv.logit(logit(gamma) + x %*% beta)

vmat <- matrix(rbeta(nl*nc, eta, 1-eta), nl, nc)
vmat[which(vmat==0)] <- 1e-6
vmat[which(vmat==1)] <- 1 - 1e-6
vmat[,nc] <- 1
theta <- vmat2theta(vmat)

vmat.expect <- eta
vmat.expect[,nc] <- 1
theta.expect <- vmat2theta(vmat.expect)

par(mfrow=c(1,nx))
par(mar=c(0,0,0,0))

for (j in 1:nx) {
  plot(eta[,1] ~ x[,j], type='n', ylim=c(0,1), axes=F)
  for (i in 1:nc) {
    points(theta.expect[,i] ~ x[,j], col=i, pch=3)
    box()
  } # i
} # l

colMeans(vmat)
colMeans(theta)
table(rowSums(theta))
table(rowSums(theta.expect))

phi <- matrix(rbeta(nc*ns, .01, .99), nc, ns)
for (i in 1:nc) {
  phi[i,] <- phi[i,] / sum(phi[i,])
}

pi <- theta %*% phi

y <- matrix(rbinom(nl*ns, nmat, pi), nl, ns)
table(y)
table(rowSums(y))
miss <- sample(1:nl, nl.miss, replace=F)
miss <- sort(miss)
y.temp <- y[miss,]
y.temp[which(y.temp==0)] <- NA
y[miss,] <- y.temp
table(y)
length(which(is.na(y)))
table(rowSums(y, na.rm=T))

x <- x[,1:nx]

#======
# test
#======
param <- list()

param$phi <- rbind(phi, matrix(1/ns, nc.max-nc, ns))
param$vmat <- cbind(vmat[,-nc], matrix(.9, nl, nc.max-nc+1))
param$vmat[,nc.max] <- 1
param$theta <- vmat2theta(param$vmat)
param$gamma <- gamma
param$beta <- cbind(beta, matrix(0, nx.true, nc.max-nc))[1:nx,]
param$jump <- jump.beta
param$accept <- rep(NA, niter)
param$accept[1] <- 0

theta.post <- vmat.post <- array(, dim=c(nl, nc.max, niter))
theta.post[,,1] <- param$theta
vmat.post[,,1] <- param$vmat

phi.post <- array(, dim=c(nc.max, ns, niter))
phi.post[,,1] <- param$phi

gamma.post <- numeric(niter)
gamma.post[1] <- param$gamma

beta.post <- array(, dim=c(nx, nc.max, niter))
beta.post[,,1] <- param$beta

jump.post <- numeric(niter)
jump.post[1] <- jump.beta

for (i in 2:niter) {
  theta.up <- update.theta(param, jump=.1, nl, nc.max, y, x, nmat, gamma=gamma)
  param$theta <- theta.up$theta
  param$vmat <- theta.up$vmat
  theta.post[,,i] <- param$theta
  vmat.post[,,i] <- param$vmat

  param$phi <- update.phi(param, jump=.1, nc.max, ns, y, nmat, a.phi=.01, b.phi=.99)
  phi.post[,,i] <- param$phi

  beta.up <- update.beta(param, nl, nc.max, x)
  param$gamma <- beta.up$gamma
  param$beta <- beta.up$beta
  param$accept[i] <- beta.up$accept
  gamma.post[i] <- param$gamma
  beta.post[,,i] <- param$beta

  param$jump <- jumpTune(param, adapt=2000, low=.3, high=.8)
  jump.post[i] <- param$jump
}

theta.est <- apply(theta.post[,,nburn:niter], 1:2, median)
phi.est <- apply(phi.post[,,nburn:niter], 1:2, median)

#===============
# graph results
#===============
pdf(file=paste(c('c:/A. UFL/beta regression_nx=', nx, '-', nx.true, '_nm=', nmat, '_nl[miss]=', nl, '[', nl.miss, ']_', k, '.pdf'), collapse=''), 
    width=7, height=7)
par(mfrow=c(nx.true,nc))
par(mar=c(5,4.5,3,.5))

boxplot(theta.est, border=c(rep('red',nc),rep('blue',nc.max-nc)))

plot(gamma.post, type='l', 
     xlab='Iteration', ylab='Esimated Gamma')
abline(h=gamma, col=2) 

for (i in 1:nx) {
  for (j in 1:(nc-1)) {
    plot(beta.post[i,j,], type='l', 
       xlab='Iteration', ylab='Esimated Beta', 
       main=paste(c('X',i,' Comm',j), collapse=''))
    abline(h=beta[i,j], col=2)
  } # i
} # j

plot(jump.post, type='l'); abline(v=2000, col=2)

dev.off()

out <- list(nc=nc, nc.max=nc.max, nl=nl, nl.miss=nl.miss, 
            ns=ns, nmat=nmat, x=x, y=y, 
            vmat=vmat, theta=theta, phi=phi, 
            gamma=gamma, beta=beta, 
            vmat.post=vmat.post, theta.post=theta.post, phi.post=phi.post, 
            gamma.post=gamma.post, beta.post=beta.post, jump.post=jump.post)

save(out, file=paste(c('c:/A. UFL/beta regression_nx=', nx, '-', nx.true, '_nm=', nmat, '_nl[miss]=', nl, '[', nl.miss, ']_', k, '.RData'), collapse=''))

} # k


