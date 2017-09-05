rm(list=ls())
library(boot)

#==============
# basic values
#==============
nsim <- 3

nc <- 6
nc.max <- 10
ns <- 50
nb <- 20
nms <- 1
nmb <- 1000

nl <- 2000
nl.miss <- round(nl * 0)
niter <- 10000
nburn <- round(niter * .5)

gamma <- .75

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
  return(z)
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

update.theta <- function(param, jump, nl, nc, y, w, nms, nmb, gamma) {
  vmat.old <- param$vmat
  vmat.tmp <- tnorm(nl * (nc-1), lo = 0, hi = 1, mu = vmat.old[,-nc], sig = jump)
  proposed <- cbind(matrix(vmat.tmp, nl, nc-1), 1)

  for (j in 1:(nc-1)) {
    # last column has to be 1
    vmat.new <- vmat.old
    vmat.new[,j] <- proposed[,j]

    theta.old <- vmat2theta(vmat = vmat.old)
    theta.new <- vmat2theta(vmat = vmat.new)

    adj <- fix.MH(lo = 0, hi = 1, old1 = vmat.old[,j], new1 = vmat.new[,j], jump = jump)

    prob_y.old <- get.logl(theta = theta.old, phi = param$phi, y = y, nmat = nms)
    prob_y.new <- get.logl(theta = theta.new, phi = param$phi, y = y, nmat = nms)

    prob_w.old <- get.logl(theta = theta.old, phi = param$omega, y = w, nmat = nmb)
    prob_w.new <- get.logl(theta = theta.new, phi = param$omega, y = w, nmat = nmb)

    pold <- rowSums(prob_y.old, na.rm=T) + 
            rowSums(prob_w.old, na.rm=T) + 
            dbeta(vmat.old[,j], 1, gamma, log = T)
    pnew <- rowSums(prob_y.new, na.rm=T) + 
            rowSums(prob_w.new, na.rm=T) + 
            dbeta(vmat.new[,j], 1, gamma, log = T)

    k <- acceptMH(p0 = pold, p1 = pnew + adj, x0 = vmat.old[,j], x1 = vmat.new[,j], BLOCK = F)
    vmat.old[,j] <- k$x
  }
  vmat <- vmat.old
  theta <- vmat2theta(vmat=vmat)
  list(theta=theta, vmat=vmat)
} # update.theta

update.phi <- function(param, jump, nc, ns, y, nms, a.phi, b.phi) {
  phi.old <- param$phi
  proposed <- matrix(tnorm(nc*ns, lo = 0, hi = 1, mu = phi.old, sig = jump), nc, ns)
    
  for (i in 1:nc) {
    phi.new <- phi.old
    phi.new[i,] <- proposed[i,]

    adj <- fix.MH(lo = 0, hi = 1, old1 = phi.old[i,], new1 = phi.new[i,], jump = jump)

    prob.old <- get.logl(theta = param$theta, phi = phi.old, y = y, nmat = nms)
    prob.new <- get.logl(theta = param$theta, phi = phi.new, y = y, nmat = nms)

    pold <- colSums(prob.old, na.rm=T) + dbeta(phi.old[i,], a.phi, b.phi, log = T)
    pnew <- colSums(prob.new, na.rm=T) + dbeta(phi.new[i,], a.phi, b.phi, log = T)

    k <- acceptMH(p0 = pold, p1 = pnew + adj, x0 = phi.old[i,], x1 = phi.new[i,], BLOCK = F)
    phi.old[i,] <- k$x
  }
  phi <- phi.old
  return(phi)
} # update.phi

update.omega <- function(param, jump, nc, nb, w, nmb, a.phi, b.phi) {
  omega.old <- param$omega
  proposed <- matrix(tnorm(nc*nb, lo = 0, hi = 1, mu = omega.old, sig = jump), nc, nb)
    
  for (i in 1:nc) {
    omega.new <- omega.old
    omega.new[i,] <- proposed[i,]

    adj <- fix.MH(lo = 0, hi = 1, old1 = omega.old[i,], new1 = omega.new[i,], jump = jump)

    prob.old <- get.logl(theta = param$theta, phi = omega.old, y = w, nmat = nmb)
    prob.new <- get.logl(theta = param$theta, phi = omega.new, y = w, nmat = nmb)

    pold <- colSums(prob.old, na.rm=T) + dbeta(omega.old[i,], a.phi, b.phi, log = T)
    pnew <- colSums(prob.new, na.rm=T) + dbeta(omega.new[i,], a.phi, b.phi, log = T)

    k <- acceptMH(p0 = pold, p1 = pnew + adj, x0 = omega.old[i,], x1 = omega.new[i,], BLOCK = F)
    omega.old[i,] <- k$x
  }
  omega<- omega.old
  return(omega)
} # update.phi

for (k in 1:nsim) {
#===============
# simulate data
#===============
vmat <- matrix(rbeta(nl*nc, 1, gamma), nl, nc)
vmat[,nc] <- 1

theta <- vmat2theta(vmat)
colMeans(vmat)
colMeans(theta)
table(rowSums(theta))

phi <- matrix(rgamma(nc*ns, 5, 1), nc, ns)
for (i in 1:nc) {
  phi[i,] <- phi[i,] / sum(phi[i,])
}

theta.x.phi <- theta %*% phi
y <- matrix(rbinom(nl*ns, nms, theta.x.phi), nl, ns)
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

omega <- matrix(rgamma(nc*nb, 5, 1), nc, nb)
for (i in 1:nc) {
  omega[i,] <- omega[i,] / sum(omega[i,])
}

theta.x.omega <- theta %*% omega
w <- matrix(rbinom(nl*nb, nmb, theta.x.omega), nl, nb)

#======
# test
#======
param <- list()
param$vmat <- cbind(vmat[,-nc], matrix(.9, nl, nc.max-nc+1))
param$vmat[,nc.max] <- 1
param$theta <- vmat2theta(param$vmat)
param$phi <- rbind(phi, matrix(1/ns, nc.max-nc, ns))
param$omega <- rbind(omega, matrix(1/nb, nc.max-nc, nb))

theta.post <- vmat.post <- array(, dim=c(nl, nc.max, niter))
theta.post[,,1] <- param$theta
vmat.post[,,1] <- param$vmat

phi.post <- array(, dim=c(nc.max, ns, niter))
phi.post[,,1] <- param$phi

omega.post <- array(, dim=c(nc.max, nb, niter))
omega.post[,,1] <- param$omega

for (i in 2:niter) {
  theta.up <- update.theta(param, jump=.1, nl, nc.max, y, w, nms, nmb, gamma=.1)
  param$theta <- theta.up$theta
  param$vmat <- theta.up$vmat
  theta.post[,,i] <- param$theta
  vmat.post[,,i] <- param$vmat

  param$phi <- update.phi(param, jump=.1, nc.max, ns, y, nms, a.phi=.01, b.phi=.99)
  phi.post[,,i] <- param$phi

  param$omega <- update.omega(param, jump=.1, nc.max, nb, w, nmb, a.phi=.01, b.phi=.99)
  omega.post[,,i] <- param$omega
}

theta.est <- apply(theta.post[,,nburn:niter], 1:2, median)
phi.est <- apply(phi.post[,,nburn:niter], 1:2, median)
omega.est <- apply(omega.post[,,nburn:niter], 1:2, median)

pdf(file=paste(c('c:/A. UFL/remote sensing_', nl, '-', nl.miss, '_', nms, '-', nmb, '_', k, '.pdf'), collapse=''), width=9, height=7)
par(mfrow=c(4,nc))
par(mar=c(1,1,1,1))

for (i in 1:nc) {
  plot(theta.est[,i] ~ theta[,i], xlab='', ylab='', axes=F)
  abline(0, 1, col=2)
  box()
} # i
for (i in 1:nc) {
  plot(phi.est[i,] ~ phi[i,], xlab='', ylab='', axes=F)
  abline(0, 1, col=2)
  box()
} # i
for (i in 1:nc) {
  plot(omega.est[i,] ~ omega[i,], xlab='', ylab='', axes=F)
  abline(0, 1, col=2)
  box()
} # i
boxplot(theta.est, xlab='', ylab='', axes=F, 
        col=c(rep('yellow',nc), rep(4,nc.max-nc)), 
        border=c(rep(2,nc), rep(4,nc.max-nc)))
box()

dev.off()

out <- list(nc=nc, nc.max=nc.max, nl=nl, nl.miss=nl.miss, 
            ns=ns, nb=nb, nms=nms, nmb=nmb, y=y, w=w, 
            vmat=vmat, theta=theta, phi=phi, omega=omega, 
            vmat.post=vmat.post, theta.post=theta.post, 
            phi.post=phi.post, omega.post=omega.post)

save(out, file=paste(c('c:/A. UFL/remote sensing_', nl, '-', nl.miss, '_', nms, '-', nmb, '_', k, '.RData'), collapse=''))

} # k




