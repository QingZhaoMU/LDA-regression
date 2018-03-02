library(boot)

#==================
# define functions
#==================
#generates multivariate normal with mean 0 and Sigma=sqrtVar1%*%t(sqrtVar1)
rmvnorm1 <- function (sqrtVar1) {
  sqrtVar1 %*% rnorm(ncol(sqrtVar1))
} # rmvnorml

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

  z.ori <- z.old <- param$z
  vmat.old=cbind(exp(z.old),1)
  
  z.tmp = rnorm(nl*(nc-1), mean=z.old, sd=jump)
  z.proposed = matrix(z.tmp,nrow=nl,ncol=nc-1)
  vmat.proposed <- cbind(exp(z.proposed),1)

  for (j in 1:(nc-1)) {
    # last column has to be 1
    vmat.new <- vmat.old
    vmat.new[,j] <- vmat.proposed[,j]

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
  z = log(vmat)[,-nc]
  list(theta=theta, z=z, accept=z.ori!=log(vmat.old[,-nc]))
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

update.beta <- function(param, nx, nc, tx, sqrtVar1, var1) {
  z <- param$z
  
  beta=matrix(, nx, nc-1)
  for (i in 1:(nc-1)){
    pmean <- tx%*%z[,i]
    beta[,i] <- rmvnorm1(sqrtVar1) + var1%*%pmean
  }
  beta 
} # update.beta

jumpTune <- function(accept, jump, ni, low=.3, high=.8) {
  nstart <- ifelse(ni>=100, ni-99, 1)
  accept.rate <- apply(accept[,,nstart:ni], 1:2, mean)
  jump[accept.rate < low]  <- jump[accept.rate < low] * 0.5 
  jump[accept.rate > high] <- jump[accept.rate > high] * 2
  jump
} # jumpTune
