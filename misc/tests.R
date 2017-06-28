p1 <- 10; nsample <- 10; n <- p1 * nsample
d <- data.frame(x1=0, x2=0, y=0, u1=rep(1:p1, each=nsample), u2=rep(1:p1, times=nsample))
d$u1 <- as.factor(d$u1); d$u2 <- as.factor(d$u2)
b1 <- 1; b2 <- -1; sd1 <- 1.5
set.seed(123)
d$x1 <- rnorm(n=n); d$x2 <- rnorm(n=n)
d$y <- b1 * d$x1 + b2 * d$x2 + rep(rnorm(n=p1, sd=sd1), each=nsample) +
  rep(rnorm(n=p1, sd=sd1), times=nsample) + rnorm(n=n)
mod = lme4::lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = F)
mod.r = lm(y ~ 1, data = d)

microbenchmark::microbenchmark(R2.ce.lmerMod(mod, mod.r),
                               R2.ce.lmerMod.cpp(mod, mod.r),
                               R2.ce.lmerMod.cpp2(mod, mod.r),
                               times = 10)

library(profvis)
profvis( {Y <- model.frame(mod)[, 1]
          X <- model.matrix(mod)
          n <- dim(X)[1]
          p <- dim(X)[2]

          X.r <- model.matrix(mod.r)
          p.r <- dim(X.r)[2]

          D <- attr(mod, "pp")$LamtUt
          V <- crossprod(D) + diag(n)

          if (class(mod.r) == "lmerMod") {
            D.r <- attr(mod.r, "pp")$LamtUt
            V.r <- crossprod(D.r) + diag(n)
          }
          if (class(mod.r) == "lm") {
            V.r <- diag(n)
          }

          iV <- solve(V)

          bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
          fitted.values <- X %*% bhat
          R <- Y - X %*% bhat

          Rhat <- matrix(0, nrow = n, ncol = 1)

          for (j in 1:n) {
            r <- R[-j]
            VV <- V[-j, -j]
            iVV <- solve(VV)
            # version using global mean: This works best
            bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]
            v <- V[j, -j]
            Rhat[j,] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
          }

          Yhat <- as.numeric(fitted.values + Rhat)
          SSE.ce <- var(Y - Yhat)

          # reduced model
          iV.r <- solve(V.r)

          bhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV.r %*% Y)
          fitted.values <- X.r %*% bhat
          R <- Y - X.r %*% bhat

          Rhat <- matrix(0, nrow = n, ncol = 1)

          for (j in 1:n) {
            r <- R[-j]
            VV.r <- V.r[-j, -j]
            iVV.r <- solve(VV.r)
            # version using global mean: This works best
            bbhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV %*% R)[1]
            v <- V.r[j, -j]
            Rhat[j] <- as.numeric(bbhat + v %*% iVV.r %*% (r - bbhat))
          }

          Yhat.r <- as.numeric(fitted.values + Rhat)
          SSE.ce.r <- var(Y - Yhat.r)

          R2.ce <- 1 - SSE.ce/SSE.ce.r})

profvis({
  Y <- model.frame(mod)[, 1]
  X <- model.matrix(mod)
  n <- dim(X)[1]
  p <- dim(X)[2]

  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]

  D <- attr(mod, "pp")$LamtUt
  V <- crossprod(D) + diag(n)

  if (class(mod.r) == "lmerMod") {
    D.r <- attr(mod.r, "pp")$LamtUt
    V.r <- crossprod(D.r) + diag(n)
  }
  if (class(mod.r) == "lm") {
    V.r <- diag(n)
  }

  iV <- solve(V)

  bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
  fitted.values <- X %*% bhat
  R <- Y - X %*% bhat

  Rhat = loop_cpp_lmerMod(as.matrix(R), as.matrix(V), as.matrix(iV), X)

  Yhat <- as.numeric(fitted.values + Rhat)
  SSE.ce <- var(Y - Yhat)

  # reduced model
  iV.r <- solve(V.r)

  bhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV.r %*% Y)
  fitted.values <- X.r %*% bhat
  R <- Y - X.r %*% bhat

  Rhat = loop_cpp_lmerMod(as.matrix(R), as.matrix(V.r), as.matrix(iV.r), X.r)

  Yhat.r <- as.numeric(fitted.values + Rhat)
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r
})

## pgls ----
n <- p1 * nsample
d <- data.frame(x=array(0, dim=n), y=0)

b1 <- 1.5
signal <- 0.7

phy <- ape::compute.brlen(ape::rtree(n=n), method = "Grafen", power = 1)
phy.x <- phy

x <- ape::rTraitCont(phy.x, model = "BM", sigma = 1)

e <- signal^0.5 * ape::rTraitCont(phy, model = "BM", sigma = 1) + (1-signal)^0.5 * rnorm(n=n)
d$x <- x[match(names(e), names(x))]

d$y <- b1 * x + e
rownames(d) <- phy$tip.label

z.x <- phylolm::phylolm(y ~ 1, phy=phy, data=d, model="lambda")
lam.x <- round(z.x$optpar, digits=4)
mod <- phylolm::phylolm(y ~ x, phy=phy, data=d, model="lambda", starting.value=.98*lam.x+.01)
mod.r <- lm(y ~ x, data=d)

microbenchmark::microbenchmark(R2.ce.phylolm(mod, mod.r, phy = phy),
                               R2.ce.phylolm.cpp(mod, mod.r, phy = phy),
                               times = 10)

## for rcpp profile ---

# library(lme4)
mod <- lme4::lmer(Reaction ~ Days + (1 | Subject), data=sleepstudy, REML = F)
mod.r <- lme4::lmer(Reaction ~ 1 + (1 | Subject), data=sleepstudy, REML = F)


Y <- model.frame(mod)[, 1]
X <- model.matrix(mod)
n <- dim(X)[1]; p <- dim(X)[2]
D <- attr(mod, "pp")$LamtUt
V <- crossprod(D) + diag(n)
iV <- solve(V)
bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
fitted.values <- X %*% bhat
R <- Y - X %*% bhat

loop_R = function(R, V, iV, X){
  bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]
  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)
    # version using global mean: This works best
    v <- V[j, -j]
    Rhat[j,] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
  }
  Rhat
}

# because I do not know how to deal with sparseMatrix in RcppArmadillo
R = as.matrix(R)
V = as.matrix(V)
iV = as.matrix(iV)

system.time(a <- loop_R(R, V, iV, X)) # 0.616
system.time(b <- loop_cpp_lmerMod(R, V, iV, X)) # 3.857
head(a); head(b) # same results

microbenchmark::microbenchmark(loop_R(R, V, iV, X),
                               loop_cpp(R, V, iV, X),
                               # loop_cpp2(R, V, iV, X),
                               times = 10)


## glmm ----
p1 <- 10; nsample <- 10; n <- p1 * nsample

d <- data.frame(x=0, y=0, u=rep(1:p1, each=nsample))
d$u <- as.factor(d$u)

b1 <- 1
sd1 <- 1.5

d$x <- rnorm(n=n)
prob <- inv.logit(b1 * d$x + rep(rnorm(n=p1, sd=sd1), each=nsample))
d$y <- rbinom(n=n, size=1, prob=prob)

mod <- lme4::glmer(y ~ x + (1 | u), data=d, family="binomial")
mod.r <- lme4::glmer(y ~ 1 + (1 | u), data=d, family="binomial")
z.v <- glm(y ~ x, data=d, family="binomial")


system.time(a <- R2.ce.glmerMod(mod, mod.r)) # 0.616
system.time(b <- R2.ce.glmerMod.cpp(mod, mod.r)) # 3.857
a; b

microbenchmark::microbenchmark(R2.ce.glmerMod(mod, mod.r),
                               R2.ce.glmerMod.cpp(mod, mod.r), times = 100)
