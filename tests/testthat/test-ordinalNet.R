MASSinst <- require(MASS)
glmnetInst <- require(glmnet)
glmnet201Inst <- if (glmnetInst) packageVersion("glmnet") >= "2.0-1" else FALSE
penalizedInst <- require(penalized)
ordinalgmifsInst <- require(ordinalgmifs)
nnetInst <- require(nnet)

set.seed(4)
n <- 50
k <- 3
p <- 5
betaNonzeroVals <- c(1, 1)
zeta <- seq(-1, 1, length.out=k-1)
beta <- c(betaNonzeroVals, rep(0, p-2))
x <- matrix(rnorm(n*p, sd=2), nrow=n)
eta <- as.vector(x%*%beta) + matrix(rep(zeta, each=n), ncol=(k-1))
cdf <- exp(eta) / (1+exp(eta))
prob <- cbind(cdf, 1) - cbind(0, cdf)
y <- as.factor(sapply(1:n, function(i) sample(1:k, 1, prob=prob[i,])))

yy <- y
levels(yy) <- c(1, 1, 2)
xSD <- apply(x, 2, sd) * sqrt((n-1)/n)
xScaled <- x / rep(xSD, each=n)

###############################################################################
## Make multinomial matrices and link function
###############################################################################
makeMultinomialXLS <- function(x, nLev)
{
    splitX <- split(x, row(x))
    xLS <- lapply(splitX, function(x) {
        xx <- c(1, x)
        m1 <- lapply(1:(k-1), function(i) do.call(rbind, c(rep(list(0), i-1), list(xx), rep(list(0), k-1-i))))
        m2 <- list(do.call(rbind, rep(list(-xx), k-1)))
        mat <- do.call(cbind, c(m1, m2))
        mat
    })
    xLS
}
xLS <- makeMultinomialXLS(xScaled, k)

yMat <- matrix(0, nrow=n, ncol=k)
yMat[cbind(1:n, y)] <- 1

makeMultinomialLink <- function()
{
    lf <- list()
    lf$g <- function(p) sapply(p, function(pp) log(pp / (1-sum(p))))
    lf$h <- function(eta) sapply(eta, function(ee) exp(ee) / (1 + sum(exp(eta))))
    lf$getQ <- function(eta) {
        p <- lf$h(eta)
        q <- -tcrossprod(p)
        diag(q) <- p * (1-p)
        q
    }
    lf
}

multinomialLink <- makeMultinomialLink()

###############################################################################

test_that("Unpenalized ordinal logit matches MASS::polr", {
    if (!MASSinst) skip("MASS not available")
    o2 <- ordinalNet(x, y, lambdaVals=0)
    p2 <- MASS::polr(y~x)
    expect_equal(coef(o2), c(p2$zeta, -p2$coef), check.attributes=F, tolerance=1e-3)
})

test_that("Elastic net binary logistic regression matches glmnet", {
    if (!glmnetInst) skip("glmnet not installed")
    o3 <- ordinalNet(x, yy, alpha=.5, lambdaVals=.1, epsIn=1e-8, epsOut=1e-8)
    g3 <- glmnet(x, yy, family="binomial", alpha=.5, lambda=.1)
    expect_equal(coef(o3), c(-g3$a0, -as.vector(g3$beta)), check.attributes=F, tolerance=1e-3)
})

test_that("Elastic net binary logistic regression matches penalized", {
    if (!penalizedInst) skip("penalized not installed")
    o3 <- ordinalNet(x, yy, alpha=.5, lambdaVals=.1, epsIn=1e-8, epsOut=1e-8)
    p3 <- penalized(yy, x, model="logistic", lambda1=.5*.1*n, lambda2=.5*.1*n, standardize=T, trace=F)
    expect_equal(coef(o3), c(-p3@unpenalized, -p3@penalized), check.attributes=F, tolerance=1e-3)
})

test_that("Elastic net binary logistic regression with positve constraints and some unpenalized terms matches penalized", {
    if (!penalizedInst) skip("penalized not installed")
    o4 <- ordinalNet(x, yy, standardize=F, lambdaVals=.1, alpha=.5,
                     penalizeID=c(F, rep(T, p-1)), positiveID=rep(T, p), epsIn=1e-8, epsOut=1e-8)
    p4 <- penalized(yy, penalized=-x[,-1], unpenalized=cbind(Intercept=-1, -x[,1]), model="logistic",
                    lambda1=.5*.1*n, lambda2=.5*.1*n, positive=T, standardize=F, trace=F)
    expect_equal(coef(o4), c(p4@unpenalized, p4@penalized), check.attributes=F, tolerance=1e-3)
})

test_that("Best AIC ordinal logit matches ordinalgmifs", {
    if (!ordinalgmifsInst) skip("ordinalgmifs not installed")
    o6 <- ordinalNet(x, y, standardize=F)
    og6 <- ordinal.gmifs(y ~ 1, x=names(data.frame(x)), data.frame(x), eps=.01, scale=F)
    expect_equal(coef(o6), coef(og6), check.attributes=F, tolerance=.05)
})

test_that("Elastic net multinomial logistic regression matches glmnet", {
    if (!glmnet201Inst) skip("glmnet (>= 2.0-1) not installed")
    m7 <- mirlsNet(xLS, yMat, alpha=.5, lambdaVals=.1, linkfun=multinomialLink,
                   penalizeID=rep(c(F, rep(T, p)), k), betaStart=rep(0, (p+1)*k),
                   epsIn=1e-8, epsOut=1e-8)
    g7 <- glmnet(xScaled, y, alpha=.5, lambda=.1, family="multinomial", standardize=F)
    g7_coef <- as.vector(Reduce(rbind, coef(g7)))
    shift <- coef(m7)[1] - g7_coef[1]
    g7_coef[seq(1, (p+1)*(k-1)+1, length.out=k)] <- g7_coef[seq(1, (p+1)*(k-1)+1, length.out=k)] + shift
    expect_equal(coef(m7), g7_coef, check.attributes=F, tolerance=1e-3)
})

test_that("Unpenalized multinomial logistic regression matches nnet::multinom", {
    if (!nnetInst) skip("nnet not installed")
    m8 <- mirlsNet(xLS, yMat, lambdaVals=0, linkfun=multinomialLink,
                   penalizeID=rep(c(F, rep(T, p)), k), betaStart=rep(0, (p+1)*k),
                   epsIn=1e-8, epsOut=1e-8)
    mu8 <- multinom(y ~ xScaled)
    m8_coef <- coef(m8)[-(1:(p+1))] - rep(coef(m8)[1:(p+1)], k-1)
    expect_equal(m8_coef, c(t(coef(mu8))), check.attributes=F, tolerance=1e-3)
})
