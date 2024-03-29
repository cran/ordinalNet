context("Compare ordinalNet results against other packages.")
MASSInst <- require(MASS)
glmnetInst <- require(glmnet)
penalizedInst <- require(penalized)
VGAMInst <- require(VGAM)
rmsInst <- require(rms)

# Simulate data from cumulative logit model
set.seed(1)
n <- 100
zeta <- c(-1, 1)
beta <- c(1, 1, 0, 0, 0)
k <- length(zeta)
p <- length(beta)
x <- matrix(runif(n*p, -1, 1), nrow=n)
eta <- c(x %*% beta) + matrix(rep(zeta, each=n), ncol=k)
cdf <- exp(eta) / (1+exp(eta))
prob <- cbind(cdf, 1) - cbind(0, cdf)
y <- factor(sapply(1:n, function(i) sample(1:(k+1), 1, prob=prob[i, ])))
yy <- factor(ifelse(y==1, 1, 2))  # yy is binary

# Function to compare unpenalized parallel/nonparallel and forward/backward models against VGAM::vglm
# Note: VGAM family objects must be passed directly because substitute() is called on link function;
#       therefore, link function cannot be passed as an argument
vglmtest <- function(family, link, f1, f2, f3, f4, tol)
{
    # parallel - check coefficients in vector form
    o <- ordinalNet(x, y, lambdaVals=0, family=family, link=link,
                    reverse=FALSE, parallelTerms=TRUE, nonparallelTerms=FALSE)
    v <- vglm(as.ordered(y)~x, family=f1)
    expect_equal(coef(o), coef(v), check.attributes=FALSE, tolerance=tol)
    rm(o, v)

    # nonparallel - check coefficients in matrix form
    o <- ordinalNet(x, y, lambdaVals=0, family=family, link=link,
                    reverse=FALSE, parallelTerms=FALSE, nonparallelTerms=TRUE,
                    warn=FALSE)
    v <- vglm(as.ordered(y)~x, family=f2)
    coefo <- coef(o, matrix=TRUE)
    coefv <- coef(v, matrix=TRUE)
    expect_equal(coefo, coefv, check.attributes=FALSE, tolerance=tol)
    rm(o, v, coefo, coefv)

    # parallel, reverse - check coefficients in matrix form
    o <- ordinalNet(x, y, lambdaVals=0, family=family, link=link,
                    reverse=TRUE, parallelTerms=TRUE, nonparallelTerms=FALSE)
    v <- vglm(as.ordered(y)~x, family=f3)
    coefo <- coef(o, matrix=TRUE)
    coefv <- coef(v, matrix=TRUE)
    coefvr <- coefv[, ncol(coefv):1]
    expect_equal(coefo, coefvr, check.attributes=FALSE, tolerance=tol)
    rm(o, v, coefo, coefv, coefvr)

    # nonparallel, reverse - check coefficients in matrix form
    o <- ordinalNet(x, y, lambdaVals=0, family=family, link=link,
                    reverse=TRUE, parallelTerms=FALSE, nonparallelTerms=TRUE,
                    warn=FALSE)
    v <- vglm(as.ordered(y)~x, family=f4)
    coefo <- coef(o, matrix=TRUE)
    coefv <- coef(v, matrix=TRUE)
    coefvr <- coefv[, ncol(coefv):1]
    expect_equal(coefo, coefvr, check.attributes=FALSE, tolerance=tol)
    rm(o, v, coefo, coefv, coefvr)

    return(NULL)
}

###############################################################################
## Make multinomial matrices and link function
###############################################################################

test_that("Unpenalized cumulative logit matches MASS::polr", {
    if (!MASSInst) skip("MASS not installed")
    o <- ordinalNet(x, y, lambdaVals=0, family="cumulative", link="logit")
    p <- polr(y~x)
    expect_equal(coef(o), c(p$zeta, -p$coef), check.attributes=FALSE, tolerance=1e-4)
})

test_that("Unpenalized cumulative logit matches VGAM::vglm", {
    if (!VGAMInst) skip("VGAM not installed")
    family <- "cumulative"
    link <- "logit"
    f1 <- cumulative(link="logitlink", parallel=TRUE, reverse=FALSE)
    f2 <- cumulative(link="logitlink", parallel=FALSE, reverse=FALSE)
    f3 <- cumulative(link="logitlink", parallel=TRUE, reverse=TRUE)
    f4 <- cumulative(link="logitlink", parallel=FALSE, reverse=TRUE)
    tol <- 1e-3
    vglmtest(family, link, f1, f2, f3, f4, tol)
})

test_that("Unpenalized sratio probit matches VGAM::vglm", {
    if (!VGAMInst) skip("VGAM not installed")
    family <- "sratio"
    link <- "probit"
    f1 <- sratio(link="probitlink", parallel=TRUE, reverse=FALSE)
    f2 <- sratio(link="probitlink", parallel=FALSE, reverse=FALSE)
    f3 <- sratio(link="probitlink", parallel=TRUE, reverse=TRUE)
    f4 <- sratio(link="probitlink", parallel=FALSE, reverse=TRUE)
    tol <- 1e-3
    vglmtest(family, link, f1, f2, f3, f4, tol)
})

test_that("Unpenalized cratio cloglog matches VGAM::vglm", {
    if (!VGAMInst) skip("VGAM not installed")
    family <- "cratio"
    link <- "cloglog"
    f1 <- cratio(link="clogloglink", parallel=TRUE, reverse=FALSE)
    f2 <- cratio(link="clogloglink", parallel=FALSE, reverse=FALSE)
    f3 <- cratio(link="clogloglink", parallel=TRUE, reverse=TRUE)
    f4 <- cratio(link="clogloglink", parallel=FALSE, reverse=TRUE)
    tol <- 1e-3
    vglmtest(family, link, f1, f2, f3, f4, tol)
})

# Could not obtain VGAM::vglm convergence with any link function and acat family
# Also, VGAM::vglm column names of coefficient matrix are written with loge link instead of logit
test_that("Unpenalized acat logit/loge matches VGAM::vglm", {
    if (!VGAMInst) skip("VGAM not installed")
    family <- "acat"
    link <- "logit"
    # link="loge" is deprecated, use "loglink" instead.
    f1 <- acat(link="loglink", parallel=TRUE, reverse=FALSE)
    f2 <- acat(link="loglink", parallel=FALSE, reverse=FALSE)
    f3 <- acat(link="loglink", parallel=TRUE, reverse=TRUE)
    f4 <- acat(link="loglink", parallel=FALSE, reverse=TRUE)
    tol <- 1e-3
    vglmtest(family, link, f1, f2, f3, f4, tol)
})

test_that("Unpenalized acat logit/loge matches VGAM::vglm - multiple response values per observation", {
    if (!VGAMInst) skip("VGAM not installed")
    pneumo <- transform(pneumo, let = log(exposure.time))
    pneumox <- as.matrix(pneumo["let"])
    pneumoyMat <- as.matrix(pneumo[c("normal", "mild", "severe")])
    # this is the category order used in the VGAM::acat example: is it really correct?

    # nonparallel - check coefficients in matrix form
    o <- ordinalNet(pneumox, pneumoyMat, lambdaVals=0, family="acat", link="logit",
                    parallelTerms=FALSE, nonparallelTerms=TRUE)
    v <- vglm(cbind(normal, mild, severe) ~ let, family=acat, data=pneumo)
    coefo <- coef(o, matrix=TRUE)
    coefv <- coef(v, matrix=TRUE)
    expect_equal(coefo, coefv, check.attributes=FALSE, tolerance=1e-3)
})

test_that("Elastic net binary logistic regression matches glmnet and penalized", {
    if (!glmnetInst) skip("glmnet not installed")
    if (!penalizedInst) skip("penalized not installed")
    # Test for standardized covariates
    o <- ordinalNet(x, yy, alpha=.5, lambdaVals=.1, family="cumulative", link="logit", standardize=TRUE)
    g <- glmnet(x, yy, family="binomial", alpha=.5, lambda=.1, standardize=TRUE)
    pp <- penalized(yy, x, model="logistic", lambda1=.5*.1*n, lambda2=.5*.1*n, trace=FALSE, standardize=TRUE)
    expect_equal(coef(o), c(-g$a0, -as.vector(g$beta)), check.attributes=FALSE, tolerance=1e-4)
    expect_equal(coef(o), c(-pp@unpenalized, -pp@penalized), check.attributes=FALSE)
    rm(o, g, pp)
    # Test for unstandardized covariates
    o <- ordinalNet(x, yy, alpha=.5, lambdaVals=.1, family="cumulative", link="logit", standardize=FALSE)
    g <- glmnet(x, yy, family="binomial", alpha=.5, lambda=.1, standardize=FALSE)
    pp <- penalized(yy, x, model="logistic", lambda1=.5*.1*n, lambda2=.5*.1*n, trace=FALSE, standardize=FALSE)
    expect_equal(coef(o), c(-g$a0, -as.vector(g$beta)), check.attributes=FALSE, tolerance=1e-4)
    expect_equal(coef(o), c(-pp@unpenalized, -pp@penalized), check.attributes=FALSE)
})

test_that("Elastic net binary logistic regression with penalty factors matches glmnet", {
    if (!glmnetInst) skip("glmnet not installed")
    # Test for standardized covariates
    penaltyFactors <- seq(0, 1, length.out=p)
    penaltyFactorsStd <- penaltyFactors * p / sum(penaltyFactors)
    o <- ordinalNet(x, yy, alpha=.5, lambdaVals=.1, family="cumulative", link="logit", penaltyFactors=penaltyFactorsStd)
    g <- glmnet(x, yy, family="binomial", alpha=.5, lambda=.1, penalty.factor=penaltyFactorsStd)
    expect_equal(coef(o), c(-g$a0, -as.vector(g$beta)), check.attributes=FALSE, tolerance=1e-4)
})

test_that("Binary logistic regression with positive constraints matches penalized", {
    if (!penalizedInst) skip("penalized not installed")
    # Note: constraints are in play for X3 and X4
    positiveID <- rep(TRUE, p)
    o <- ordinalNet(x, yy, alpha=.5, lambdaVals=0, family="cumulative", link="logit", positiveID=positiveID)
    pp <- penalized(yy, -x, model="logistic", trace=FALSE, standardize=TRUE, positive=TRUE)
    expect_equal(coef(o), c(-pp@unpenalized, pp@penalized), check.attributes=FALSE)
})

test_that("Cumulative logit ridge matches rms::lrm", {
    if (!rmsInst) skip("rms not installed.")
    m1 <- ordinalNet(x, y, alpha=0, lambdaVals=.1)
    m2 <- lrm(y~., data=data.frame(y, x), penalty=.1*n)
    expect_equal(coef(m1), -coef(m2), check.attributes=FALSE, tolerance=1e-2)
})
