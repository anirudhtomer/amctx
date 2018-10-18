library(MASS)
library(splines)

n <- 300 # number of subjects
K <- 10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################

# parameters for the linear mixed effects model
betas <- c("Group0" = 3.3493, "Group1" = 2.8444, "Group0:Time1" = 1.5125,
    "Group1:Time1" = 1.2533, "Group0:Time2" = 2.4739, "Group1:Time2" = 2.2994,
    "Group0:Time3" = 2.5379, "Group1:Time3" = 1.8489)
sigma.y <- 0.8034 # measurement error standard deviation


# parameters for the survival model
gammas <- c("(Intercept)" = -5.7296, "Group" = 0.4092) # coefficients for baseline covariates
alpha <- 0.3917 # association parameter
phi <- 1.6458 # shape for the Weibull baseline hazard
mean.Cens <- 12 # mean of the exponential distribution for the censoring mechanism

D <- diag(c(0.6968, 2.1259, 1.5260, 1.2329)^2)

################################################

Bkn <- c(0, 19.5)
kn <- c(2.1, 5.5)

# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
DF <- data.frame(year = times, drug = factor(rep(group, each = K)))
X <- model.matrix(~ 0 + drug + drug:ns(year, knots = kn, Boundary.knots = Bkn), data = DF)
Z <- model.matrix(~ ns(year, knots = kn, Boundary.knots = Bkn), data = DF)

# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Group" = group)

################################################

#simulate random effects
b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) # linear predictor
y <- rnorm(n * K, eta.y, sigma.y)

# simulate event times
eta.t <- as.vector(W %*% gammas)
invS <- function (t, u, i) {
    h <- function (s) {
        group0 <- 1 - group[i]
        group1 <- group[i]
        NS <- ns(s, knots = kn, Boundary.knots = Bkn)
        XX <- cbind(group0, group1, group0*NS[, 1], group1*NS[, 1],
            group0*NS[, 2], group1*NS[, 2], group0*NS[, 3], group1*NS[, 3])
        ZZ <- cbind(1, NS)
        f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
        exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha)
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
}
u <- runif(n)
trueTimes <- numeric(n)
for (i in 1:n) {
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
        tries <- tries - 1
        Up <- Up + 200
        Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
}
na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind]
W <- W[na.ind, , drop = FALSE]
long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]
n <- length(trueTimes)

# simulate censoring times from an exponential distribution,
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- runif(n, 0, 2 * mean.Cens)
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator


################################################


# keep the nonmissing cases, i.e., drop the longitudinal measurements
# that were taken after the observed event time for each subject.
ind <- times[long.na.ind] <= rep(Time, each = K)
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))

dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]
dat.id <- data.frame(Time = Time, event = event, group = W[, 2])
names(dat) <- c("time", "group", "id", "y", "Time", "event")

#summary(tapply(id, id, length))
#table(event)
#n
#mean(event)


# delete all unused objects
rm(y, X, Z, id, n, na.ind, long.na.ind, ind, Ctimes, Time, event, W, 
    betas, sigma.y, gammas, alpha, eta.t, eta.y, phi, mean.Cens, t.max,
    trueTimes, u, Root, invS, D, b, K, 
    times, group, i, tries, Up, Bkn, kn, DF)


#####################################################################

# Fit the corresponding joint model using package JM

library(JM)

lmeFit <- lme(y ~ 0 + group + group:ns(time, knots = c(2.1, 5.5), Boundary.knots = c(0, 19.5)),
        random = list(id = pdDiag(form = ~ ns(time, knots = c(2.1, 5.5), 
            Boundary.knots = c(0, 19.5)))), data = dat)
coxFit <- coxph(Surv(Time, event) ~ group, data = dat.id, x = TRUE)

jmFit <- jointModel(lmeFit, coxFit, timeVar = "time", method = "weibull-PH-aGH")
summary(jmFit)
