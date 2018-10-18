simulateJoint <- function (n = 500, n.test = 50) {
    K <- 10  # number of planned repeated measurements per subject, per outcome
    t.max <- 19.5 # maximum follow-up time
    
    ################################################
        
    # parameters for the linear mixed effects model
    betas <- c("Intercept" = 2.94, "Time1" = 1.30, "Time2" = 1.84, "Time3" = 1.82)
    sigma.y <- 0.6 # measurement error standard deviation
    
    # parameters for the survival model
    gammas <- c("(Intercept)" = -6.7, "Group" = 0.5)
    alpha <- 0.191
    Dalpha <- -1.064
    phi <- 2 #2.2
    mean.Cens <- 14
    
    D <- matrix(0, 4, 4)
    D[lower.tri(D, TRUE)] <- c(0.71, 0.33, 0.07, 1.26, 2.68, 3.81, 4.35, 7.62, 5.4, 8)
    D <- D + t(D)
    diag(D) <- diag(D) * 0.5
    
    ################################################
    
    Bkn <- c(0, 19.5)
    kn <- c(2.1, 5.5)
    
    # design matrices for the longitudinal measurement model
    times <- c(replicate(n, c(0, 0.5, 1, sort(runif(K - 3, 1, t.max))))) 
    group <- rep(0:1, each = n/2)
    DF <- data.frame(time = times)
    X <- model.matrix(~ ns(time, knots = kn, Boundary.knots = Bkn), 
                      data = DF)
    Z <- model.matrix(~ ns(time, knots = kn, Boundary.knots = Bkn), data = DF)
    
    # design matrix for the survival model
    W <- cbind("(Intercept)" = 1, "Group" = group)
    
    ################################################
    
    # simulate random effects
    b <- mvrnorm(n, rep(0, nrow(D)), D)
    
    # simulate longitudinal responses
    id <- rep(1:n, each = K)
    eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ]))
    y <- rnorm(n * K, eta.y, sigma.y)
    
    # simulate event times
    eta.t <- as.vector(W %*% gammas)
    invS <- function (t, u, i) {
        h <- function (s) {
            NS <- ns(s, knots = kn, Boundary.knots = Bkn)
            DNS <- dns(s, knots = kn, Boundary.knots = Bkn)
            XX <- cbind(1, NS)
            ZZ <- cbind(1, NS)
            XXd <- DNS
            ZZd <- DNS
            f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
            f2 <- as.vector(XXd %*% betas[2:4] + rowSums(ZZd * b[rep(i, nrow(ZZd)), 2:4]))
            exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha + f2 * Dalpha)
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
            Up <- Up + 50
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
    DF <- DF[long.na.ind, , drop = FALSE]
    n <- length(trueTimes)
    
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
    
    dat <- DF[ind, , drop = FALSE]
    dat$id <- id
    dat$y <- y
    dat$Time <- Time[id]
    dat$event <- event[id]
    dat.id <- data.frame(Time = Time, event = event, group = W[, 2])
    dat <- dat[c("id", "y", "time", "Time", "event")]
    dat$group <- dat.id$group[id]
    #summary(tapply(id, id, length))
    #table(event)
    #n
    #mean(event)
    #summary(dat.id$Time)
    #summary(dat$time)
    
    # extract 100 subjects for whom predictions will be be made
    ids <- unique(dat$id[dat$Time > 1])
    samp.id <- sample(ids, n.test)
    predDat <- dat[dat$id %in% samp.id, ]
    dat <- dat[!dat$id %in% samp.id, ]
    dat.id <- dat.id[-samp.id, ]
    
    # true values for parameters and random effects
    trueValues <- list(betas = betas, tau = 1/sigma.y^2, gammas = gammas, 
                       alphas = alpha, Dalphas = Dalpha, sigma.t = phi, inv.D = solve(D), 
                       b = b[-samp.id, , drop = FALSE])
    
    # return list
    list(DF = dat, DF.id = dat.id, predDat = predDat, 
         b = b[samp.id, , drop = FALSE], trueValues = trueValues)
}

simulateLong <- function (time, b) {
    # parameters for the linear mixed effects model
    betas <- c("Intercept" = 2.94, "Time1" = 1.30, "Time2" = 2.84, "Time3" = 2.82)
    sigma.y <- 0.6    
    # design matrices for the longitudinal measurement model
    Bkn <- c(0, 19.5)
    kn <- c(2.1, 5.5)
    DF <- data.frame(time = time)
    X <- model.matrix(~ ns(time, knots = kn, Boundary.knots = Bkn), 
                      data = DF)
    Z <- model.matrix(~ ns(time, knots = kn, Boundary.knots = Bkn), data = DF)
    # simulate
    eta.y <- as.vector(X %*% betas + Z %*% b)
    rnorm(1, eta.y, sigma.y)
}

optimalTime <- function (screening, Data, jointFit, cutoff, gap, Dt) {
    # Test data
    predDat <- Data$predDat[Data$predDat$time %in% c(0, 0.5, 1), ]
    b <- Data$b
    rownames(b) <- unique(predDat$id)
    optimal_times <- numeric(nrow(b))
    names(optimal_times) <- rownames(b)
    visit_times <- vector("list", nrow(b))
    names(visit_times) <- rownames(b)
    it <- 0
    while (nrow(predDat)) {
        survPreds <- sapply(split(predDat, factor(predDat$id)), function (d) {
            last.time <- tail(d$time, 1)
            sfit <- survfitJM(jointFit, d, survTimes = last.time + Dt, seed = 2015 + nsub)
            c(sfit$summaries[[1]][, 2], sfit$last.time)
        })
        ind <- survPreds[1, , drop = FALSE] < cutoff
        ids_out <- colnames(ind[, ind, drop = FALSE])
        if(!is.null(ids_out)) {
            optimal_times[ids_out] <- survPreds[2, ids_out]
            visit_times[ids_out] <- with(predDat[predDat$id %in% ids_out, ], split(time, id))
            predDat <- predDat[!predDat$id %in% ids_out, , drop = FALSE]
            b <- b[!rownames(b) %in% ids_out, , drop = FALSE]
        }
        if (nrow(b)) {
            predDat <- do.call(rbind, mapply(function (d, b) {
                k <- nrow(d)
                d_out <- rbind(d, d[k, ])
                if (screening == "fixed") {
                    d_out$time[k + 1] <- d_out$time[k + 1] + gap
                    d_out$y[k + 1] <- simulateLong(d_out$time[k + 1], b)            
                } else {
                    sfit <- if (max(d$time) > quantile(jointFit$y$Time, 0.9) + 0.01) {
                        sss <- seq(max(d$time), max(d$time) * 1.3, len = 35)
                        survfitJM(jointFit, newdata = d, survTimes = sss, 
                                  seed = 2015 + nsub)$summaries[[1]]
                    } else {
                        survfitJM(jointFit, newdata = d, seed = 2015 + nsub)$summaries[[1]]
                    }
                    tup <- min(gap, max(sfit[sfit[, "Mean"] > cutoff, "times"] - max(d$time)))
                    v <- dynInfo(jointFit, newdata = d, Dt = tup, K = 5, 
                                 seed = 2015 + nsub)$summary
                    tt <- v$time[which.max(v$Info)]
                    d_out$time[k + 1] <- tt
                    d_out$y[k + 1] <- simulateLong(d_out$time[k + 1], b)            
                }
                d_out
            }, d = split(predDat, factor(predDat$id)), b = split(b, row(b)), SIMPLIFY = FALSE))
        }
        it <- it + 1
        if (it > 7) {
            optimal_times[optimal_times == 0] <- NA
            break
        }
    }
    list(optimal_times = optimal_times, visit_times = visit_times)
}

true_optimal_times <- function (b, predDat, trueValues) {
    betas <- trueValues$betas
    gammas <- trueValues$gammas
    alpha <- trueValues$alphas
    Dalpha <- trueValues$Dalphas
    phi <- trueValues$sigma.t
    W <- model.matrix(~ group, data = predDat[!duplicated(predDat$id), ])
    eta.t <- c(W %*% gammas)
    Bkn <- c(0, 19.5)
    kn <- c(2.1, 5.5)
    cumHaz <- function (t) {
        h <- function (s) {
            NS <- ns(s, knots = kn, Boundary.knots = Bkn)
            DNS <- dns(s, knots = kn, Boundary.knots = Bkn)
            XX <- cbind(1, NS)
            ZZ <- cbind(1, NS)
            XXd <- DNS
            ZZd <- DNS
            f1 <- as.vector(XX %*% betas + ZZ %*% b.i)
            f2 <- as.vector(XXd %*% betas[2:4] + ZZd %*% b.i[2:4])
            exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha + f2 * Dalpha)
            exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha)
        }
        integrate(h, lower = 0, upper = t)$value + log(0.8)
    }
    n <- nrow(W)
    out <- numeric(n)
    for (i in seq_len(n)) {
        b.i <- b[i, ]
        out[i] <- uniroot(cumHaz, interval = c(1e-05, 20))$root
    }
    out
}

true_optimal_times <- function (b, predDat, trueValues, Dt, cutoff) {
    betas <- trueValues$betas
    gammas <- trueValues$gammas
    alpha <- trueValues$alphas
    Dalpha <- trueValues$Dalphas
    phi <- trueValues$sigma.t
    W <- model.matrix(~ group, data = predDat[!duplicated(predDat$id), ])
    eta.t <- c(W %*% gammas)
    Bkn <- c(0, 19.5)
    kn <- c(2.1, 5.5)
    surv <- function (t) {
        h <- function (s) {
            NS <- ns(s, knots = kn, Boundary.knots = Bkn)
            DNS <- dns(s, knots = kn, Boundary.knots = Bkn)
            XX <- cbind(1, NS)
            ZZ <- cbind(1, NS)
            XXd <- DNS
            ZZd <- DNS
            f1 <- as.vector(XX %*% betas + ZZ %*% b.i)
            f2 <- as.vector(XXd %*% betas[2:4] + ZZd %*% b.i[2:4])
            exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha + f2 * Dalpha)
            exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha)
        }
        exp(- integrate(h, lower = 0, upper = t)$value)
    }
    n <- nrow(W)
    times <- round(seq(0, 21, 0.1), 1)
    ntimes <- length(times)
    ind <- seq_along(times)
    ind2 <- seq(which(times == Dt), ntimes)
    ind3 <- cbind(ind2, seq_along(ind2))
    out <- numeric(n)
    for (i in seq_len(n)) {
        b.i <- b[i, ]
        survs <- sapply(times, surv)
        pi.u.t <- sapply(ind, function (k) survs / survs[k])
        pi.u.t[pi.u.t > 1] <- 1
        ii <- min(which.min(abs(pi.u.t[ind3] - cutoff)) + which(times == Dt), length(times))
        out[i] <- times[ii]
    }
    out
}

