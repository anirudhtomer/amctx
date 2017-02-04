x = seq(0,4,0.02)
y = 60 + 6.4*x -5.3*x*x + 1.2*x^3 + rnorm(n = length(x), sd = 1)

training = sample(1:201, size = 100)
real = setdiff(1:201, training)

plot(y=y[real],x=x[real])

model = lm(y~x+I(x^2) + I(x^3))
summary(model)

model = lm(y~poly(x,3))
summary(model)


glmnet_fit = glmnet(cbind(x[training], x[training]^2, x[training]^3), y[training], family="gaussian", alpha = 1, nlambda=100000)
cv_glmnet_fit = cv.glmnet(cbind(x[training], x[training]^2, x[training]^3), y[training], nfolds=10, type.measure = "mse", alpha=1)
yhat_test=predict(glmnet_fit, cbind(x[real], x[real]^2, x[real]^3), s=cv_glmnet_fit$lambda.min)
c(sum((y[real]-yhat_test)^2), cv_glmnet_fit$lambda.min)
glmnet(cbind(rep(1, length(training)),x[training], x[training]^2, x[training]^3), y[training], family="gaussian", alpha = 1, lambda = cv_glmnet_fit$lambda.min)$beta

glmnet_fit = glmnet(poly(x[training],3), y[training], family="gaussian", alpha = 1, nlambda=100000)
cv_glmnet_fit = cv.glmnet(poly(x[training],3), y[training], nfolds=10, type.measure = "mse", alpha=1)
yhat_test=predict(glmnet_fit, poly(x[real],3), s=cv_glmnet_fit$lambda.min)
c(sum((y[real]-yhat_test)^2), cv_glmnet_fit$lambda.min)
tt=glmnet(poly(x[training],3), y[training], family="gaussian", alpha = 1, lambda = cv_glmnet_fit$lambda.min)


error = rnorm(nrow(x))
x = matrix(c(c(1,1,1,0,0,0,0),c(1,1,1,0,1,0,1)), ncol = 2,byrow = F)
y = 10 + x[,1]*5 + 6*x[,2] + error
summary(lm(y~x))


####################
# Band of residuals
####################
x = 1:100
x3 = (11:110) * 10
x2 = factor(rep(c("A","B", "C", "D"), each = 25))
x2_val = rep(c(40, 60, 100, 110), each = 25)
y = 60 + 20*x + rnorm(length(x), sd= 10)
y = y + x2_val
y = y + 10.5 * log(x3) + 6 * x2_val*log(x3)

plot(y~x)
mod = lm(y~x2 * I(log(x3)))
qplot(y=mod$residuals, mod$fitted.values, geom=c("point"))

