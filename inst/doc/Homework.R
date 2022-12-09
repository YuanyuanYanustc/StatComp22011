## -----------------------------------------------------------------------------
y<-c(1,2,3,4,5,6)
plot(y)

## -----------------------------------------------------------------------------
plot(rnorm(10))

## -----------------------------------------------------------------------------
cat("I'm raw content.\n")

## -----------------------------------------------------------------------------
knitr::kable(head(iris))

## -----------------------------------------------------------------------------
n<-500
u<-runif(n)
x<-2/((1-u)^0.5)
hist(x, prob = TRUE, main = expression(f(x)==8*x^(-3)))
y <- seq(0, 40, 0.1)
lines(y, 8*y^(-3))

## -----------------------------------------------------------------------------
##generation function
f<-function(n,a,b,c){
  k<-0
  j<-0
  y<-numeric(n)
  while(k<n){
    u<-runif(1)
    j<-j+1
    x<-runif(1)
    if(x^(a-1)*(1-x)^(b-1)/beta(a,b)/c>=u){
      k<-k+1
      y[k]<-x
    }
  }
j
y
}

##Beta(3,2)
f(1000,3,2,3)
hist(f(1000,3,2,3), prob = TRUE)
y <- seq(0, 1, 0.01)
lines(y, dbeta(y,3,2))

## -----------------------------------------------------------------------------
n<-1000
r<-4
beta312<-2
lambda<-rgamma(n,r,beta312)
x<-rexp(n,lambda)
x
hist(x,prob=T)

## -----------------------------------------------------------------------------
n<-1000
r<-4
beta313<-2
u<-runif(n)
x<-beta313/((1-u)^(1/r))-beta313
hist(x, prob = TRUE)
y <- seq(0, 15, 0.01)#f(x)=r*(beta313^r)/(beta313+y)^(r+1)
lines(y, r*(beta313^r)/(beta313+y)^(r+1))

## -----------------------------------------------------------------------------
MC.Phi <- function(R = 10000, antithetic = TRUE) { u <- runif(R/2)
if (!antithetic) v <- runif(R/2) else v <- 1 - u
u <- c(u, v)
theta.hat<-mean(exp(u))
theta.hat }

m <- 1000
MC1 <- MC2 <- numeric(m) 
for (i in 1:m) {
MC1[i] <- MC.Phi(R = 1000, anti = FALSE)
MC2[i] <- MC.Phi(R = 1000) }

print(sd(MC1))
print(sd(MC2))
print((var(MC1) - var(MC2))/var(MC1))

## -----------------------------------------------------------------------------
m<-10000
theta.hat<-se<-numeric(2)
g<-function(x){
  1/((2*pi)^(0.5))*(x^2)*exp(-(x^2)/2)*(x>1)*(x<2)
}

x<-rweibull(m,2,2^(0.5))  #usingf1
i<-c(which(x<1),which(x>2))
x[i]<-0.5 #to catch overflow errors in g(x)
fg<-g(x)/dweibull(x,2,2^(0.5))
theta.hat[1]<-mean(fg)
se[1]<-sd(fg)

x<-rnorm(m)  #usingf2
i<-c(which(x<1),which(x>2))
x[i]<-0.5 #to catch overflow errors in g(x)
fg<-g(x)/dnorm(x)
theta.hat[2]<-mean(fg)
se[2]<-sd(fg)

rbind(theta.hat, se)

## -----------------------------------------------------------------------------
M <- 10000; k <-10 
r <- M/k #replicates per stratum
N <- 50 #number of times to repeat the estimation 
T2 <- numeric(k)
est <- matrix(0, N, 1) 
g<-function(x){
  1/((2*pi)^(0.5))*(x^2)*exp(-(x^2)/2)*(x>1)*(x<2)
}
for (i in 1:N) {
for(j in 1:k){T2[j]<-mean(g(runif(M/k,1+(j-1)/k,1+j/k)))} 
est[i, 1] <- mean(T2)
}
mean(est)
sd(est)

## -----------------------------------------------------------------------------
M<-10000
g<-function(x){
     exp(-x)/(1+x^2)
}
u<-runif(M/2)
v<-1-u
u<-c(u,v)
mean(g(u))
sd(g(u))

## -----------------------------------------------------------------------------
calcCI <- function(n, alpha) {
x <- rlnorm(n, 2, 2)
return(mean(x)+sd(x)/(n^(0.5)) * qt(1-alpha, df=n-1))
}
UCL <- replicate(1000, expr = calcCI(n = 20, alpha = .05))
sum(UCL>2)
mean(UCL>2)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y)) 
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(as.integer(max(c(outx, outy)) > 5))
}
# return 1 (reject) or 0 (do not reject H0) 


m<-10000
n<-c(20,200,2000) #small,medium,large sample size
power1<-numeric(3)
power2<-numeric(3)

# generate samples under H1 
sigma1 <- 1
sigma2 <- 1.5
for (i in 1:3) {
#estimate power of count five test
power1[i] <- mean(replicate(m, expr={ 
x <- rnorm(n[i], 0, sigma1)
y <- rnorm(n[i], 0, sigma2) 
count5test(x, y)
}))

#estimate power of F test
pvalues <- replicate(m, expr = {
x <- rnorm(n[i], 0, sigma1)
y <- rnorm(n[i], 0, sigma2) 
vartest <- var.test(x,y)
vartest$p.value } )
power2[i]<- mean(pvalues <= .05)  
}

print(power1)
print(power2)

## -----------------------------------------------------------------------------
library(boot)
lambda.mle <-function(x,ind) 
{1/mean(x[ind])}
x<-c(3,5,7,18,43,85,91,98,100,130,230,487)
obj<-boot(data=x,statistic=lambda.mle,R=10000) 
obj
round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,se=sd(obj$t)),3)

## -----------------------------------------------------------------------------
library(boot)
ci.norm<-ci.basic<-ci.perc<-ci.bca<-numeric(2) 
de <- boot(data=x,statistic=lambda.mle, R = 999)
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
ci.norm<-ci$norm[2:3]
ci.basic<-ci$basic[4:5]
ci.perc<-ci$percent[4:5]
ci.bca<-ci$bca[4:5]
print(rbind(ci.norm,ci.basic,ci.perc,ci.bca))

## -----------------------------------------------------------------------------
mu<-0
n<-1e1
m<-100
library(boot);
set.seed(12345) 
boot.mean <- function(x,i) mean(x[i]) 
ci.norm<-ci.basic<-ci.perc<-matrix(NA,m,2) 
for(i in 1:m){
R<-rnorm(n) 
de <- boot(data=R,statistic=boot.mean, R = 999)
ci <- boot.ci(de,type=c("norm","basic","perc")) 
ci.norm[i,]<-ci$norm[2:3]
ci.basic[i,]<-ci$basic[4:5]
ci.perc[i,]<-ci$percent[4:5]
}
cat('norm =',mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),'basic =',mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),'perc =',mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu),"left.norm=",sum(ci.norm[,1]>=mu),"right.norm=",sum(ci.norm[,2]<=mu),"left.basic=",sum(ci.basic[,1]>=mu),"right.basic=",sum(ci.basic[,2]<=mu),"left.perc=",sum(ci.perc[,1]>=mu),"right.perc=",sum(ci.perc[,2]<=mu))

## -----------------------------------------------------------------------------
data(patch, package = "bootstrap")
n <- nrow(patch)
cormatrix<-cov(patch[,2:6])
lambda<-eigen(cormatrix)$val
lambda1<-max(lambda)
theta.hat <-lambda1 / sum(lambda) 
print (theta.hat)

#compute the jackknife replicates, leave-one-out estimates theta.jack <- numeric(n)
theta.jack<-numeric(n)
for (i in 1:n){
  cormatrix.jack<-cov(patch[-i,2:6])
  lambda.jack<-eigen(cormatrix.jack)$val
  lambda1.jack<-max(lambda.jack)
  theta.jack[i] <-lambda1.jack / sum(lambda.jack) 
}

bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack<-sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack),3)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag 
e1 <- e2 <- e3 <- e4 <- matrix(0,(choose(n,2)),2)
m<-combn(n,2)
# for n-fold cross validation
# fit models on leave-two-out samples 
for (k in 1:choose(n,2)) {
y <- magnetic[-m[,k]] 
x <- chemical[-m[,k]]

J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[m[,k]] 
e1[k,] <- magnetic[m[,k]] - yhat1

J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[m[,k]] +J2$coef[3] * chemical[m[,k]]^2
e2[k,] <- magnetic[m[,k]] -yhat2
  
J3 <- lm(log(y) ~ x) 
logyhat3 <- J3$coef[1]+ J3$coef[2] * chemical[m[,k]] 
yhat3 <- exp(logyhat3) 
e3[k,] <- magnetic[m[,k]] -yhat3

J4 <- lm(log(y) ~ log(x))
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[m[,k]]) 
yhat4 <- exp(logyhat4)
e4[k,] <- magnetic[m[,k]] - yhat4
}

c(mean(rbind(e1[,1],e1[,2])),mean(rbind(e2[,1],e2[,2])),mean(rbind(e3[,1],e3[,2])),mean(rbind(e4[,1],e4[,2])))

#plot
a <- seq(10, 40, .1) #sequence for plotting fits
L2 <- lm(magnetic ~ chemical + I(chemical^2)) 
plot(chemical, magnetic, main="Quadratic", pch=16) 
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2 
lines(a, yhat2, lwd=2)

summary(L2)

plot(L2$fit, L2$res) 
abline(0, 0) 
qqnorm(L2$res) 
qqline(L2$res) 


## -----------------------------------------------------------------------------
x<-rnorm(20,1,2)
y<-rnorm(20,3,7)
R <- 999
z <- c(x, y)
K <- 1:40
n<-length(x) 
set.seed(12345)
reps <- numeric(R)
t0 <- cor.test(x, y,method = "spearman")$statistic 
for (i in 1:R) {
k <- sample(K, size = n, replace = FALSE) 
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- cor.test(x1, y1,method = "spearman")$statistic}
p <- mean(abs(c(t0, reps)) >= abs(t0)) 
round(c(p,cor.test(x,y)$p.value),3)

## -----------------------------------------------------------------------------
##the code below is from book
rw.Metropolis <- function(sigma, x0, N) { 
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= ((0.5*exp(-abs(y)))/(0.5*exp(-abs(x[i-1]))))) 
        x[i]<-y else{
          x[i] <- x[i-1]
          k <- k + 1 
          }
    }

return(list(x=x, k=k)) 
}

N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 1
rw1 <- rw.Metropolis(sigma[1], x0, N) 
rw2 <- rw.Metropolis(sigma[2], x0, N) 
rw3 <- rw.Metropolis(sigma[3], x0, N) 
rw4 <- rw.Metropolis(sigma[4], x0, N) 
#number of candidate points rejected
print(c(rw1$k, rw2$k, rw3$k, rw4$k))

## -----------------------------------------------------------------------------
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
    }

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j]) 
# for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  psi.w <- apply(psi, 1, "var") 
  W <- mean(psi.w)
  v.hat <- W*(n-1)/n + (B/n)
  r.hat <- v.hat / W 
  return(r.hat)
  }

sigma<-2
k <- 4
n <- 15000 
b <- 1000 #burn in length
x0 <- c(-10, -5, 5, 10)
#generate the chain
X <- matrix(0, nrow=k, ncol=n) 
for (i in 1:k) X[i, ]<-rw.Metropolis(sigma, x0[i], n)$x

#compute diagnostic statistics 
psi <- t(apply(X, 1, cumsum)) 
for (i in 1:nrow(psi)) 
  psi[i,] <- psi[i,] / (1:ncol(psi)) 
print(Gelman.Rubin(psi))

#plot psi for the four chains par(mfrow=c(2,2))
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics 
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j]) 
plot(rhat[(b+1):n], type="l", xlab="", ylab="R") 
abline(h=1.2, lty=2)


## -----------------------------------------------------------------------------
##the code below is from book
#initialize constants and parameters
N <- 5000
burn <- 1000
X <- matrix(0, N, 2)

rho <- 0.9
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1 
s2 <- sqrt(1-rho^2)*sigma2

###### generate the chain #####
X[1, ] <- c(mu1, mu2)

for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2 
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1 
  X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1 
x <- X[b:N, ]

#plot the generated sample after discarding a suitable burn-in sample
plot(x, main="", cex=.5, xlab=bquote(X), ylab=bquote(Y), ylim=range(x[,2]))

#Fit a simple linear regression model 
fit<-lm(X[,2]~X[,1])
summary(fit)
plot(x, main="", cex=.5, xlab=bquote(X), ylab=bquote(Y), ylim=range(x[,2])) #画散点图
abline(fit) # 添加回归线

# check the residuals of the model for normality and constant
hist(fit$residuals)   #for normality
qqnorm(fit$residuals)
qqline(fit$residuals)

plot(X[,1],fit$residuals) #for constant

## -----------------------------------------------------------------------------
data.gene<-function(a,b,N){
  set.seed(12345)
  e.m<-rnorm(N)
  e.y<-rnorm(N)
  a.M<-1
  a.y<-2
  x<-rnorm(N,1,2)
  M<-a.M+a*x+e.m
  y<-a.y+b*M+1*x+e.y
  return(cbind(x,M,y))
}

#h=1,2,3分别表示在（1）alpha=0（2）beta=0（3）alpha=0&beta=0原假设下的置换检验
permutation<-function(x,M,y,h){
  fit1<-lm(M~x)
  fit2<-lm(y~M+x)
  alpha.hat<-fit1$coefficients[2]
  beta.hat<--fit2$coefficients[2]
  sab<-((alpha.hat^2)*(summary(fit2)$coefficients[2,2]^2)+(beta.hat^2)*(summary(fit1)$coefficients[2,2]^2))^0.5
  t0<-alpha.hat*beta.hat/sab
  R <- 999
  reps <- numeric(R)
  
  set.seed(12345)
  for(i in 1:R){
     if(h==1){
      z<-c(x,M)
      k <- sample(length(z), size = N, replace = FALSE) 
      x1 <- z[k]
      M1 <- z[-k]
      y1<-y
      }else if(h==2){
        z<-c(y,M)
         k <- sample(length(z), size = N, replace = FALSE) 
         x1 <- x
         y1<- z[k] 
         M1 <- z[-k]
         }else{
           z<-c(x,y,M)
           k <- sample(3*N,replace = FALSE) 
           x1 <- z[k[1:100]]
           y1<- z[k[101:200]]
           M1 <- z[k[201:300]]
         }
  fit1<-lm(M1~x1)
  fit2<-lm(y1~M1+x1)
  alpha.hat<-fit1$coefficients[2]
  beta.hat<--fit2$coefficients[2]
  sab<-((alpha.hat^2)*(summary(fit2)$coefficients[2,2]^2)+(beta.hat^2)*(summary(fit1)$coefficients[2,2]^2))^0.5
  reps[i]<-alpha.hat*beta.hat/sab  
   }
  p<-mean(abs(c(t0,reps))>=abs(t0)) 
  return(p)
  }

## -----------------------------------------------------------------------------
N<-100
A<-data.gene(0,0,N)
B<-data.gene(0,1,N)
C<-data.gene(1,0,N)
round(c(permutation(A[,1],A[,2],A[,3],1),permutation(B[,1],B[,2],B[,3],1),permutation(C[,1],C[,2],C[,3],1)),3)

## -----------------------------------------------------------------------------
N<-100
A<-data.gene(0,0,N)
B<-data.gene(0,1,N)
C<-data.gene(1,0,N)
round(c(permutation(A[,1],A[,2],A[,3],2),permutation(B[,1],B[,2],B[,3],2),permutation(C[,1],C[,2],C[,3],2)),3)

## -----------------------------------------------------------------------------
N<-100
A<-data.gene(0,0,N)
B<-data.gene(0,1,N)
C<-data.gene(1,0,N)
round(c(permutation(A[,1],A[,2],A[,3],3),permutation(B[,1],B[,2],B[,3],3),permutation(C[,1],C[,2],C[,3],3)),3)

## -----------------------------------------------------------------------------
#极大似然法
u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)
n<-length(u)
mlogL <- function(lambda=1) {
# minus log-likelihood
  l<-0
  for(i in 1:n){
    l<-l+((u[i]*exp(-lambda*u[i])-v[i]*exp(-lambda*v[i]))/(exp(-lambda*v[i])-exp(-lambda*u[i])))
  }
return(l) }
library(stats4)
fit <- mle(mlogL)
as.numeric(fit@coef)

## -----------------------------------------------------------------------------
# atomic
a <- 1:5
b <- letters[1:5]
# list
df <-data.frame(a,b)

# as.vector作用于list无效
vdf <- as.vector(df)
attributes(vdf) # 属性没有改变
is.vector(vdf) # FALSE
vdf

## -----------------------------------------------------------------------------
# 转化数据框后向量带有名字
udf <- unlist(df)
attributes(udf) # 只有一个names属性
udf
vudf <- as.vector(udf)
attributes(vudf) # NULL

## -----------------------------------------------------------------------------
as.numeric("1" )
as.numeric(FALSE )
as.numeric("one" )

## -----------------------------------------------------------------------------
a<-rep(1:9)
dim(a)<-c(3,3)
a

## -----------------------------------------------------------------------------
is.matrix(a)
is.array(a)

## -----------------------------------------------------------------------------
df

## -----------------------------------------------------------------------------
#example1
as.matrix(df)
#example2
b2<-c(TRUE,F,T,F,F)
a<-1:5
df2<-data.frame(a,b2)
as.matrix(df2)

## -----------------------------------------------------------------------------
df<-df[-a,]
df

## -----------------------------------------------------------------------------
scale1<- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}
data(iris)
iris
for(i in 1:(dim(iris)[2])){
  if(is.numeric(iris[,i])) iris[,i]<-scale1(iris[,i]) 
} 
  
iris

## -----------------------------------------------------------------------------
#a)
vapply(iris[,-5], sd,numeric(1))
#b)
vapply(iris, is.numeric, logical(1))
vapply(iris[,vapply(iris, is.numeric, logical(1))], sd, numeric(1))

## -----------------------------------------------------------------------------
##the code below is from book
#initialize constants and parameters
N <- 5000
burn <- 1000
X <- matrix(0, N, 2)

rho <- 0.9
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1 
s2 <- sqrt(1-rho^2)*sigma2

###### generate the chain #####
X[1, ] <- c(mu1, mu2)

for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2 
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1 
  X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1 
x <- X[b:N, ]

qqplot(x[,1],x[,2])


