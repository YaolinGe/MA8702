K <- function(p){
    return(sum(t(p) %*% p) / 2)
}

# Auxilliary variables for visualisation
x=seq(-5,5,by=0.1)
y=seq(-5,5,by=0.1)

z1 = matrix(0, ncol=101, nrow=101)
z2 = matrix(0, ncol=101, nrow=101)
z3 = matrix(0, ncol=101, nrow=101)


##############################################
# 1. Standard Gaussian with Correlation

# Parameters
Sigma = matrix(c(1,0.9,0.9,1), nrow=2)
SigmaInv = solve(Sigma)
SigmaSqrt = sqrt(det(Sigma))

# Density
myGaussian <- function(x){
    f = 1/(2*pi*SigmaSqrt)*exp(-0.5*t(x)%*%SigmaInv%*%x)
    return(f)
}

library(numDeriv)
dU <- function(x){
    return(grad(myGaussian, x))
}

x0 = c(1, 2)
dx = dU(x0)


##############################################
# 2. Multimodal

# Parameters
mu1 = c(-1.5,-1.5)
Sigma1 = matrix(c(1,0,0,1), ncol=2)
Sigma1Inv = solve(Sigma1)
Sigma1Sqrt = sqrt(det(Sigma1))

mu2 = c(1.5,1.5)
Sigma2 = matrix(c(1,0,0,1), ncol=2)
Sigma2Inv = solve(Sigma2)
Sigma2Sqrt = sqrt(det(Sigma2))

mu3 = c(-2,2)
Sigma3 = matrix(c(0.8,0,0,0.8), ncol=2)
Sigma3Inv = solve(Sigma3)
Sigma3Sqrt = sqrt(det(Sigma3))

# Density
myMultimodal <- function(x){
    f = 1/3 * 1/(2*pi*Sigma1Sqrt)*exp(-0.5*t(x-mu1)%*%Sigma1Inv%*%(x-mu1))
    f = f + 1/3 * 1/(2*pi*Sigma2Sqrt)*exp(-0.5*t(x-mu2)%*%Sigma2Inv%*%(x-mu2))
    f = f + 1/3 * 1/(2*pi*Sigma3Sqrt)*exp(-0.5*t(x-mu3)%*%Sigma3Inv%*%(x-mu3))
    return(f)
}

##############################################
# 3. Volcano

# Density
myVolcano <- function(x){
    f = 1/(2*pi)*exp(-0.5*t(x)%*%x)*(t(x)%*%x+0.25)
    return(f)
}


