## M: number of MCMC iterations to go
## sd: standard deviation for the random walk proposal.

## Random walk Metropolis Hastings
mcmc_rw <- function(M, sd){
    
    # store samples here
    xsamples <- rep(NA, M)
    # specify a starting value
    xsamples[1] <- 0.0
    
    # store acceptance rates for sets of iterations
    plotAcc <- 10
    acceptrates <- rep(NA, (M/plotAcc))
    
    # count acceptances and rejections
    yes <- 0
    no <- 0
    
    # MH-Iteration
    for(k in 2:M){
        # value of the past iteration
        old <- xsamples[k-1]
        
        # propose a new value based on the old one
        proposal <- rnorm(1, mean=old, sd=sd)
        
        # compute acceptance ratio (Note the log-scale!!!)
        log.target.ratio <- dnorm(proposal, mean=0, sd=1, log=TRUE) - 
            dnorm(old, mean=0, sd=1, log=TRUE)
        
        # the proposal ratio is equal to 1 as we have 
        #a symmetric proposal distribution
        log.proposal.ratio <- 0
        
        # get the acceptance probability (on log scale)
        alpha <- log.target.ratio+log.proposal.ratio
        
        # accept-reject step
        if(log(runif(1)) <= alpha){	
            # accept the proposed value
            xsamples[k] <- proposal
            # increase counter of accepted values
            yes <- yes + 1
        }
        else{
            # stay with the old value
            xsamples[k] <- old
            no <- no + 1
        }
        
        # every plotAcc iterations, record the acceptance rate
        if(k%%plotAcc == 0) 
            acceptrates[k/plotAcc] <- yes/(yes+no)*100
    }
    
    # acceptance rate
    cat("The acceptance rate is: ", round(yes/(yes+no)*100,2), "%\n", sep="")
    return(list(xsamples=xsamples, acceptrates=acceptrates))
}

set.seed(03122010)
s1 <- mcmc_rw(1000, sd=0.24)
s2 <- mcmc_rw(1000, sd=2.4)
s3 <- mcmc_rw(1000, sd=24)

par(mfrow=c(3,3), mar=c(4,4,2,0.5))
plot(s1$xsam, ylim=c(-3.5,3.5), type="l", 
     main=expression(paste(sigma,"=0.24")))
plot(s2$xsam, ylim=c(-3.5,3.5), type="l", 
     main=expression(paste(sigma,"=2.4")))
plot(s3$xsam, ylim=c(-3.5,3.5), type="l", 
     main=expression(paste(sigma,"=24")))
acf(s1$xsam, lag.max=40, main="", ylab="s1 ACF")
acf(s2$xsam, lag.max=40, main="", ylab="s2 ACF")
acf(s3$xsam, lag.max=40, main="", ylab="s3 ACF")
plot(s1$acc, type="l", ylim=c(0,100))
plot(s2$acc, type="l", ylim=c(0,100))
plot(s3$acc, type="l", ylim=c(0,100))

