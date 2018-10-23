# Clayton-Weibull model full likelihood function
# p include Weibull baseline parameters (lambda, k), regression parameter(s) (beta), copula parameter (eta)
# x1, x2 are design matrices of two margins
# status1, status2: censoring status for two margins, with 0 for right censor and 1 for event
# var_list: a vector of covariates to be fitted in the model
# return -loglik

lik <-function(p, x1, x2,status1,status2, var_list)
{
    
    lambda <- p[1]
    k <- p[2]
    beta<- p[3:(length(p)-1)] # coefficients
    eta <- p[length(p)]
    
    t1<-indata1[,"event_time"]
    t2<-indata2[,"event_time"]
    x1<-as.matrix(indata1[,var_list])
    x2<-as.matrix(indata2[,var_list])
    
    u1<-exp(-(t1/lambda)^k*exp(x1%*%beta))
    u2<-exp(-(t2/lambda)^k*exp(x2%*%beta))
    u1_t1<- -u1*k*(1/lambda)*(t1/lambda)^(k-1)*exp(x1%*%beta)
    u2_t2<- -u2*k*(1/lambda)*(t2/lambda)^(k-1)*exp(x2%*%beta)
    
    c_val<-clt_f(u1,u2,eta)
    c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
    
    
    term1 <- log(C_val)
    term1 <- ifelse((status1 == 0) & (status2 == 0), term1, 0)
    
    term2 <- log(c_u1_val)+log(-u1_t1)
    term2 <- ifelse((status1 == 1) & (status2 == 0), term2, 0)
    
    term3 <- log(c_u2_val)+log(-u2_t2)
    term3 <- ifelse((status1 == 0) & (status2 == 1), term3, 0)
    
    term4 <- log(c_val)+log(-u1_t1)+log(-u2_t2)
    term4 <- ifelse((status1 == 1) & (status2 == 1), term4, 0)
    
    logL<-sum( term1 + term2 + term3 + term4 )
    return(-1*logL)
}
