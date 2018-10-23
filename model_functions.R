# Density function of Clayton copula
# u1, u2: marginal distribution
# eta: copula parameter
clt_f<-function(u1,u2,eta)
{
    c_val<-(1+eta)*(u1*u2)^(-1-eta)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)
    return(c_val)
}

