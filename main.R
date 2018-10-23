#################################################################################################
### Estimation and Inference under the Clayton-Weibull model for bivariate right-censored data###
##################################################################################################
source("model_functions.R") 
source("step1b.R") 
source("lik.R") 
library("stats") 
library("pracma")
library("survival")


# data input
load("example_data.RData")
var_list= c("var1","var2")
n <- dim(data)[1]
p <- length(var_list)

# separate into two margins
indata1 <- data[data[,"ind"]==1, ]
indata2 <- data[data[,"ind"]==2, ]
status1 <- indata1$status
status2 <- indata2$status

x1 <- as.matrix(indata1[,var_list],nrow = n)
x2 <- as.matrix(indata2[,var_list],nrow = n)
x <- data.frame(id=c(indata1$id,indata2$id), event_time = c(indata1$event_time,indata2$event_time), status = c(indata1$status,indata2$status), rbind(x1,x2))



##################################################################################################
### Fitting Clayton-Weibull using covariate var1 only
##################################################################################################
#### step 1a: marginal model
M <- survreg(Surv(event_time, status) ~ var1 + cluster(id), data=x, dist="weibull")
lambda_ini <- exp(M$coef[1]) # scale
k_ini <- 1/M$scale # shape
beta_ini <- -1*coef(M)[-1]*k_ini # coefficients

####  step 1b: update eta
eta_ini <- 2
fit0<-nlm(step1b, p=eta_ini, p2 = c(lambda_ini,k_ini,beta_ini),
          x1=x1, x2=x2,status1=status1,status2=status2, var_list=c("var1"),
          iterlim=500, steptol = 1e-6)
eta_ini<-fit0$estimate

### step 2, with all starting values of p, optimize all parameters at the same time 
model_step2 <- nlm(lik, p=c(lambda_ini,k_ini,beta_ini,eta_ini),
                   x1=x1, x2=x2,status1=status1,status2=status2, var_list=c("var1"),
                   iterlim = 500, steptol = 1e-6, hessian = T)
                   
inv_info = solve(model_step2$hessian) 
se = sqrt(diag(inv_info))
beta = model_step2$estimate # contains lambda, k, beta and eta
stat = (beta-0)^2/se^2
pvalue = pchisq(stat,1,lower.tail=F)
summary = cbind(beta, se, stat, pvalue)
rownames(summary) = c("lambda","k", "var1", "eta")
colnames(summary) = c("estimate","SE","stat","pvalue")
summary
# estimate         SE      stat        pvalue
# lambda 9.80414809 0.66057049 220.28320  7.845309e-50
# k      1.93680321 0.07098497 744.45460 6.444329e-164
# var1   0.09564934 0.01956758  23.89404  1.017863e-06
# eta    0.58768426 0.12204221  23.18825  1.468920e-06


############################################################################
### Generalized score test Under Null (H0: effect of var2 =0)
############################################################################
estimates = c(model_step2$estimate[1:2], # lambda, k
              model_step2$estimate[3:(2+p-1)], # var1
              0, # var2
              model_step2$estimate[length(model_step2$estimate)]) # eta
score = grad(lik, x0 = estimates,x1=x1, x2=x2,status1=status1,status2=status2, var_list=c("var1", "var2")) # score function by numerical approximation
hes = hessian(lik, x0 = estimates,x1=x1, x2=x2,status1=status1,status2=status2, var_list=c("var1", "var2")) # hessian matrix by numeircal approximation
test_stat = t(score) %*% solve(hes) %*% score  # score test statistics
test_stat
# 0.05225796
p_value_score <- pchisq(test_stat,1, lower.tail = F) # score test p value
p_value_score
# 0.8191798  

    

