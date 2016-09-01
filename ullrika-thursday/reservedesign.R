#########################################################
## Exercise for the environmental risk anlaysis lecture
## Nature reserve design using PVA and info-gap analysis
## Model from Halpern et al Ecol Let
#########################################################

## Function to calculate persistences in a nature reserve consisting
## of two patches separated by a distance d km

calc_q <- function(d,beta,pc,alpha){
  #d is distance between reserves
  #1/beta is mean dispersal distance
  #pc is local extinction probability
  #alpha is the covariation in condistios between the two patches
  beta = max(0,beta)#beta must be positive
  0.5*
    (exp(-beta*d)*(2*pc-1)-(pc-1)*(2+(exp(-alpha*d)-1)*pc))+
    0.5*sqrt(
      4*(pc-1)*((exp(-beta*d)+pc-1)*(pc-1)-
                  exp(-alpha*d)*pc*(pc-exp(-beta*d)-1))+
        (2-3*pc-exp(-alpha*d)*pc*(pc-1)+pc^2+exp(-beta*d)*(2*pc-1))^2
      
    )
}

## Apply the function for given values on the parameters
q = calc_q(d = 10,beta=0.05,pc = 0.01,alpha = 0.05)

## Apply the function for a vector of distances
find_opt_and_plot <- function(beta,pc){
  d = 0:100
  q = sapply(d,calc_q,beta=0.05,pc=pc,alpha = 0.05)
  plot(d,q,type='l',ylim=c(pc,1),xlab='distance between reserves (km)',
       ylab = 'population persistence')
  
  ## Find optimal d for a given beta (i.e. u -> 0)
  opt = optimize(f=calc_q,interval=c(0,1000),beta=0.05,pc = pc,alpha = 0.05,maximum=TRUE)
  
  abline(v=opt$maximum,col='red')
  abline(h=opt$objective,col=1,lty=2)
  return(list(d=d,q=q))
}

saveforlater = find_opt_and_plot(beta=0.05,pc=0.5)
abline(h=0.65,col='blue',lwd=3)
#################################################
## ## Choose d when u is bounded by u_plus

## function to derive the interval expressing uncertainty in beta
unc_beta = function(u,beta_tilde){
  c((1-u)*beta_tilde,(1+u)*beta_tilde)
}

## calculate an interval for beta given an upper bound on u-value
persist_over_d_unc <- function(u_plus = 0.3,beta_tilde = 0.05,pc = 0.5,color = 1){
  beta_interval = unc_beta(u=u_plus,beta_tilde = beta_tilde)
  find_opt_and_plot(beta=beta_tilde,pc=pc)
  beta_val = seq(min(beta_interval),max(beta_interval),by=0.001)
  for(i in 1:length(beta_val)){
    q = sapply(d,calc_q,beta=beta_val[i],pc = pc,alpha = 0.05)
    lines(d,q,col=color)
  }
}

# Plot the persistence over distance (decision alternatives) and uncertainty in mean disperal distance
persist_over_d_unc(u_plus=0.4,beta_tilde = 0.05,pc = 0.5,color = 'black')
abline(h=0.65,col='blue',lwd=3)# a possible threshold for what is acceptable
lines(saveforlater$d,saveforlater$q,col='green',lwd=3,lty=2)#the model when no unc in mean disperal distance

#  A function to calcuate the largest uncertainty allowed and still have a performance larger than the criterion set by Q
info_gap <- function(Q = 0.65,d = 10,beta_tilde = 0.05,pc = 0.5){
  u_val = seq(0.001,1,by=0.01)
  dQ <- 1:length(u_val)
  for(j in 1:length(u_val)){
    beta_interval = unc_beta(u=u_val[j],beta_tilde = beta_tilde)
    beta_val = seq(min(beta_interval),max(beta_interval),by=0.001)
    q = beta_val
    for(i in 1:length(beta_val)){
      q[i] = calc_q(d,beta=beta_val[i],pc = pc,alpha = 0.05)
    }
    dQ[j] <- sum(min(q)>Q)
  }
  if(sum(dQ)>0){
    u_val[max(which(dQ==1))]}else{u_hat=0}
}

## Calculate robustness for different decision alternatives (distances) and performance levels
d = 1:50
u_hat = d
for(k in 1:length(d)){
  u_hat[k] = info_gap(Q = 0.65,d = d[k],beta_tilde = 0.05,pc = 0.5)
}
plot(d,u_hat,type='l')
#repeat but with another critera for what is acceptable
u_hat2 = u_hat
for(k in 1:length(d)){
  u_hat2[k] = info_gap(Q = 0.6,d = d[k],beta_tilde = 0.05,pc = 0.5)
}
# plot and compare the consequence of lowering the performance criteria
plot(d,u_hat,type='l',ylab='level of uncertainty',xlab='distance/decision',
     main='Info-gap analysis')
lines(d,u_hat2,col='red')
legend('topright',c('Q=0.65','Q=0.6'),col=c(1,2),lty=c(1,1),bty='n')
