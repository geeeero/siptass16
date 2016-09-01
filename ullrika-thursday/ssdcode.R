#############################################
## Useful functions for the SSD exercises 
## Ullrika Sahlin 2016
############################################



###############################
## functions in this file:
# HPD_func - derives highest posterior density based on a sample from the posterior
# generate_data  - generates artifical toxicity data to use in hazard assessment
# ppcheck_plot_toxicity  - plots some results from the standard SSD model
# E_LINEX_loss  - Loss function for the regulatory hazardous concentration in light of a SSD analysis
# find_LHCp_hat - function to find the Bayes optimal hazardous concentration
# generate_training_data  - generate training data for a simple regression
# ppcheck_plot_regression - plots some resluts from the regression

#############################
## simple function to calculate highest posterior density intervals
## taken from http://stats.stackexchange.com/questions/18533/find-probability-density-intervals
HPD_func <- function(x,conf=0.95){
  dx <- density(x)
  dn <- cumsum(dx$y)/sum(dx$y)
  li <- which(dn>=(1-conf)/2)[1]
  ui <- which(dn>= (1-(1-conf)/2))[1]
  dx$x[c(li,ui)]
}

#####################################
## Code to generate artificial SSD_data
## K is the number of species
## s_sizes is the possible number of replicates for each species, sampled with probability s_prob
## mu and sigma are the parameters of the true SSD
## uncertainty in estimates of logEC50 (toxicity) are generated in the code
## this uncertainty is set to be smaller than sigma and can be more similar within a species that between

generate_data <- function(mu = 2,sigma = 1,K = 4,s_sizes=1:3,
                          seed = 1988){
  s_prob=rep(1,length(s_sizes))
  set.seed(seed)
  alpha = rnorm(K,0,sigma)
  t = mu + alpha #logEC50 values for K species
  s = sample(s_sizes,K,replace = TRUE,prob = s_prob) # number of replicated for each species
  
  y = vector('list',K)
  sigma_y = vector('list',K)
  for(i in 1:K){
    #species specific uncertainty x indivdual study uncertainty x SSD sigma
    sigma_y[[i]] = runif(1)*runif(s[i])*sigma 
    y[[i]] = rnorm(s[i],t[i],sigma_y[[i]])
  }
  
  return(list(y=y,sigma_y=sigma_y,s=s,
              alpha=alpha,mu=mu,sigma=sigma,
              K=K,s_sizes=s_sizes,s_prob=s_prob,seed=seed))
}

#####################################
## Plot the SSD when
## stanmodel is model_ssd
plot_SSD <- function(stanmodel,dat,how='pdf'){
  ext <- extract(stanmodel)
  xx <- seq(mean(ext$mu)-3*mean(ext$sigma),mean(ext$mu)+3*mean(ext$sigma),by=0.01)#possible values of LHCp_hat
  plot(xx,pnorm(xx,mean(ext$mu),mean(ext$sigma)),
       xlab='LogEC50',ylab='cdf',main='SSD with species toxicity information',
       type='l',lwd=2)
  col_ind <- unlist(lapply(1:length(dat$s),function(i,s) rep(i,s[i]),dat$s))
  
  if(!is.null(dat$sigma_y)){
    if(how=='pdf'){
      pdf_scale <- max(dnorm(xx,dat$y[1],mean(dat$sigma_y)))*1.3
      for(i in 1:length(dat$y)){
        lines(xx,dnorm(xx,dat$y[i],dat$sigma_y[i])/pdf_scale,col=col_ind[i])
      }
    }else{
      for(i in 1:length(dat$y)){
        lines(xx,pnorm(xx,dat$y[i],dat$sigma_y[i]),col=col_ind[i])
      }
    }}
  else{
    points(dat$y,rep(0,dat$N),col=col_ind)
  }
  
}


#####################################
## Quick plot for posterior predictive check and compare to true values
## stanmode$l is model_1, model_2 or model_3
ppcheck_plot_toxicity <- function(stanmodel,ssd_data){
  ext <- extract(stanmodel)
  t_pred <- ext$alpha
  for(i in 1:ssd_data$K){
    t_pred[,i] <- ext$mu+ext$alpha[,i]
  }
  true_toxicity <- ssd_data$mu+ssd_data$alpha
  est_toxicity <- colMeans(t_pred)
  obs_toxicity <- unlist(ssd_data$y)
  replicated_obs_toxicity <- colMeans(ext$y_rep)
  
  plot(true_toxicity, est_toxicity, col=1:length(ssd_data$alpha),
       xlab = 'toxicity', ylab = 'predicted toxicity',
       xlim = range(true_toxicity,obs_toxicity),
       ylim = range(est_toxicity,replicated_obs_toxicity))
  abline(0,1)
  
  points(obs_toxicity,replicated_obs_toxicity , 
         col=unlist(lapply(1:length(ssd_data$s),function(i,s) rep(i,s[i]),ssd_data$s)),
         pch='+')
  legend('topleft',c('true','observed'),pch=c('o','+'),bty='n')
}

######################
## LINEX loss function 
## LHCp denotes the log10 of HCp
## E_LINEX_loss calculates the posterior mean of loss given that LHCp_hat = d
## ext is a posterior sample from the SSD sigma and the p-th percentile of the SSD, LHCp
## LINEX_alpha is parameter to control how we value over versus under estimation

E_LINEX_loss <- function(d,ext,LINEX_alpha){
  LHCp <- ext$LHCp
  sigma <- ext$sigma
  mean(exp(LINEX_alpha*(d-LHCp)/sigma) - LINEX_alpha*(d-LHCp)/sigma - 1)
}

################################################
## Find the Bayes optimal hazardous concentration and plot results
## stanmodel is model_3 or model_4
## The dat should be corresponding to the stanmodel used
## p is the target percentile in the hazard assessment

find_LHCp_hat <- function(stanmodel=model_3,LINEX_alpha = 2.5,p=0.05,
                          dat=dat_3,makeplot=FALSE,how='pdf'){
  ext <- extract(stanmodel)
  #Find LHCp_hat that minimize expected loss
  opt_out <- nlm(f = E_LINEX_loss,p = mean(ext$LHCp),
                 ext = ext,LINEX_alpha = LINEX_alpha)
  LHCp_hat <- opt_out$estimate
  
  if(makeplot){
    #Plot expected loss 
    d <- sort(c(LHCp_hat-0.02,seq(mean(ext$mu)-3*mean(ext$sigma),mean(ext$mu),by=0.01)))
    ELoss <- unlist(lapply(d,FUN=E_LINEX_loss,ext=ext,LINEX_alpha=LINEX_alpha))
    
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    xx <- sort(c(LHCp_hat-0.02,
                 seq(mean(ext$mu)-3*mean(ext$sigma),mean(ext$mu)+3*mean(ext$sigma),by=0.01)))#possible values of LHCp_hat
    plot(xx,pnorm(xx,mean(ext$mu),mean(ext$sigma)),
         xlab='LogEC50',ylab='cdf',main='SSD with species toxicity information',
         type='l',lwd=2)
    col_ind <- unlist(lapply(1:length(dat$s),function(i,s) rep(i,s[i]),dat$s))
    
    if(how=='pdf'){
      pdf_scale <- max(dnorm(xx,dat$y[1],mean(dat$sigma_y)))*1.3
      for(i in 1:length(dat$y)){
        lines(xx,dnorm(xx,dat$y[i],dat$sigma_y[i])/pdf_scale,col=col_ind[i])
      }
      if(!is.null(dat$y_reg)){
        lines(xx,dnorm(xx,mean(ext$alpha[,dat$K+1]+ext$mu),mean(ext$sigma_reg))/pdf_scale,col=max(col_ind)+1)
        legend('topright',c('QSAR prediction'),
               col=max(col_ind)+1,lty=c(1),bty='n')
      }
    }else{
      for(i in 1:length(dat$y)){
        lines(xx,pnorm(xx,dat$y[i],dat$sigma_y[i]),col=col_ind[i])
      }
      if(!is.null(dat$y_reg)){
        lines(xx,pnorm(xx,mean(ext$alpha[,dat$K+1]+ext$mu),mean(ext$sigma_reg)),col=max(col_ind)+1)
        legend('topright',c('QSAR prediction'),
               col=max(col_ind)+1,lty=c(1),bty='n')
        }
    } 
  
  abline(v = mean(ext$LHCp),col='blue',lty=2)
  abline(v = LHCp_hat,col='red',lty=2)
  
  plot(density(ext$LHCp),xlab='LogEC50',ylab='density',
       main='Uncertainty in LHCp')
  abline(v = LHCp_hat,col='red',lty=2)
  abline(v = mean(ext$LHCp),col='blue',lty=2)
  legend('topleft',c('Post mean LHC_p','LHC_p^hat'),
         col=c('blue','red'),lty=c(2,2),bty='n')
  
  plot(d,ELoss,type='l',ylab='E(loss)',xlab='LHCp_hat',
       main='Loss function')
  mtext(paste('LINEX alpha =',eval(LINEX_alpha)),3)
  abline(v = LHCp_hat,col='red',lty=2)
  abline(v = mean(ext$LHCp),col='blue',lty=2)
}

return(list(LHCp_hat=LHCp_hat,ELoss=opt_out$minimum,
            LINEX_alpha=LINEX_alpha))
}


############################
## A simple regression is added to the SSD to predict toxicity to fill data gap
## This code generates a training data of sample size n
## Toxicity is a linear function of a molecular descriptor x = 1:n
## Uncertainty in the estimates of logEC50 is considered
## x_pred is the value of the predictor that will be used to fill data gap in the SSD analysis

generate_training_data <- function(n,a,b,sigma_reg,seed=100,
                                   x_pred=7){
  set.seed(seed)
  x=runif(n,1,n)
  t_reg <- a + b*x + rnorm(n,0,sigma_reg)
  sigma_y_reg <- runif(n)*sigma_reg*2
  y_reg <- rnorm(n,t_reg,sigma_y_reg)
  return(list(n=n,y_reg=y_reg,x=x,x_pred=x_pred))
}

## plot output from sampling of simple regression
ppcheck_plot_regression <- function(stanmodel=model_reg,dat_reg){
  ext <- extract(stanmodel)
  xx <- as.list(seq(min(dat_reg$x,dat_reg$x_pred),max(dat_reg$x,dat_reg$x_pred),length.out=20))
  xL <- lapply(xx,function(x) ext$a+x*ext$b + rnorm(length(ext$a),0,ext$sigma_reg))
  m <- unlist(lapply(xL,mean))
  CI <- array(unlist(lapply(xL,HPD_func)),c(2,length(xx)))
  replicated_y <- colMeans(ext$y_reg_rep)
  #q <- apply(ext$y_reg_rep,2,'quantile',probs=c(0.025,0.975))
  q <- apply(ext$y_reg_rep,2,'HPD_func')# Highest posterior density intervals
  
  par(mfrow = c(2,1))
  plot(dat_reg$x,dat_reg$y_reg,
       xlab ='x',ylab ='y',
       ylim = range(q,CI),xlim = range(dat_reg$x_pred,dat_reg$x),
       main = 'QSAR model, simple regression')
  #points(dat_reg$x,replicated_y,col='red',pch='+')
  segments(x0=dat_reg$x,y0=q[1,],y1=q[2,],lty=1,lwd=2)
  
  lines(xx,m,col='red')
  lines(xx,CI[1,],lty=2)
  lines(xx,CI[2,],lty=2)
  points(dat_reg$x_pred,mean(ext$y_pred),col='blue',pch=3)
  lines(rep(dat_reg$x_pred,2),HPD_func(ext$y_pred),col='blue',lwd=2)
  plot(dat_reg$y_reg, replicated_y,
       xlab = 'observed', ylab = 'predicted',
       ylim = range(q))
  abline(0,1)
}



## perform sensitivity analysis

robust_hazardassessment = function(p = 0.05,LINEX_alpha = 0.15,
                        mean_mu = -5, sig_mu = 10,upper_sigma,
                        w = unlist(sapply(ssd_data$s,FUN=function(x) rep(1/x,x))),
                        ssd_data,
                        iter=5000,chains=3){
                        
#p is the largest fraction species we allow to be affected by the substance
## based on p we derive lambda which is the pth percentile of the normal distribution (which we use for the SSD)
dat = list(N=sum(ssd_data$s),K=ssd_data$K,y=unlist(ssd_data$y),s=ssd_data$s,
           sigma_y=unlist(ssd_data$sigma_y),#toxicity data
           mean_mu = mean_mu, sig_mu = sig_mu,upper_sigma = upper_sigma,#hyper priors
           w=w,
           lambda=qnorm(p))#quantile for calculating loss function

model = stan(model_name="model", model_code = code_ssd, data=dat,
             iter = iter, chains = chains, verbose = FALSE,thin=1)
find_LHCp_hat(stanmodel=model,LINEX_alpha = LINEX_alpha,p = p,
                     dat=dat,makeplot=FALSE,how='cdf')

}