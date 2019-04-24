library(foreign)
library(survival)
library(stats4)
#library(rjags)
library(mcmc)
library(MCMCpack)
library(cubature)
#library(R2Cuba)
#library(R2jags)
library(Hmisc)
library(gtools)
library(R.utils)
library(plyr)
library(dplyr)
library(parallel)
library(rstan)

N=40
#nsample=2000
#sim.seed=1
#delta = 100

##########Adult
#adult_ctrl<- rnorm(25,mean=0,sd=22) 
#adult_ctrl<- adult_ctrl-mean(adult_ctrl)
#adult_case<- rnorm(25,mean=30,sd=22) #*(19/18) 
#adult_case<- adult_case-mean(adult_case)+30
#adult<-data.frame(y=matrix(c(adult_ctrl,adult_case),ncol=1),k=c(rep(1,25),rep(2,25)),hist=rep(1,50))
#save(adult,file="adult.Rdata")

#designgen(&delta, &m, IAonly, allonly, IAcomb, allcomb); *delta=20, and IA at 15x2 (m=15) pts;
designgen<-function(sample,delta, m,sim_sd = 22){
  set.seed(sample)
  ###########make simple mean 0 for current contrls
  y_ctrl = rnorm(N/2,mean=0,sd=sim_sd) ##################new study control
  
  theta<-rnorm(N/2,mean=delta,sd=1)
  y_case = sapply(1:(N/2), function(i){rnorm(1,theta[i],sd=sim_sd)})
  
  allonly<-data.frame(y=matrix(c(y_ctrl,y_case),ncol=1),k=c(rep(1,N/2),rep(2,N/2)),hist=rep(0,N))
  IAonly<-allonly[c(1:m,(N/2+1):(N/2+m)),]
  IAcomb<- rbind(IAonly,adult)
  allcomb<-rbind(allonly,adult)
  
  return(list(IAonly, allonly, IAcomb, allcomb))
}

datgen<-function(sample,delta, m,sim_sd = 22){
  set.seed(sample)
  y_ctrl = rnorm(N/2,mean=0,sd=sim_sd) ##################new study control
  
  y_case = rnorm(N/2,mean=delta,sd=sim_sd) 
  
  allonly<-data.frame(y=matrix(c(y_ctrl,y_case),ncol=1),k=c(rep(1,N/2),rep(2,N/2)),hist=rep(0,N))
  IAonly<-allonly[c(1:m,(N/2+1):(N/2+m)),] ###########if subj <= &m or subj > &N-&m then output; * select pts for IA analysis m/arm;
  IAcomb<- rbind(IAonly,adult)
  allcomb<-rbind(allonly,adult)
  
  return(list(IAonly, allonly, IAcomb, allcomb))
}

noborrow<-function(IAonly){
  dataTemp<- list("Y"=IAonly$y,"k"=IAonly$k,"N"=nrow(IAonly))

  result_noborrow<-sampling(noborrow_stan, data = dataTemp, iter = niter, warmup = 0.2*niter, thin = 1,chains = 1)
              
  summary(result_noborrow) 
  
  return(list(var10=summary(result_noborrow)$summary[4,3]^2,var20=summary(result_noborrow)$summary[5,3]^2))
}

#macro sim1(dat,  cut1, cut2, time, out, out1, out2 );
#sim1(IAcomb, &pu, &pl, 1, IAout1, m1v, m2v); *run commensurate prior at IA;
#sim1(allcomb, &p0, &p0, 2, IAout2, m1v2, m2v2); *run commensurate prior at final;
sim1<-function(IAcomb,time=1, p_u=0.99,p_l=0.3){
  dataTemp<- list("Y"=IAcomb$y,"hist"=IAcomb$hist,"k"=IAcomb$k,"N"=nrow(IAcomb))
  
  result_commensurate<-sampling(commensurate_stan,data = dataTemp, iter = niter, warmup = 0.2*niter, thin = 1,
                                chains = 1)
  
  summary(result_commensurate)
  
  pos = as.numeric(rstan::extract(result_commensurate)[["mudiff"]]>0)
  min_eff = as.numeric(rstan::extract(result_commensurate)[["mu2"]]>15)
  
  success = as.numeric(mean(pos) > p_u)
  #fail = min(mean(min_eff) < p_l,1-success)
  #index = c("IA","FA")[time]
  status = c(success,min(mean(min_eff) < p_l,1-success),c("IA","FA")[time], success + (min(mean(min_eff) < p_l,1-success)))
  names(status) = c("success","fail","index","overlap")
  
  return(list(var1=summary(result_commensurate)$summary[5,3]^2,var2=summary(result_commensurate)$summary[6,3]^2,
              status))
  
  
}

finalrun<-function(sample,p_u,p_l,p_0,delta,m,payoff=FALSE,sim_sd=22){
  if (payoff){
  design_data<-designgen(sample,delta, m,sim_sd = sim_sd)
  IAonly<-design_data[[1]]
  IAcomb<-design_data[[3]]
  allcomb<-design_data[[4]]
  }else{
  dat_data<-datgen(sample,delta, m,sim_sd = sim_sd)
  IAonly<-dat_data[[1]]
  IAcomb<-dat_data[[3]]
  allcomb<-dat_data[[4]]
  
}
  noborrow_result<-noborrow(IAonly)
  var10<-noborrow_result[[1]]
  var20<-noborrow_result[[2]]
  
  sim1_result<-sim1(IAcomb,time = 1, p_u=p_u,p_l=p_l)
  var1<-sim1_result[[1]]
  var2<-sim1_result[[2]]
  IAout<-sim1_result[[3]]
  
  ess1<- min(25,max(25*(var10/var1-1),0)) ##############why not -1
  ess2<- min(25,max(25*(var20/var2-1),0))
  
  IAout2<-sim1(allcomb,time = 2, p_u=p_0,p_l=p_0)[[3]]
  
  success = IAout2["success"]
  IAstop = as.numeric( IAout["success"]==1 | IAout["fail"] ==1)
  win = as.numeric(IAout["success"]==1 | (IAstop==0 & IAout2["success"]==1 ))
  overlap = IAout2["overlap"]
  
  return(as.numeric(c(win,IAout["success"],IAout["fail"],IAstop,ess1,ess2,overlap)))
}


noborrow_sim1<-function(IAonly,time=1, p_u=0.99,p_l=0.3){
  dataTemp<- list("Y"=IAonly$y,"k"=IAonly$k,"N"=nrow(IAonly))
  
  result_noborrow<-sampling(noborrow_stan, data = dataTemp, iter = niter, warmup = 0.2*niter, thin = 1,chains = 1)
  
  summary(result_noborrow) 
  
  pos = as.numeric(rstan::extract(result_noborrow)[["mudiff"]]>0)
  min_eff = as.numeric(rstan::extract(result_noborrow)[["mu2"]]>15)
  
  success = as.numeric(mean(pos) > p_u)
  #fail = min(mean(min_eff) < p_l,1-success)
  #index = c("IA","FA")[time]
  status = c(success,min(mean(min_eff) < p_l,1-success),c("IA","FA")[time])
  names(status) = c("success","fail","index")
  
  return(list(status))
}

noborrow_run<-function(sample,p_u,p_l,p_0,delta,m,payoff=FALSE,sim_sd=22){
  if (payoff){
    design_data<-designgen(sample,delta, m,sim_sd = sim_sd)
    IAonly<-design_data[[1]]
    #IAcomb<-design_data[[3]]
    allonly<-design_data[[2]]
  }else{
    dat_data<-datgen(sample,delta, m,sim_sd = sim_sd)
    IAonly<-dat_data[[1]]
    #IAcomb<-dat_data[[3]]
    allonly<-dat_data[[2]]
    
  }
  noborrow_result<-noborrow_sim1(IAonly,time = 1, p_u=p_u,p_l=p_l)
  IAout<-noborrow_result[[1]]

  
  IAout2<-noborrow_sim1(allonly,time = 2, p_u=p_0,p_l=p_0)[[1]]
  
  success = IAout2["success"]
  IAstop = as.numeric( IAout["success"]==1 | IAout["fail"] ==1)
  win = as.numeric(IAout["success"]==1 | (IAstop==0 & IAout2["success"]==1 ))
  
  return(as.numeric(c(win,IAout["success"],IAout["fail"],IAstop)))
}

#TypeI<-mcmapply(function(i){noborrow_run(sample=i,p_u=p_u,p_l=p_l,p_0=p_0,delta=0,m=m_IA,payoff=FALSE,sim_sd=22)},1:10,mc.cores = 1)


