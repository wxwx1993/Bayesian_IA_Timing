#' Optimizing Interim Analysis Timing for Bayesian Adaptive Commensurate Designs
#'
#' The function to find the optimized timing for the Bayesian adaptive commensurate designs.
#' @export
#' @param historical The historical (adult) data used for borrowing information under commensurate designs.
#' @param w The prespecified expert opinion weight in the payoff function that trades off the competing decisions, between 0 and 1, default 0.5.
#' @param IA_num All possible interim analysis timings to choose from, maintain balanced arm at IA.
#' @param N The maximum total sample size if the current trial runs to completion.
#' @param p_l Early futility rule, that is the required lowest probability at the interim look,
#' that the novel treatment arm is better than prespecified minimally efficacy \code{delta_min} to declare the trial continue
#' @param p_0 Final winner rule, that is the required lowest probability at the trial completion,
#' that the novel treatment arm is better than the control arm to declare the trial succeed
#' @param p_u Early winner rule, that is the required lowest probability at the interim look,
#' that the novel treatment arm is better than the control arm to declare the trial succeed.
#' @param delta_0 The efficacy under null hypohtesis, default 0.
#' @param delta_1 The efficacy under targeted alternative hypohtesis.
#' @param delta_design The efficacy under a selected design prior.
#' @param delta_min The prespecified minimally efficacy to declare the trial continue in early futility rule.
#' @param type1_level The targeted level for the calibrated Type I error, default 0.05.
#' @param sim_sd The expected standard deviation under a selected design prior.
#' @param niter The number of interations for the Bayesian design MCMC steps, default 5000.
#' @param nrep The number of repeated calibration for the type I errors and powers calculation, default 5000.
#' @param mc.cores The number of computing cores used, default 6
#' @return The function returns a list contains:
#'
#' \code{result} The optimized timing of interim analysis under the Bayesian adaptive commensurate designs.
#'
#' \code{data} Additional data table for all calibrated designs under all possible IA timings for further study;
#' the table contains \code{PP1} and \code{PP2} are the marginal probabilities of early futility under null hypothesis and early efficacy under alternative hypothesis;
#' \code{cost} and \code{benefit} are the denominator and numerator of payoff function and \code{payoff} is the corresponding value of payoff. We consider two types of payoffs;
#' \code{var_cost} is the variance of costs;
#' \code{IA_stop} is the overall probability of stopping at IA;
#' \code{p1fin} and \code{p2fin} are the calibrated overall type I errors and powers;
#' \code{ess1} and \code{ess2} are the expected historical sample sizes for two arms respectively;
#' \code{PP1_Bayes} and \code{PP2_Bayes} are the marginal probabilities of early futility and early efficacy under selected design prior.)

optimIA <-function(historical,
                   w=0.5,
                   IA_num,
                   N,
                   p_l,
                   p_0,
                   p_u,
                   delta_0,
                   delta_1,
                   delta_design,
                   delta_min,
                   type1_level=0.05,
                   sim_sd,
                   niter=5000,
                   nrep=5000,
                   mc.cores=6){


  searchcriterion <- function(m_IA,
                              p_l=p_l,
                              p_0=p_0,
                              p_u=p_u){
    #p_u up type I down, p_l type I down, p_0 up type I down
    TypeI<-mcmapply(function(i){finalrun(sample=i,p_u=p_u,p_l=p_l,p_0=p_0,delta=delta_0,delta_min=delta_min,m=m_IA,N=N,payoff=FALSE,sim_sd=sim_sd, historical=historical,niter=niter)},1:nrep,mc.cores = mc.cores)

    Power<-mcmapply(function(i){finalrun(sample=i,p_u=p_u,p_l=p_l,p_0=p_0,delta=delta_1,delta_min=delta_min,m=m_IA,N=N,payoff=FALSE,sim_sd=sim_sd, historical=historical,niter=niter)},1:nrep,mc.cores = mc.cores)

    Payoff1<-mcmapply(function(i){finalrun(sample=i,p_u=p_u,p_l=p_l,p_0=p_0,delta=delta_design,delta_min=delta_min,m=m_IA,N=N,payoff=TRUE,sim_sd=sim_sd, historical=historical,niter=niter)},1:nrep,mc.cores = mc.cores)

    ##############PP1 and PP2 only related to p_u and p_l
    PP1 = mean(TypeI[3,]) ###########IA claim futility under H0
    PP2 = mean(Power[2,]) ###########IA success under Ha

    cost<-(2*m_IA*mean(Payoff1[4,])+N*(1-mean(Payoff1[4,])))
    var_cost<-var(2*m_IA*(Payoff1[4,])+N*(1-(Payoff1[4,])))

    benefit<-(w*PP1+(1-w)*PP2)
    benefit2<-(w*mean(Payoff1[3,])+(1-w)*mean(Payoff1[2,]))
    Payoff_vals<-benefit/cost*N
    Payoff_vals2<-benefit2/cost*N

    output<-c(w,2*m_IA,cost,PP1,PP2,benefit,benefit2,Payoff_vals,Payoff_vals2,mean(Payoff1[4,]),mean(TypeI[1,]),mean(Power[1,]),p_u,p_l,p_0,mean(Payoff1[5,]),mean(Payoff1[6,]),var_cost,mean(Payoff1[3,]),mean(Payoff1[2,]))
    names(output) = c("w","IA_num","cost","PP1","PP2","benefit","benefit2","payoff","payoff2","IA_stop","p1fin","p2fin","pu","pl","p0","ess1","ess2","var_cost","PP1_Bayes","PP2_Bayes")
    return(output)
  }
args <-expand.grid(p_l = p_l, p_0 = p_0, p_u = p_u)
pu25<-lapply(IA_num/2, function(m_IA){
  t(mapply(searchcriterion,p_l=args$p_l, p_0=args$p_0, p_u =args$p_u, m_IA=m_IA))
})

pu25.optimal<-data.frame(t(sapply(1:length(pu25),function(i){pu25[[i]][which.min(abs(data.frame(pu25[[i]])$p1fin-type1_level)),]})))
pu25.optimal<-pu25.optimal[order(pu25.optimal$w, pu25.optimal$IA_num),]

optimIA = pu25.optimal[which.max(pu25.optimal$payoff),]$IA_num

return(list(
            result = paste0("The interim analysis timing is optimzied when ", optimIA, " (", percent(optimIA/N),") patients are enrolled." ),
            data = pu25.optimal))
}


designgen<-function(historical,
                    sample,
                    N,
                    m,
                    delta,
                    sim_sd){
  set.seed(sample)
  ###########make simple mean 0 for current controls
  y_ctrl = rnorm(N/2, mean = 0, sd = sim_sd) ##################new study control

  theta = rnorm(N/2, mean = delta ,sd = 1)
  y_case = sapply(1:(N/2), function(i){rnorm(1, theta[i],sd = sim_sd)})

  allonly = data.frame(y=matrix(c(y_ctrl, y_case), ncol=1), k = c(rep(1,N/2), rep(2,N/2)), hist=rep(0,N))
  IAonly = allonly[c(1:m,(N/2+1):(N/2+m)),]
  IAcomb = rbind(IAonly, historical)
  allcomb = rbind(allonly, historical)

  return(list(IAonly, allonly, IAcomb, allcomb))
}

datgen<-function(historical,
                 sample,
                 N,
                 m,
                 delta,
                 sim_sd){
  set.seed(sample)
  y_ctrl = rnorm(N/2, mean = 0, sd = sim_sd) ##################new study control

  y_case = rnorm(N/2, mean = delta, sd = sim_sd)

  allonly = data.frame(y = matrix(c(y_ctrl, y_case),ncol=1),k = c(rep(1,N/2), rep(2,N/2)), hist = rep(0, N))
  IAonly = allonly[c(1:m, (N/2+1):(N/2+m)),] ###########if subj <= &m or subj > &N-&m then output; * select pts for IA analysis m/arm;
  IAcomb = rbind(IAonly, historical)
  allcomb = rbind(allonly, historical)

  return(list(IAonly, allonly, IAcomb, allcomb))
}

noborrow<-function(IAonly,
                   niter){
  dataTemp<- list("Y"=IAonly$y,"k"=IAonly$k,"N"=nrow(IAonly))

  result_noborrow<-rstan::sampling(stanmodels$commensurate_noborrow, data = dataTemp, iter = niter, warmup = 0.2*niter, thin = 1,chains = 1,
                                   show_messages = FALSE,verbose=FALSE)

  summary(result_noborrow)

  return(list(var10=rstan::summary(result_noborrow)$summary[4,3]^2,var20=rstan::summary(result_noborrow)$summary[5,3]^2))
}

commensurate_borrow<-function(IAcomb,
                              time,
                              p_u,
                              p_l,
                              delta_min,
                              niter){
  dataTemp<- list("Y"=IAcomb$y,"hist"=IAcomb$hist,"k"=IAcomb$k,"N"=nrow(IAcomb))

  result_commensurate<-rstan::sampling(stanmodels$commensurate, data = dataTemp, iter = niter, warmup = 0.2*niter, thin = 1,
                                       chains = 1, show_messages = FALSE,verbose=FALSE)

  summary(result_commensurate)

  pos = as.numeric(rstan::extract(result_commensurate)[["mudiff"]] > 0)
  min_eff = as.numeric(rstan::extract(result_commensurate)[["mu2"]] > delta_min)

  success = as.numeric(mean(pos) > p_u)
  #fail = min(mean(min_eff) < p_l,1-success)
  #index = c("IA","FA")[time]
  status = c(success,min(mean(min_eff) < p_l,1-success),c("IA","FA")[time], success + (min(mean(min_eff) < p_l,1-success)))
  names(status) = c("success","fail","index","overlap")

  return(list(var1=rstan::summary(result_commensurate)$summary[5,3]^2, var2 = rstan::summary(result_commensurate)$summary[6,3]^2,
              status))

}

finalrun<-function(historical,
                   sample,
                   N,
                   m,
                   p_u,
                   p_l,
                   p_0,
                   delta,
                   delta_min,
                   payoff=FALSE,
                   sim_sd,
                   niter){
  if (payoff){
    design_data<-designgen(historical=adult, sample,N, delta=delta, m,sim_sd = sim_sd)
    IAonly<-design_data[[1]]
    IAcomb<-design_data[[3]]
    allcomb<-design_data[[4]]
  }else{
    dat_data<-datgen(historical, sample,N, delta=delta, m, sim_sd = sim_sd)
    IAonly<-dat_data[[1]]
    IAcomb<-dat_data[[3]]
    allcomb<-dat_data[[4]]
  }
  noborrow_result<-noborrow(IAonly,niter=niter)
  var10<-noborrow_result[[1]]
  var20<-noborrow_result[[2]]

  borrow_result<-commensurate_borrow(IAcomb, time = 1, p_u=p_u, p_l=p_l,delta_min=delta_min, niter=niter)
  var1<-borrow_result[[1]]
  var2<-borrow_result[[2]]
  IAout<-borrow_result[[3]]

  ess1<- min(nrow(subset(historical, k==1)),max(nrow(subset(historical, k==1))*(var10/var1-1),0))
  ess2<- min(nrow(subset(historical, k==2)),max(nrow(subset(historical, k==2))*(var20/var2-1),0))

  IAout2<-commensurate_borrow(allcomb, time = 2, p_u=p_0, p_l=p_0,delta_min=delta_min, niter=niter)[[3]]

  success = IAout2["success"]
  IAstop = as.numeric( IAout["success"]==1 | IAout["fail"] ==1)
  win = as.numeric(IAout["success"]==1 | (IAstop==0 & IAout2["success"]==1 ))
  overlap = IAout2["overlap"]

  return(as.numeric(c(win, IAout["success"], IAout["fail"], IAstop, ess1, ess2, overlap)))
}
