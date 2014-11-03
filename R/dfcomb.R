logistic_sim = function(ndose_a1, ndose_a2, p_tox, target, target_min, target_max, prior_tox_a1, prior_tox_a2, n_cohort, 
                        cohort, tite=FALSE, time_full=0, poisson_rate=0, nsim, c_e=0.85, c_d=0.45, c_stop=0.95, n_min=4, 
                        seed = 14061991){
  n_min=n_min-1                      
  c_d = 1-c_d
  dim_ptox = dim(p_tox)
  p_tox_c = t(p_tox)[,dim_ptox[1]:1]
  
  if(dim_ptox[2] != ndose_a1 || dim_ptox[1] != ndose_a2){
    stop("Wrong dimension of the matrix for true toxicity probabilities.")
  }
  n_prior_tox_a1 = length(prior_tox_a1)
  if(n_prior_tox_a1 != ndose_a1){
    stop("The entered vector of initial guessed toxicity probabities for agent 1 is of wrong length.")
  }
  n_prior_tox_a2 = length(prior_tox_a2)
  if(n_prior_tox_a2 != ndose_a2){
    stop("The entered vector of initial guessed toxicity probabities for agent 2 is of wrong length.")
  }
  
  
  ndose_a1 = as.integer(ndose_a1)[1]
  ndose_a2 = as.integer(ndose_a2)[1]
  target = as.double(target)[1]
  target_min = as.double(target_min)[1]
  target_max = as.double(target_max)[1]
  prior_tox_a1 = as.double(prior_tox_a1)
  prior_tox_a2 = as.double(prior_tox_a2)
  n_cohort = as.integer(n_cohort)[1]
  cohort = as.integer(cohort)[1]
  tite = as.logical(tite)[1]
  time_full = as.double(time_full)[1]
  poisson_rate = as.double(poisson_rate)[1]
  nsim = as.integer(nsim)[1]
  c_e = as.double(c_e)[1]
  c_d = as.double(c_d)[1]
  c_stop = as.double(c_stop)[1]
  n_min = as.integer(n_min)[1]
  seed = as.integer(seed)[1]
  
  for(a1 in 1:ndose_a1){
    if(prior_tox_a1[a1] < 0 || prior_tox_a1[a1] > 1){
      stop("At least one of the initial guessed toxicity probability for agent 1 is not comprised between 0 and 1.")
    }
  }
  for(a2 in 1:ndose_a2){
    if(prior_tox_a2[a2] < 0 || prior_tox_a2[a2] > 1){
      stop("At least one of the initial guessed toxicity probability for agent 2 is not comprised between 0 and 1.")
    }
  }
  for(a1 in 1:ndose_a1){
    for(a2 in 1:ndose_a2){
      if(p_tox_c[a1,a2] < 0 || p_tox_c[a1,a2] > 1){
        stop("At least one of the initial guessed toxicity probability is not comprised between 0 and 1.")
      }
    }
  }
  p_tox_c_na = matrix(NA, nrow=ndose_a1+1, ncol=ndose_a2+1)
  p_tox_c_na[1:ndose_a1, 1:ndose_a2] = p_tox_c
  for(a1 in 1:ndose_a1){
    for(a2 in 1:ndose_a2){
      if(p_tox_c[a1,a2] >
         min(1,p_tox_c_na[a1+1,a2],p_tox_c_na[a1,a2+1],p_tox_c_na[a1+1,a2+1],na.rm=TRUE)){
        stop("The partial ordering between true toxicity probabilities is not satisfied.")
      }
    }
  }

  p_tox_c = as.double(p_tox_c)  
  inconc = as.double(numeric(1))
  n_pat_dose = as.double(numeric(ndose_a1*ndose_a2))
  rec_dose = as.double(numeric(ndose_a1*ndose_a2))
  n_tox_dose = as.double(numeric(ndose_a1*ndose_a2))

  
  # Appeler fonction C
  logistic = .C(C_logistic_sim, tite=tite, ndose_a1=ndose_a1, ndose_a2=ndose_a2, time_full=time_full, poisson_rate=poisson_rate,
    p_tox_c=p_tox_c, target=target, target_max=target_max, target_min=target_min, prior_tox_a1=prior_tox_a1, prior_tox_a2=prior_tox_a2,
    n_cohort=n_cohort, cohort=cohort, nsim=nsim, c_e=c_e, c_d=c_d, c_stop=c_stop, n_min=n_min, seed=seed,
    rec_dose=rec_dose, n_pat_dose=n_pat_dose, n_tox_dose=n_tox_dose, inconc=inconc)

  nsim = logistic$nsim

  inconc=logistic$inconc*100
  rec_dose=logistic$rec_dose*100
  n_pat_dose=logistic$n_pat_dose
  n_tox_dose=logistic$n_tox_dose

  res = list(call = match.call(),
             tite=tite, 
             ndose_a1=ndose_a1, 
             ndose_a2=ndose_a2, 
             time_full=time_full, 
             poisson_rate=poisson_rate,
             p_tox=p_tox,
             p_tox_c=matrix(p_tox_c,nrow=ndose_a1), 
             target=target, 
             target_min=target_min,
             target_max=target_max, 
             prior_tox_a1=prior_tox_a1, 
             prior_tox_a2=prior_tox_a2,
             n_cohort=n_cohort, 
             cohort=cohort, 
             nsim=nsim, 
             c_e=c_e, 
             c_d=c_d, 
             c_stop=c_stop, 
             n_min=n_min, 
             seed=seed,
             rec_dose_c=matrix(rec_dose,nrow=ndose_a1), 
             n_pat_dose_c=matrix(n_pat_dose,nrow=ndose_a1), 
             n_tox_dose_c=matrix(n_tox_dose,nrow=ndose_a1), 
             inconc=inconc)
             
  class(res) = "logistic_sim"

  return(res)
}



print.logistic_sim = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call) 
  
  rec_dose = t(x$rec_dose_c)[x$ndose_a2:1,]
  n_pat_dose = t(x$n_pat_dose_c)[x$ndose_a2:1,]
  n_tox_dose = t(x$n_tox_dose_c)[x$ndose_a2:1,]
  p_tox = x$p_tox
  dimnames(p_tox) = list("Agent 2" = x$ndose_a2:1, "Agent 1" = 1:x$ndose_a1)
  dimnames(rec_dose) = list("Agent 2 " = x$ndose_a2:1, "Agent 1" = 1:x$ndose_a1)
  dimnames(n_pat_dose) = list("Agent 2"=x$ndose_a2:1, "Agent 1" = 1:x$ndose_a1)
  dimnames(n_tox_dose) = list("Agent 2" = x$ndose_a2:1, "Agent 1" = 1:x$ndose_a1)

  print_rnd= function (hd, x) {cat(hd, "\n"); print(round(x, digits = dgt)); cat("\n")}
  print_rnd("True toxicities:", p_tox)
  print_rnd("Percentage of Selection:", rec_dose)
  print_rnd("Number of patients:" , n_pat_dose)
  print_rnd("Number of toxicities:", n_tox_dose)

  cat(paste("Percentage of inconclusive trials:\t",x$inconc,"\n",sep=""), sep="")
  nmin = which(x$cohort*(0:x$n_cohort) >= x$n_min)[1]
  cat("The minimum number of cohorts to stop the trial is:\t", nmin-1, "\n")
  cat("\n", "\n")              
  cat("Number of simulations:\t", x$nsim, "\n")
  cat("Cohort size:\t", x$cohort, "\n")
  cat("Number of cohort planned:\t", x$n_cohort, "\n")
  cat("Total patients accrued:\t", round(sum(x$n_pat_dose),1), "\n")
  cat("Toxicity target:\t", x$target, "\n")
  cat("Targeted toxicity interval:\t [", x$target_min, ",", x$target_max, "]\n")
  cat("Prior toxicity probabilities for agent 1:\n")
  print(round(x$prior_tox_a1, digits = dgt))
  cat("Prior toxicity probabilities for agent 2:\n")
  print(round(x$prior_tox_a2, digits = dgt))
  cat("Escalation threshold:\t", x$c_e, "\n")
  cat("Deescalation threshold:\t", 1-x$c_d, "\n")
  cat("Stopping threshold:\t", x$c_stop, "\n")
  if (x$tite) {
    cat("Toxicity is a time-to-event \n")
    cat("Full follow-up time:\t", x$time_full, "\n")
    cat("Patient arrival is modeled as a Poisson process with rate:", x$poisson_rate, "\n") 
  }
  else{
    cat("Toxicity is not a time-to-event but binary \n")
  }
}



logistic_next = function(ndose_a1, ndose_a2, target, target_min, target_max, prior_tox_a1, prior_tox_a2, final, 
                         pat_incl, dose_adm1, dose_adm2, tite=FALSE, time_full=0, toxicity=0, time_tox=0, time_follow=0,
                         c_e=0.85, c_d=0.45, c_stop=0.95, n_min){  
                         
  if(tite == TRUE) {
    toxicity = as.numeric(time_tox < time_follow)
  }  
  cdose1 = dose_adm1[pat_incl]
  cdose2 = dose_adm2[pat_incl]               
  in_startup = ifelse( sum(toxicity)>1 || (cdose1==ndose_a1 && cdose2==ndose_a2), FALSE, TRUE)
 
  n_prior_tox_a1 = length(prior_tox_a1)
  if(n_prior_tox_a1 != ndose_a1){
    stop("The entered vector of initial guessed toxicity probabities for agent 1 is of wrong length.")
  }
  n_prior_tox_a2 = length(prior_tox_a2)
  if(n_prior_tox_a2 != ndose_a2){
    stop("The entered vector of initial guessed toxicity probabities for agent 2 is of wrong length.")
  }
  
  n_toxicity = length(toxicity)
  n_time_follow = length(time_follow)
  n_time_tox = length(time_tox)
  n_dose_adm1 = length(dose_adm1)
  n_dose_adm2 = length(dose_adm2)
  if(tite==FALSE && n_toxicity != pat_incl){
    stop("The entered vector of observed toxicities is of wrong length.")
  }
  if(tite==TRUE && n_time_follow != pat_incl){
    stop("The entered vector for patients' follow-up time is of wrong length.")
  }
  if(tite==TRUE && n_time_tox != pat_incl){
    stop("The entered vector for patients' time-to-toxicity is of wrong length.")
  }
  if(n_dose_adm1 != pat_incl){
    stop("The entered vector for patients' dose of agent 1 is of wrong length.")
  }
  if(n_dose_adm2 != pat_incl){
    stop("The entered vector for patients' dose of agent 2 is of wrong length.")
  }
  
  tite = as.logical(tite)
  ndose_a1 = as.integer(ndose_a1)[1]
  ndose_a2 = as.integer(ndose_a2)[1]
  time_full = as.double(time_full)[1]
  target = as.double(target)[1] 
  target_max = as.double(target_max)[1]
  target_min = as.double(target_min)[1]
  prior_tox_a1 = as.double(prior_tox_a1)
  prior_tox_a2 = as.double(prior_tox_a2)
  final = as.logical(final)
  c_e = as.double(c_e)[1]
  c_d = as.double(c_d)[1]
  c_stop = as.double(c_stop)[1]
  n_min = as.integer(n_min)[1]
  pat_incl = as.integer(pat_incl)[1]
  cdose1 = as.integer(cdose1-1)
  cdose2 = as.integer(cdose2-1)
  in_startup = as.logical(in_startup)[1]
  dose_adm1 = as.integer(dose_adm1-1)
  dose_adm2 = as.integer(dose_adm2-1)
  time_tox = as.double(time_tox)
  time_follow = as.double(time_follow)
  toxicity = as.logical(toxicity)
  
  for(i in 1:ndose_a1){
    if(prior_tox_a1[i] < 0 || prior_tox_a1[i] > 1){
      stop("At least one of the initial guessed toxicity for agent 1 is not comprised between 0 and 1.")
    }
  }
  for(i in 1:ndose_a2){
    if(prior_tox_a2[i] < 0 || prior_tox_a2[i] > 1){
      stop("At least one of the initial guessed toxicity for agent 2 is not comprised between 0 and 1.")
    }
  }
  if(target < 0 || target > 1){stop("The toxicity target is not comprised between 0 and 1.")}
  if(target_max < 0 || target_max > 1){stop("The maximum of the targeted toxicity interval is not comprised between 0 and 1.")}
  if(target_min < 0 || target_min > 1){stop("The minimum of the targeted toxicity interval is not comprised between 0 and 1.")}
  
                
  inconc = as.logical(numeric(1))
  pi = as.double(numeric(ndose_a1*ndose_a2))
  ptox_inf = as.double(numeric(ndose_a1*ndose_a2))
  ptox_inf_targ = as.double(numeric(ndose_a1*ndose_a2))
  ptox_targ = as.double(numeric(ndose_a1*ndose_a2))
  ptox_sup_targ = as.double(numeric(ndose_a1*ndose_a2))
  
  logistic = .C(C_logistic_next, tite=tite, ndose_a1=ndose_a1, ndose_a2=ndose_a2, time_full=time_full, target=target, 
                target_max=target_max, target_min=target_min, prior_tox_a1=prior_tox_a1, prior_tox_a2=prior_tox_a2,
                final=final, c_e=c_e, c_d=c_d, c_stop=c_stop, n_min=n_min, pat_incl=pat_incl, cdose1=cdose1,
                cdose2=cdose2, in_startup=in_startup, dose_adm1=dose_adm1, dose_adm2=dose_adm2, time_tox=time_tox,
                time_follow=time_follow, toxicity=toxicity, inconc=inconc, pi=pi, ptox_inf=ptox_inf, 
                ptox_inf_targ=ptox_inf_targ, ptox_targ=ptox_targ, ptox_sup_targ=ptox_sup_targ,
                NAOK=TRUE)

  res = list(call = match.call(),
             tite=tite, 
             ndose_a1=ndose_a1, 
             ndose_a2=ndose_a2, 
             time_full=time_full, 
             target=target, 
             target_max=target_max, 
             target_min=target_min, 
             prior_tox_a1=prior_tox_a1, 
             prior_tox_a2=prior_tox_a2,
             final=final, 
             c_e=c_e, 
             c_d=c_d, 
             c_stop=c_stop, 
             n_min=n_min, 
             pat_incl=pat_incl, 
             cdose1=cdose1+1,
             cdose2=cdose2+1, 
             in_startup=in_startup, 
             dose_adm1=dose_adm1+1, 
             dose_adm2=dose_adm2+1, 
             time_tox=time_tox,
             time_follow=time_follow, 
             toxicity=toxicity, 
             inconc=logistic$inconc, 
             pi=matrix(logistic$pi, nrow=ndose_a1),
             ptox_inf=matrix(logistic$ptox_inf, nrow=ndose_a1), 
             ptox_inf_targ=matrix(logistic$ptox_inf_targ, nrow=ndose_a1), 
             ptox_targ=matrix(logistic$ptox_targ, nrow=ndose_a1),
             ptox_sup_targ=matrix(logistic$ptox_sup_targ, nrow=ndose_a1))
             
  class(res) = "logistic_next"

  return(res)
}




print.logistic_next = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call) 
  n_pat_comb = matrix(0, nrow=x$ndose_a1, ncol=x$ndose_a2)
  n_tox_comb = matrix(0, nrow=x$ndose_a1, ncol=x$ndose_a2)
  for(i in 1:x$pat_incl){
    n_pat_comb[x$dose_adm1[i],x$dose_adm2[i]] = n_pat_comb[x$dose_adm1[i],x$dose_adm2[i]]+1
    n_tox_comb[x$dose_adm1[i],x$dose_adm2[i]] = n_tox_comb[x$dose_adm1[i],x$dose_adm2[i]]+x$toxicity[i]
  }
  n_pat_comb = t(n_pat_comb)[x$ndose_a2:1,]
  n_tox_comb = t(n_tox_comb)[x$ndose_a2:1,]
  pi = t(x$pi)[x$ndose_a2:1,]
  ptox_inf = t(x$ptox_inf)[x$ndose_a2:1,]
  ptox_inf_targ = t(x$ptox_inf_targ)[x$ndose_a2:1,]
  ptox_targ = t(x$ptox_targ)[x$ndose_a2:1,]
  ptox_sup_targ = t(x$ptox_sup_targ)[x$ndose_a2:1,]   
  dimnames(n_pat_comb) = list("Agent 2"=x$ndose_a2:1, "Agent 1"=1:x$ndose_a1)
  dimnames(n_tox_comb) = list("Agent 2"=x$ndose_a2:1, "Agent 1"=1:x$ndose_a1)
  dimnames(pi) = list("Agent 2"=x$ndose_a2:1, "Agent 1"=1:x$ndose_a1)
  dimnames(ptox_inf) = list("Agent 2"=x$ndose_a2:1, "Agent 1"=1:x$ndose_a1)
  dimnames(ptox_inf_targ) = list("Agent 2"=x$ndose_a2:1, "Agent 1"=1:x$ndose_a1)
  dimnames(ptox_targ) = list("Agent 2"=x$ndose_a2:1, "Agent 1"=1:x$ndose_a1)
  dimnames(ptox_sup_targ) = list("Agent 2"=x$ndose_a2:1, "Agent 1"=1:x$ndose_a1)

  print_rnd= function (hd, x) {cat(hd, "\n"); print(round(x, digits = dgt)); cat("\n")}
  print_rnd("Number of patients:" , n_pat_comb)
  print_rnd("Number of toxicities:", n_tox_comb)
  print_rnd("Toxicity prob:", pi)
  print_rnd("P(toxicity prob < target):", ptox_inf)
  print_rnd("Prob underdosing:", ptox_inf_targ)
  print_rnd("Prob targeted interval:", ptox_targ)
  print_rnd("Prob overdosing:", ptox_sup_targ)
  
  cat("Start-up phase ended:\t", ifelse(x$in_startup==0, "YES", "NO"),"\n")
  if(!x$inconc){
    if(x$final){
      cat(paste("RECOMMENDED COMBINATION at the end of the trial:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
    }
    else{
      cat(paste("NEXT RECOMMENDED COMBINATION:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
    }
  }
  else{
    cat(paste("THE DOSE-FINDING PROCESS SHOULD BE STOPPED WITHOUT COMBINATION RECOMMENDATION\n",sep=""), sep="")
  }
  cat("\n", "\n")
  
  cat("Number of patients included:\t", round(sum(x$pat_incl),1), "\n")
  cat("Toxicity target:\t", x$target, "\n")
  cat("Targeted toxicity interval:\t [", x$target_min, ",", x$target_max, "]\n")
  cat("Prior toxicity probabilities for agent 1:\n")
  print(round(x$prior_tox_a1, digits = dgt))
  cat("Prior toxicity probabilities for agent 2:\n")
  print(round(x$prior_tox_a2, digits = dgt))
  cat("The minimum number of patients to stop the trial is:\t", x$n_min, "\n")
  cat("Escalation threshold:\t", x$c_e, "\n")
  cat("Deescalation threshold:\t", 1-x$c_d, "\n")
  cat("Stopping threshold:\t", x$c_stop, "\n")
  if (x$tite) {
    cat("Toxicity is a time-to-event \n")
    cat("Full follow-up time:\t", x$time_full, "\n")
  }
  else{
    cat("Toxicity is not a time-to-event but binary \n")
  }
}
