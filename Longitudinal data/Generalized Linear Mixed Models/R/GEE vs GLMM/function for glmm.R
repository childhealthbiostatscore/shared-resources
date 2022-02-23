logit_gee<- function(var1,outcome) {
  #var1<-dat$Lactate...initial.within.8.hours
  #outcome<-dat$CDE.within.48.hrs_01
  # if (is.factor(var1)){
  #   #var1_use<-rescale(as.numeric(var1),binary.inputs="center")
  #   var1_use<-var1
  #   label(var1_use)<-label(var1)
  # }
  # 
  # if (!(is.factor(var1))) {
  #   var1<-scale(var1)
  # }
  # 
  #####use this portion (copy outside of function) for multivariable model#####
  mod<-gee(outcome ~ var1,data = dat, id = id_n, family = binomial,
           corstr = "exchangeable")
  
  mod_sum<-data.frame(summary(mod)$coefficients)
  mod_sum$p<-2 * pnorm(abs(mod_sum$Robust.z), lower.tail = FALSE)
  mod_sum$Var<-row.names(mod_sum)
  
  mod_sum<-mod_sum[-1,-c(2,3,5)]
  p<-round(mod_sum$p,4)
  p[p<0.001]<-"<0.001"
  
  
  #table
  if (is.factor(var1)){
    tab.1<-data.frame(Var=paste0(label(var1),": ",levels(as.factor(var1))[-1]),
                      ref=levels(as.factor(var1))[1],
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      nobs=mod$nobs,
                      row.names=NULL)
  }
  if (!is.factor(var1)){
    tab.1<-data.frame(Var=label(var1),
                      ref=NA,
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      nobs=mod$nobs,
                      row.names=NULL)
  }
  #}
  #####copy until here#####
  
  return(tab.1)
}

cont_logged_gee<- function(var1,outcome) {
  #var1<-dat$Lactate...initial.within.8.hours
  #outcome<-dat$CDE.within.48.hrs_01
  # if (is.factor(var1)){
  #   #var1_use<-rescale(as.numeric(var1),binary.inputs="center")
  #   var1_use<-var1
  #   label(var1_use)<-label(var1)
  # }
  # 
  # if (!(is.factor(var1))) {
  #   var1<-scale(var1)
  # }
  # 
  #####use this portion (copy outside of function) for multivariable model#####
  mod<-gee(outcome ~ var1,data = dat, id = id_n, family = gaussian,
           corstr = "exchangeable")
  
  mod_sum<-data.frame(summary(mod)$coefficients)
  mod_sum$p<-2 * pnorm(abs(mod_sum$Robust.z), lower.tail = FALSE)
  mod_sum$Var<-row.names(mod_sum)
  
  mod_sum<-mod_sum[-1,-c(2,3,5)]
  p<-round(mod_sum$p,4)
  p[p<0.001]<-"<0.001"
  
  
  #table - EXPONENTIATED VALUES, DUE TO LOG TRANSFORMED OUTCOME
  if (is.factor(var1)){
    tab.1<-data.frame(Var=paste0(label(var1),": ",levels(as.factor(var1))[-1]),
                      ref=levels(as.factor(var1))[1],
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      nobs=mod$nobs,
                      row.names=NULL)
  }
  if (!is.factor(var1)){
    tab.1<-data.frame(Var=label(var1),
                      ref=NA,
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      nobs=mod$nobs,
                      row.names=NULL)
  }
  #}
  #####copy until here#####
  
  return(tab.1)
}




logit_mixed<- function(var1,outcome) {
  #var1<-dat$Lactate...initial.within.8.hours
  #outcome<-dat$CDE.within.48.hrs_01
  # if (is.factor(var1)){
  #   #var1_use<-rescale(as.numeric(var1),binary.inputs="center")
  #   var1_use<-var1
  #   label(var1_use)<-label(var1)
  # }
  
  if (!(is.factor(var1))) {
    var1<-scale(var1)
  }
  # 
  #####use this portion (copy outside of function) for multivariable model#####
  mod<-glmer(outcome ~ dat$encounter_num+var1 +(1 | Subject.ID), data = dat, family = binomial)
  mod_sum<-data.frame(summary(mod)$coefficients)
  mod_sum$Var<-row.names(mod_sum)
  mod_sum<-mod_sum[-c(1,2),]
  p<-round(mod_sum$Pr...z..,4)
  p[p<0.001]<-"<0.001"
  
  #table
  if (is.factor(var1)){
    tab.1<-data.frame(Var=paste0(label(var1),": ",levels(as.factor(var1))[-1]),
                      ref=levels(as.factor(var1))[1],
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      n_encounters=nobs(mod),
                      n_subjects=summary(mod)$ngrps[1],
                      row.names=NULL)
  }
  if (!is.factor(var1)){
    tab.1<-data.frame(Var=label(var1),
                      ref=NA,
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      n_encounters=nobs(mod),
                      n_subjects=summary(mod)$ngrps[1],
                      row.names=NULL)
  }
  #}
  #####copy until here#####
  
  return(tab.1)
}

logit_mod<- function(var1,outcome) {
  #var1<-dat$Age..years.
  #outcome<-dat$CDE.within.48.hrs
  # if (is.factor(var1)){
  #   #var1_use<-rescale(as.numeric(var1),binary.inputs="center")
  #   var1_use<-var1
  #   label(var1_use)<-label(var1)
  # }
  if (!(is.factor(var1))) {
    var1<-scale(var1)
  }
  # 
  #####use this portion (copy outside of function) for multivariable model#####
  mod<-glm(outcome ~ dat$encounter_num +var1, data = dat, family = binomial)
  mod_sum<-data.frame(summary(mod)$coefficients)
  mod_sum$Var<-row.names(mod_sum)
  # mod_sum<-subset(mod_sum,mod_sum$Var!="(Intercept)")
  mod_sum<-mod_sum[-c(1,2),]
  p<-round(mod_sum$Pr...z..,4)
  p[p<0.001]<-"<0.001"
  
  #table
  if (is.factor(var1)){
    tab.1<-data.frame(Var=paste0(label(var1),": ",levels(as.factor(var1))[-1]),
                      ref=levels(as.factor(var1))[1],
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      n_encounters=nobs(mod),
                      n_subjects=NA,
                      row.names=NULL)
  }
  if (!is.factor(var1)){
    tab.1<-data.frame(Var=label(var1),
                      ref=NA,
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      n_encounters=nobs(mod),
                      n_subjects=NA,
                      row.names=NULL)
  }
  #}
  #####copy until here#####
  
  return(tab.1)
}