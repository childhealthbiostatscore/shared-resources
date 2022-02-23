death_mod<-function(x){
  #x<-dat$vent_2
  #f<-formula(paste0(dat$death_yn,"~",paste0(x,collapse = '+')))
  model<-glm(dat$death_yn~x,dat,family='binomial')
  if (!is.factor(x)){
    lab<-label(x)
  }
  if (is.factor(x)){
    lab<-levels(x)[-1]
  }
  ret<-data.frame(Independent_variable=lab,
                  OR=paste0(round(exp(coef(model)[-1]),2)," (",
                            round(exp(confint(model)[,1]),2)[-1],", ",
                            round(exp(confint(model)[,2]),2)[-1],")"),
                  Pvalue=round(summary(model)$coef[,4],4)[-1],
                  row.names=NULL)
  
  ret$Pvalue[ret$Pvalue==0]<-'<0.0001'
  colnames(ret)<-c("Parameter","OR (95% CI)","P Value")
  return(ret)
}
death_mod_surv<-function(x){
  #x<-dat$HGG
  # x<-c("CSF_LYMPHS..","CSF_Malignant","vent_2","ATRT","Germinoma","HGG","LGG","OET")
  # f<-paste0(x,collapse = '+')
  model<-coxph(Surv(time_to_death, death_yn_num) ~ x, dat)
  sum<-summary(model)
  if (!is.factor(x)){
    lab<-label(x)
  }
  if (is.factor(x)){
    lab<-levels(x)[-1]
  }
  ret<-data.frame(Independent_variable=lab,
                  HR=paste0(round(sum$coefficients[,2],2)," (",
                            round(exp(sum$coefficients[,1]-1.96*sum$coefficients[,3]),2),", ",
                            round(exp(sum$coefficients[,1]+1.96*sum$coefficients[,3]),2),")"),
                  Pvalue=round(sum$coefficients[,5],4),
                  row.names=NULL)
  
  ret$Pvalue[ret$Pvalue==0]<-'<0.0001'
  colnames(ret)<-c("Parameter","HR (95% CI)","P Value")
  return(ret)
}

death_mod_surv_multi<-function(model){
# #model<-coxph(Surv(time_to_death, death_yn_num) ~ 
#   dat.complete$Lymphs.absolute
#   +dat.complete$HGG
#   +dat.complete$ATRT
#   +dat.complete$location_2
#   +dat.complete$OET
#   +dat.complete$CSF_Malignant,
#   dat.complete)
  sum<-summary(model)
  
  Pvalue=paste0("p=",round(sum$coefficients[,5],4))
  Pvalue[Pvalue=="p=0"]<-'p<0.0001'

  ret<-data.frame(Parameter=names(model$coefficients),
                  HR=paste0(round(sum$coefficients[,2],2)," (",
                            round(exp(sum$coefficients[,1]-1.96*sum$coefficients[,3]),2),", ",
                            round(exp(sum$coefficients[,1]+1.96*sum$coefficients[,3]),2),"); ",Pvalue),
                  row.names=NULL)
  

  return(ret)
}


zinb_mod<-function(x,y){
  #x<-dat$vent_2
  #f<-formula(paste0(dat$death_yn,"~",paste0(x,collapse = '+')))
  #model<-glm(dat$death_yn~x,dat,family='binomial')
  # y<-dat$CSF_LYMPHS..
  # x<-dat$CSF_Malignant
  model <- zeroinfl(y ~ x,
                 dist = 'negbin',
                 data = dat)
  if (!is.factor(x)){
    lab<-label(x)
  }
  if (is.factor(x)){
    lab<-levels(x)[-1]
  }
  
  sum<-summary(model)
  sum_count<-sum$coefficients$count
  sum_count<-sum_count[-c(1,nrow(sum_count)),]
  sum_zero<-sum$coefficients$zero
  sum_zero<-sum_zero[-1,]
  ret<-data.frame(Independent_variable=lab,
                  Estimate_count=paste0(round(exp(sum_count[1]),2)," (",
                            round(exp(sum_count[1]-1.96*sum_count[2]),2),", ",
                            round(exp(sum_count[1]+1.96*sum_count[2]),2),"); p=",round(sum_count[4],4)),
                  Estimate_zero=paste0(round(exp(sum_zero[1]),2)," (",
                                        round(exp(sum_zero[1]-1.96*sum_zero[2]),2),", ",
                                        round(exp(sum_zero[1]+1.96*sum_zero[2]),2),"); p=",round(sum_zero[4],4)),
                  row.names=NULL)
  return(ret)
}


zinb_mod_multi<-function(model){
  #x<-dat$vent_2
  #f<-formula(paste0(dat$death_yn,"~",paste0(x,collapse = '+')))
  #model<-lymphs_perc_multi
  # y<-dat$CSF_LYMPHS..
  # x<-dat$CSF_Malignant

  sum<-summary(model)
  sum_count<-sum$coefficients$count
  sum_count<-sum_count[-c(1,nrow(sum_count)),]
  sum_zero<-sum$coefficients$zero
  sum_zero<-sum_zero[-1,]
  ret<-data.frame(Independent_variable=names(model$coefficients$count)[-1],
                  Estimate_count=paste0(round(exp(sum_count[,1]),2)," (",
                                        round(exp(sum_count[,1]-1.96*sum_count[,2]),2),", ",
                                        round(exp(sum_count[,1]+1.96*sum_count[,2]),2),"); p=",round(sum_count[,4],4)),
                  Estimate_zero=paste0(round(exp(sum_zero[,1]),2)," (",
                                       round(exp(sum_zero[,1]-1.96*sum_zero[,2]),2),", ",
                                       round(exp(sum_zero[,1]+1.96*sum_zero[,2]),2),"); p=",round(sum_zero[,4],4)),
                  row.names=NULL)
  return(ret)
}

log_mod<-function(x,y){
  #x<-dat$CSF_Malignant
  #y<-dat$CSF_NUCLEATED_CELLS_yn
  #f<-formula(paste0(dat$death_yn,"~",paste0(x,collapse = '+')))
  model<-glm(y~x,dat,family='binomial')
  if (!is.factor(x)){
    lab<-label(x)
  }
  if (is.factor(x)){
    lab<-levels(x)[-1]
  }
  ret<-data.frame(Independent_variable=lab,
                  OR=paste0(round(exp(coef(model)[-1]),2)," (",
                            round(exp(confint(model)[,1]),2)[-1],", ",
                            round(exp(confint(model)[,2]),2)[-1],")"),
                  Pvalue=round(summary(model)$coef[,4],4)[-1],
                  row.names=NULL)
  
  ret$Pvalue[ret$Pvalue==0]<-'<0.0001'
  colnames(ret)<-c("Parameter","OR (95% CI)","P Value")
  return(ret)
}