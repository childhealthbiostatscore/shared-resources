#########Functions for Jenna's project:
fx_model_logged<-function(y,x,data){
  # y<-"sd_glucose_pre0_L"
  # x<-"bmi_per_cat"
  # data<-dat
  f<-formula(paste0(y,"~",paste0(x,collapse = '+')))
  model<-lm(f,data)
  
  ret<-data.frame(Independent_variable=paste0(unlist(lapply(data[x],label_fx)),"+"),
                  Estimate=paste0(round(exp(coef(model)[-1]),2)," (",
                                  round(exp(coef(model)-1.96*coef(summary(model))[, 2]),2)[-1],", ", #######need to fix this for logged version before using!!
                                  round(exp(coef(model)+1.96*coef(summary(model))[, 2]),2)[-1],")"),
                  Pvalue=round(summary(model)$coef[,4],4)[-1],
                  row.names=NULL)
  
  ret$Pvalue[ret$Pvalue==0]<-'<0.0001'
  colnames(ret)<-c("Predictor","Estimate (95% CI)","P Value")
  return(ret)
}


fx_model_surv<-function(model,covs){
  
  #model<-multi_infection_pre
  #covs<-cov_pre
  ret<-data.frame(Independent_variable=colnames(covs),
                  Estimate=paste0(round(exp(model$coef),2)," (",
                                  round(exp(model$coef-1.96*summary(model)$coef[, 3]),2),", ", #######need to fix this for logged version before using!!
                                  round(exp(model$coef+1.96*summary(model)$coef[, 3]),2),")"),
                  Pvalue=round(summary(model)$coef[,5],4),
                  row.names=NULL)
  
  ret$Pvalue[ret$Pvalue==0]<-'<0.0001'
  colnames(ret)<-c("Predictor","Hazard Ratio (95% CI)","P Value")
  return(ret)
}

fx_model_nonexp<-function(model,covs){
 # covs<-sbi_0_30_covs
#  model<-sbi_0_30_model
  ret<-data.frame(Independent_variable=colnames(covs),
                  Estimate=round(coef(model)[-1],2),
                  SE=round(summary(model)$coefficients[-1,2],2),
                  Pvalue=round(summary(model)$coef[,4],4)[-1],
                  row.names=NULL)
  
  ret$Pvalue[ret$Pvalue==0]<-'<0.0001'
  colnames(ret)<-c("Predictor","Raw Estimate","SE", "P Value")
  return(ret)
  
}

fx_model_nonexp_surv<-function(model,covs){
  # covs<-cov_0_30
  #  model<-multi_gvhd_0_30
  ret<-data.frame(Independent_variable=colnames(covs),
                  Estimate=round(model$coef,2),
                  SE=round(summary(model)$coef[,3],2),
                  Pvalue=round(summary(model)$coef[,5],4),
                  row.names=NULL)
  
  ret$Pvalue[ret$Pvalue==0]<-'<0.0001'
  colnames(ret)<-c("Predictor","Raw Estimate","SE", "P Value")
  return(ret)
  
}