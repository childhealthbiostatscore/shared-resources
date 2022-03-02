# For testing
load("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/BDC/Projects/Janet Snell-Bergeon/AHA collaborative grant/aha_master_data_no_snps.Rdata")

# This is a function for performing cross validation (CV) to select an optimal model 
# using the ElasticNet. By default uses leave one out (LOO) CV, but k-fold CV
# can also be used by setting cv_method = "cv" and folds = k. See the trainControl 
# function in caret for additional details. There are two options for the output:
# out = "min.error" produces the model with the lowest CV error, and 
# out = "1se.error" produces all "acceptable" models (CV error within 
# 1 standard error of the minimum).
easy_glinternet = function(df,outcome,predictors,metric = "RMSE",
                           cv_method = "LOOCV",folds = NULL,na = "na.omit",
                           out = "min.error"){
  require(caret)
  # CV
  cv = trainControl(method = cv_method,number = folds)
  # Model formula
  # Fix names if necessary
  colnames(df) = make.names(colnames(df),unique = T,allow_ = F)
  preds = make.names(predictors,unique = T,allow_ = F)
  f = as.formula(paste0(outcome,"~",paste0(preds,collapse = "+")))
  # Train model
  t = train(f,data = df,method = "glmnet",na.action = na)
  # Get parameters for refitting
  res = t$results
  min_err = min(res[,metric],na.rm = T)
  se_err = sd(res[,metric],na.rm = T)/sqrt(sum(!is.na(res[,metric])))
  good_mods = which(res[,metric] <= (min_err + se_err))
  params = res[good_mods,]
}
