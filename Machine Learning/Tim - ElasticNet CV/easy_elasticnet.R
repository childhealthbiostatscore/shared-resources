load("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/BDC/Projects/Janet Snell-Bergeon/AHA collaborative grant/aha_master_data_no_snps.Rdata")
# This is a function for performing cross validation (CV) to select an optimal model 
# using the ElasticNet. By default uses leave one out (LOO) CV, but k-fold CV
# can also be used by setting cv_method = "cv" and folds = k. See the trainControl 
# function in caret for additional details. There are two options for the output:
# out = "min.error" produces the model with the lowest CV error, and 
# out = "1se.error" produces all "acceptable" models (CV error within 
# 1 standard error of the minimum). 
easy_glinternet = function(data,outcome,predictors,
                           n_alphas = 10,n_lambdas = 100,
                           model_type = "gaussian",time = NULL,
                           cv_method = "LOOCV",folds = NULL){
  require(ensr)
  df = data
  # Fix names if necessary
  colnames(df) = make.names(colnames(df),unique = T,allow_ = F)
  preds = make.names(predictors,unique = T,allow_ = F)
  # Predictor matrix
  X = df[,preds]
  # Outcome matrix depending on model type
  if(model_type == "cox"){
    require(survival)
    # Outcome matrix
    Y = cbind(time = df[,time], status = df[,outcome])
    # Complete cases
    idx = intersect(which(complete.cases(Y)),which(complete.cases(X)))
    X = data.matrix(X[idx,])
    Y = data.matrix(Y[idx,])
  } else if (model_type == "binomial" | model_type == "gaussian"){
    # Outcome matrix
    Y = df[,outcome]
    # Complete cases
    idx = intersect(which(complete.cases(Y)),which(complete.cases(X)))
    X = data.matrix(X[idx,])
    Y = as.numeric(Y[idx])
  }
  # CV parameters
  if(cv_method == "LOOCV"){folds = nrow(X)}
  # Grid search with glmnet - super slow
  e = ensr(x = X,y = Y,alphas = seq(0, 1, length = n_alphas),nlambda = n_lambdas,
           family = model_type,nfolds = folds,grouped=FALSE)
  # Get alpha and lambdas
  res = summary(e)
  min_err = min(res$cvm,na.rm = T)
  se_err = sd(res$cvm,na.rm = T)/sqrt(sum(!is.na(res$cvm)))
  good_mods = which(res$cvm <= (min_err + se_err))
  params = data.frame(res[good_mods,])
  # Refit models to get selected parameters (the coef() function output for caret is confusing)
  mods = apply(params,1,function(r){
    a = as.numeric(r["alpha"])
    l = as.numeric(r["lambda"])
    mod = glmnet(y = Y,x = X,alpha = a,lambda = l,family = model_type)
    selected = as.matrix(coef(mod))
    selected = rownames(selected)[selected[,1] != 0]
    selected = selected[selected != "(Intercept)"]
    return(selected)
  })
  names(mods) = NULL
  return(mods)
}
