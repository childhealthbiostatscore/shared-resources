# This is a function for performing cross validation (CV) to select an optimal model 
# using the ElasticNet. By default uses leave one out (LOO) CV, but k-fold CV
# can also be used by setting cv_method = "cv" and folds = k. See the trainControl 
# function in caret for additional details. There are two options for the output:
# out = "min.error" produces the model with the lowest CV error, and 
# out = "1se.error" produces all "acceptable" models (CV error within 
# 1 standard error of the minimum). By default RMSE is the metric used for 
# continuous variables, and accuracy for binary, but different metrics can be 
# specified with the metric argument. Variables must be formatted correctly
# (e.g. sex as a factor variable) for this to work properly! The only exception 
# is that outcomes will automatically be converted to factor variables if binom = T.
easy_glinternet = function(data,outcome,predictors,binom = F,metric = NULL,
                           cv_method = "LOOCV",folds = NULL,na = "na.omit",
                           out = "min.error"){
  require(caret)
  require(glmnet)
  df = data
  # Make outcome binary if necessary
  if(binom){
    df[,outcome] = as.factor(df[,outcome])
    if(is.null(metric)){
      metric = "Accuracy"
    }
  }
  if(is.null(metric)){
    metric = "RMSE"
  }
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
  # Minimum error model
  if (out == "min.error") {
    params = t$bestTune
  } else if (out == "1se.error") {
    res = t$results
    if(binom & metric == "Accuracy"){
      res[,metric] = 1 - res[,metric]
    }
    min_err = min(res[,metric],na.rm = T)
    se_err = sd(res[,metric],na.rm = T)/sqrt(sum(!is.na(res[,metric])))
    good_mods = which(res[,metric] <= (min_err + se_err))
    params = res[good_mods,]
  }
  # Refit models to get selected parameters (the coef() function output for caret is confusing)
  # Matrices
  X = model.matrix(update(f,.~.-1),df)
  Y = df[rownames(X),outcome]
  mods = apply(params,1,function(r){
    a = as.numeric(r["alpha"])
    l = as.numeric(r["lambda"])
    if(binom){
      mod = glmnet(y = Y,x = X,alpha = a,lambda = l,family = "binomial")
    } else {
      mod = glmnet(y = Y,x = X,alpha = a,lambda = l)
    }
    selected = as.matrix(coef(mod))
    selected = rownames(selected)[selected[,1] != 0]
    selected = selected[selected != "(Intercept)"]
    return(selected)
  })
  names(mods) = NULL
  return(mods)
}
