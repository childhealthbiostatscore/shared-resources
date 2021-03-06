# This is a function for performing cross validation (CV) to select an optimal model 
# using the ElasticNet. By default uses leave one out (LOO) CV, but k-fold CV
# can also be used by setting cv_method = "kfold" and folds = k. See the trainControl 
# function in caret for additional details. There are two options for the output:
# out = "min.error" produces the model with the lowest CV error, and 
# out = "1se.error" produces all "acceptable" models (CV error within 
# 1 standard error of the minimum). 
easy_lasso = function(data,outcome,predictors,
                           n_lambdas = 100,max_coef = NULL,
                           model_type = "gaussian",time = NULL,
                           cv_method = "loo",folds = NULL,out = "min.error",
                           cores = 4,seed = 1017){
  require(glmnet)
  df = data
  # Random seed
  set.seed(seed)
  # Fix names if necessary
  colnames(df) = make.names(colnames(df),unique = T,allow_ = F)
  preds = make.names(predictors,unique = T,allow_ = F)
  outcome = make.names(outcome,unique = T,allow_ = F)
  # Predictor matrix
  X = data.frame(df[,preds])
  # Outcome matrix depending on model type
  if(model_type == "cox"){
    # Outcome matrix
    Y = cbind(time = df[,time], status = df[,outcome])
    # Complete cases
    idx = intersect(which(complete.cases(Y)),which(complete.cases(X)))
    X = data.matrix(X[idx,])
    Y = data.matrix(Y[idx,])
    # Remove variables without any variance
    near_zero = caret::nearZeroVar(X)
    if(length(near_zero)>0){
      X = X[,-near_zero]
    }
  } else if (model_type == "binomial" | model_type == "gaussian"){
    # Outcome matrix
    Y = df[,outcome]
    # Complete cases
    idx = intersect(which(complete.cases(Y)),which(complete.cases(X)))
    X = data.matrix(X[idx,])
    Y = as.numeric(Y[idx])
    # Remove variables without any variance
    near_zero = caret::nearZeroVar(X)
    if(length(near_zero)>0){
      X = X[,-near_zero]
    }
  }
  # CV parameters
  if(cv_method == "loo"){
    folds = nrow(X)
  } else if (cv_method != "kfold"){
    stop("Please select either LOO or k-fold CV. If you are selecting k-fold CV, please specify the number of folds.")
  }
  # Parallel - recommended
  if(!is.null(cores)){
    require(doParallel)
    p = TRUE
    registerDoParallel(cores)
  } else {p = FALSE}
  # glmnet
  lasso_cv = cv.glmnet(X,Y,family = model_type,nfolds = folds,parallel = p)
  # Refit models to get selected parameters
  if(out == "min.error"){
    l = lasso_cv$lambda.min
  } else if (out == "1se.error"){
    l = lasso_cv$lambda.1se
  }
  # Return coefficients
  mod = glmnet(X,Y,family = model_type,parallel = p,lambda = l)
  c = coef(mod)
  c = rownames(c)[which(c != 0 & rownames(c) != "(Intercept)")]
  return(c)
}
