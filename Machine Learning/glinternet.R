# This is a function for formatting input to and output from the glinternet 
# package. The input dataframe must be correctly formatted with categorical 
# variables stored as factors.
easy_glinternet = function(df,outcome = "Y",n_lambda = 100,
                           fam = "gaussian",cores = 1,
                           id_vars = "record_id",var_names = NULL,
                           cv_type = "k_fold",folds = 5,random_seed = 3654){
  require(glinternet)
  require(caret)
  require(parallel)
  require(pROC)
  # Remove NaNs
  df = data.frame(lapply(df,function(c){
    c[is.nan(c)] = NA
    return(c)
  }))
  # Separate outcome from predictors
  if(is.null(var_names)){
    df = df[complete.cases(df),]
  } else {
    df = df[,c(id_vars,outcome,var_names)]
    df = df[complete.cases(df),]
  }
  # Outcome and predictors
  Y = df[,outcome]
  X = df[,!names(df) %in% c(id_vars,outcome)]
  # Variables names and number of levels
  vars = names(X)
  types = sapply(X,class)
  num_levels = sapply(X,function(c){length(levels(c))})
  num_levels[num_levels == 0] = 1
  cat_vars = vars[types == "factor"]
  cont_vars = vars[types == "numeric" | types == "integer"]
  # Numeric outcome and predictors
  if(fam == "binomial"){
    Y = as.numeric(Y) - 1
  } else {
    Y = data.matrix(Y)
  }
  X = data.matrix(X)
  X[,which(types == "factor")] = X[,which(types == "factor")] - 1
  # Get lambdas to test from full dataset
  glin_full = glinternet(X,Y,num_levels,nLambda = n_lambda,family = fam,numCores = cores)
  lambdas = glin_full$lambda
  # Caret for cross validation - cv_type can be k_fold or loo
  if(cv_type == "loo"){
    train_sets = createFolds(df[,outcome],k = nrow(df),returnTrain = T)
  } else if (cv_type == "k_fold"){
    set.seed(random_seed)
    train_sets = createFolds(df[,outcome],k = folds,returnTrain = T)
  } 
  # For each lambda, perform LOO CV
  lambda_error = sapply(lambdas, function(l){
    cl = makeCluster(cores,type="PSOCK")
    clusterExport(cl,list('train_sets','l','X','Y','glinternet','fam',
                          'num_levels'),envir = environment())
    predictions = t(parSapply(cl,train_sets, function(t){
      # Split
      train_out = Y[t]
      train_pred = X[t,]
      test_out = Y[-t]
      # Some bug in glinternet prevents predict from working on a single observation, so duplicate that row
      test_pred = rbind(X[-t,],X[-t,])
      # Fit
      glin = glinternet(train_pred,train_out,num_levels,lambda = c(0.1,l),family = fam)
      # Predicted values and true value
      p = predict(glin,test_pred,type = "response",lambda = l)
      return(c(p[1],test_out))
    }))
    stopCluster(cl)
    # Get AUC or RMSE for this lambda value
    if (fam == "binomial"){
      predicted = ifelse(predictions[,1] > .5, 1, 0)
      r = sum(predictions[,2] == predicted)/length(predictions[,2]) # Accuracy
      r = 1 - r
    } else if (fam == "gaussian"){
      r = RMSE(pred = predictions[,1],obs = predictions[,2])
    }
    return(r)
  })
  # Re-fit with best lambda and pull coefficients
  lambdas = lambdas[which(lambda_error > 0)]
  lambda_error = lambda_error[which(lambda_error > 0)]
  best_l = lambdas[which.min(lambda_error)]
  fit = glinternet(X,Y,num_levels,lambda = lambdas,family = fam,numCores = cores)
  w = which.min(abs(fit$lambda - best_l))
  coefs = coef(fit)[[w]]
  # Make formula
  # Main effects
  if(!is.null(coefs$mainEffects$cont)){
    cont = cont_vars[coefs$mainEffects$cont]
  } else {cont = NULL}
  if (!is.null(coefs$mainEffects$cat)){
    cat = cat_vars[coefs$mainEffects$cat]
  } else {cat = NULL}
  # Interactions
  if(!is.null(coefs$interactions$contcont)){
    contcont = apply(coefs$interactions$contcont,1,function(r){
      paste0(cont_vars[r],collapse = ":")
    })
    contcont = paste0(contcont,collapse = "+")
  } else {contcont = NULL}
  if(!is.null(coefs$interactions$catcont)){
    catcont = apply(coefs$interactions$catcont,1,function(r){
      paste0(c(cat_vars[r[1]],cont_vars[r[2]]),collapse = ":")
    })
    catcont = paste0(catcont,collapse = "+")
  } else {catcont = NULL}
  if(!is.null(coefs$interactions$catcat)){
    catcat = apply(coefs$interactions$catcat,1,function(r){
      paste0(cat_vars[r],collapse = ":")
    })
    catcat = paste0(catcat,collapse = "+")
  } else {catcat = NULL}
  f = paste0(outcome,"~",paste0(c(cont,cat,contcont,catcont,catcat),collapse = "+"))
  return(f)
}
