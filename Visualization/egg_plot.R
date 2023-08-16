egg_plot = function(xvar,yvar,group,df,log_transform = T){
  # Format data (most likely need log transformation)
  df = data.frame(df)
  df[,xvar] = as.numeric(df[,xvar])
  df[,yvar] = as.numeric(df[,yvar])
  # Baseline model for line
  f = as.formula(paste0(yvar,"~",xvar))
  if(log_transform){
    f = update(f,log(.)~log(.))
  }
  base_model <- lm(f,data = df)
  # Ellipses
  xy = subset(df,select = c(xvar,yvar,group))
  xy = data.frame(xy[complete.cases(xy),])
  xy[,group] = factor(xy[,group])
  e = dataEllipse(xy[,xvar],xy[,yvar],xy[,group],levels = 0.95)
  # Make plot df
  plot_df = data.frame(seq(min(df[,"xvar"],na.rm = T),max(df[,"xvar"],na.rm = T),length.out = 1000))
  colnames(plot_df) = xvar
  plot_df$base_pred = predict(base_model,plot_df)
  if(log_transform){
    plot_df = exp(plot_df)
  }
  # 
  base_plot = ggplot(plot_df,aes(x=m_i,y=base_pred))+geom_line()
  
  xy_cov    <- var(xy, na.rm = T)
  n_points  <- sum(complete.cases(xy))
  
  # Confidence ellipse for baseline point
  confellipse.base <- ellipse::ellipse(xy/n_points, centre=c(m,e))
  
  # Predicted Baseline curve data
  baseline <- data.frame(x=exp(base_model$model$x),y=exp(predict(base_model))) 
  baseline <- baseline[order(baseline$x),]
}