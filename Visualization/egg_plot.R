egg_plot = function(xvar,yvar,group,df,log_transform = T,
                    xlabel,ylabel){
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
  # Large predicted value df for plotting a smooth looking line
  line_df = data.frame(seq(min(df[,xvar],na.rm = T),max(df[,xvar],na.rm = T),length.out=1000))
  colnames(line_df) = xvar
  line_df$pred = predict(base_model,line_df)
  # Plot
  p = ggplot(df,aes(x=.data[[xvar]],y=.data[[yvar]],group = .data[[group]],fill=.data[[group]],color=.data[[group]]))+
    stat_ellipse(na.rm = T,geom = "polygon",alpha=0.2) + 
    geom_point(na.rm = T) +
    geom_line(data = line_df,aes(x=m_i,y=exp(pred)),inherit.aes = F,color="black") + 
    xlab(xlabel)+ylab(ylabel)+
    theme_bw()+ 
    theme(legend.title = element_blank())
  return(p)
}
