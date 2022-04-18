#Useful Code

#label all columns with colname
for (i in c(1:ncol(dat))){
  label(dat[,i])<-paste0(colnames(dat[i]))
}  

#convert variables to 0/1 (combine with lapply to do lots of columns)
convert_num<-function(var){
  ## var<-dat.elig$pain_p
  dat.elig$temp<-NA
  dat.elig$temp[var=="No"]<-0
  dat.elig$temp[var=="Yes"]<-1
  var<-dat.elig$temp
  return(var)
}

#if zero, make missing. for things like "number of days in ICU, if ICU=1"
num_yes<-function(x){
  #x<-dat$meds_gastrointestinal_num
  temp<-x
  temp.2<-NA
  temp.2<-ifelse(x>0,temp,NA)
  return(temp.2)
}

#refactor according to number in each factor, useful for Table 1 sometimes
sortLvlsByN.fnc <- function(oldFactor, ascending = TRUE) {
  # Return factor with new level order
  newFactor <- factor(oldFactor, levels = table(oldFactor)  %>% sort(., decreasing = ascending)  %>% names())
  return(newFactor)
}

#### basic functions for summarizing variables into table format
mean_table<-function(x){
  # x<-dat$Distance.from.Hospital..miles.
  mean.1<-mean(x,na.rm=T)
  sd.1<-sd(x,na.rm=T)
  miss.1<-table(is.na(x),useNA="always")
  miss.prop<-prop.table(miss.1)
  #table for x
  tab.x<-data.frame(name=label(x),
                    Est=paste0(round(mean.1,2)," (",round(sd.1,2),")"),
                    Missing=paste0(miss.1[2]," (",miss.prop[2],"%)"))
  return(tab.x)
}

median_table<-function(x){
  # x<-dat[c("Distance.from.Hospital..miles.")]
  q.1<-quantile(x,na.rm=T)
  miss.1<-table(is.na(x),useNA="always")
  miss.prop<-prop.table(miss.1)
  #table for x
  tab.x<-data.frame(name=label(x),
                    Est=paste0(round(q.1[3],2)," (",round(q.1[2],2),", ",round(q.1[4],2),")"),
                    Missing=paste0(miss.1[2]," (",miss.prop[2],"%)"))
  return(tab.x)
}

prop_table<- function(x) {
  #x<-dat$Race..choice.Native.Hawaiian.or.Other.Pacific.Islander.
  y<-table(x,useNA="always")
  y.prop<-prop.table(y)
  
  estimate<-c("",paste(y[],paste0("(",round(y.prop*100,0),"%",")")))
  miss<-estimate[length(estimate)]
  estimate<-estimate[-length(estimate)]
  
  #table
  tab.x<-data.frame(name=c(label(x),paste(levels(as.factor(x)),sep=",")),
                    Est=estimate,
                    Missing=miss
  )
  
  return(tab.x)
}

mean_table_2<-function(x1,x2,name){
  
  mean.x1<-t.test(x1,na.rm=T)
  test.t1<-t.test(x1,conf.level = 0.95,na.rm=T)
  mean.x2<-t.test(x2,na.rm=T)
  test.t2<-t.test(x2,conf.level = 0.95,na.rm=T)
  test<-t.test(x1,x2,na.rm=T)
  pval<-as.character(round(test$p.value,2))
  pval[pval<0.001]<-"<0.001"
  #table for x
  tab.x<-data.frame(
    Summary_meas=name,
    No=paste0(round(test.t1$estimate[1],2)," ","(",round(sd(x1,na.rm=T),2),")"),
    Yes=paste0(round(test.t2$estimate[1],2)," ","(",round(sd(x2,na.rm=T),2),")"),
    Pval=pval)
  return(tab.x)
}