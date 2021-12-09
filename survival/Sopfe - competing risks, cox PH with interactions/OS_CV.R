######Overall Survival:

options(scipen=999)

dat$death_time_table<-NA
dat$death_time_table[dat$pt_alive=="Died"]<-dat$time_death[dat$pt_alive=="Died"]
label(dat$death_time_table)<-"Days to Death"
dat$death_time_table_cat<-NA
dat$death_time_table_cat[dat$death_time_table<=30]<-"<30 days"
dat$death_time_table_cat[dat$death_time_table>30]<-"30+ days"


os_table<-final_table(dat,c('pt_alive','death_time_table','death_time_table_cat'), 
                      dat$tx_type,margin=2,single=F,2,col.names=T, summary.stat='median')
os_table

hist(dat$death_time_table,main="Time to Death",xlab="Time to Death",breaks = 10,xaxt="n"
)
axis(1,at=c(0,30,100,200,400,2000),labels=c(0,30,100,200,400,2000))


################Time-to-death Models###############
dat$death_yn<-0
dat$death_yn[dat$pt_alive=="Died"]<-1
dat$tx_type<-factor(dat$tx_type,levels=c('Autologous','Allogenic'))
#dat$tx_type<-factor(dat$tx_type,levels=c('Allogenic','Autologous'))

dat.pre<-subset(dat,!is.na(dat$cv_glucose_pre0_log2))
dat.pre<-subset(dat.pre,dat.pre$cv_glucose_pre0_log2!="-Inf")
cov_pre <- model.matrix(~ cv_glucose_pre0_log2 + post_tx_steroids+gvhd_new+
                          tx_type*cv_glucose_pre0_log2
                        ,
                        data = dat.pre)[, -1]


os_pre<-coxph(Surv(time_death,death_yn)~cov_pre,data=dat.pre)

os_pre.plot1<-coxph(Surv(time_death,death_yn)~cov_pre,data=dat.pre)

exp(os_pre$coefficients)
exp(confint(os_pre))
os_pre<-fx_model_surv(os_pre,cov_pre)

#os_pre<-fx_model_nonexp_surv(os_pre,cov_pre)

dat.0_30<-subset(dat,!is.na(dat$cv_glucose_0_30_log2))
#remove those with events prior to 30 days
dat.0_30<-subset(dat.0_30,!(dat.0_30$time_death<=30 & dat.0_30$death_yn==1))
#dat.0_30$tx_type<-factor(dat.0_30$tx_type,levels=c('Autologous','Allogenic'))
#dat.0_30$tx_type<-factor(dat.0_30$tx_type,levels=c('Allogenic','Autologous'))
cov_0_30 <- model.matrix(~ cv_glucose_0_30_log2+ post_tx_steroids + gvhd_new+
                           tx_type+cv_glucose_0_30_log2*tx_type
                           ,
                         data = dat.0_30)[, -1]

os_0_30<-coxph(Surv(time_death,death_yn)~cov_0_30,data=dat.0_30)
os_0_30.plot1<-coxph(Surv(time_death,death_yn)~cov_0_30,data=dat.0_30)

exp(os_0_30$coefficients)
exp(confint(os_0_30))
os_0_30<-fx_model_surv(os_0_30,cov_0_30)

#os_0_30<-fx_model_nonexp_surv(os_0_30,cov_0_30)
# 
# 
# dat.0_100 <-subset(dat,!is.na(dat$cv_glucose_0_100_log2))
# dat.0_100$tx_type<-factor(dat.0_100$tx_type,levels=c('Autologous','Allogenic'))
# cov_0_100 <- model.matrix(~ cv_glucose_0_100_log2+ post_tx_steroids +
#                             gvhd_new+tx_type+cv_glucose_0_100_log2*tx_type
#                             ,
#                           data = dat.0_100)[, -1]
# os_0_100<-coxph(Surv(time_death,death_yn)~cov_0_100,data=dat.0_100)
# os_0_100<-fx_model_nonexp_surv(os_0_100,cov_0_100)
# 
# 

