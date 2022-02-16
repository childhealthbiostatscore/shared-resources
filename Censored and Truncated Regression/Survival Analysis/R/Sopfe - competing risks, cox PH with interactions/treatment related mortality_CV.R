#######TRM analyses#########

library(cmprsk)
options(scipen=999)

#somehow sas program missed one subject -if all cause of death are missing, it should be a "2" which is unrelated
dat$death_tx[dat$study_id==162]<-2
dat$time_death_tx[dat$study_id==162]<-dat$time_death[dat$study_id==162]


dat$trm_tab<-NA
dat$trm_tab[dat$death_tx==0]<-"Alive"
dat$trm_tab[dat$death_tx==1]<-"Treatment Related Death"
dat$trm_tab[dat$death_tx==2]<-"Other Death: primary disease or unrelated"
dat$trm_tab<-as.factor(dat$trm_tab)
dat$trm_tab<-factor(dat$trm_tab,levels=c("Treatment Related Death","Other Death: primary disease or unrelated",
                                         "Alive"))
label(dat$trm_tab)<-"Treatment Related Mortality"
dat$trm_time_table<-NA
dat$trm_time_table[dat$death_tx==1]<-dat$time_death_tx[dat$death_tx==1]
label(dat$trm_time_table)<-"Days to Treatment Related Mortality"
dat$trm_time_table_cat<-NA
dat$trm_time_table_cat[dat$death_tx==1 & dat$time_death_tx<=30]<-"0-30"
dat$trm_time_table_cat[dat$death_tx==1 & dat$time_death_tx>30]<-"30+"
dat$trm_time_table_cat<-as.factor(dat$trm_time_table_cat)
label(dat$trm_time_table_cat)<-"TRM cat"

dat$oth_time_table_cat<-NA
dat$oth_time_table_cat[dat$death_tx==2 & dat$time_death_tx<=30]<-"0-30"
dat$oth_time_table_cat[dat$death_tx==2 & dat$time_death_tx>30]<-"30+"
dat$oth_time_table_cat<-as.factor(dat$oth_time_table_cat)
label(dat$oth_time_table_cat)<-"oth cat"

dat$other_time_table<-NA
dat$other_time_table[dat$death_tx==2]<-dat$time_death_tx[dat$death_tx==2]
label(dat$other_time_table)<-"Days to Other Mortality"

trm_table<-final_table(dat,c('trm_tab','trm_time_table','other_time_table','oth_time_table_cat'), 
                       dat$tx_type,margin=2,single=T,2,col.names=T, summary.stat='median')
trm_table

################Time-to-TRM Models###############
######Pre_tx Glucose Variability:

dat.pre<-subset(dat,!is.na(dat$cv_glucose_pre0_log2))
dat.pre<-subset(dat.pre,dat.pre$cv_glucose_pre0_log2!="-Inf")
dat.pre$tx_type<-factor(dat.pre$tx_type,levels=c('Autologous','Allogenic'))
cov_pre <- model.matrix(~ cv_glucose_pre0_log2  + post_tx_steroids+gvhd_new+
                          tx_type
                          
                        ,
                        data = dat.pre)[, -1]

#dat.covs$cv_ster<-dat.covs$cv_glucose_preinf_pre0*dat.covs$post_tx_steroids
multi_trm_pre <- crr(dat.pre$time_death_tx, dat.pre$death_tx, cov_pre, failcode=1, cencode=0,
                     na.action=na.omit, variance=TRUE) 
summary(multi_trm_pre)

trm_pre<-fx_model_surv(multi_trm_pre,cov_pre)

trm_pre.plot1<-multi_trm_pre


dat.0_30<-subset(dat,!is.na(dat$cv_glucose_0_30_log2))
dat.0_30<-subset(dat.0_30,!(dat.0_30$time_death<=30 & dat.0_30$death_tx==2))
dat.0_30$tx_type<-factor(dat.0_30$tx_type,levels=c('Autologous','Allogenic'))

cov_0_30 <- model.matrix(~ cv_glucose_0_30_log2+ post_tx_steroids+gvhd_new+
                           tx_type,
                         data = dat.0_30)[, -1]

multi_trm_0_30 <- crr(dat.0_30$time_death_tx, dat.0_30$death_tx, cov_0_30, failcode=1, cencode=0,
                      na.action=na.omit, variance=TRUE) 
trm_0_30<-fx_model_surv(multi_trm_0_30,cov_0_30)
trm_0_30
trm_0_30.plot1<-multi_trm_0_30
# 
# cov_0_100 <- model.matrix(~ cv_glucose_0_100_log2+ post_tx_steroids +
#                             gvhd_new + tx_type
#                           ,
#                           data = dat)[, -1]
# dat.0_100<-subset(dat,!is.na(dat$cv_glucose_0_100_log2))
# 
# multi_trm_0_100 <- crr(dat.0_100$time_death_tx, dat.0_100$death_tx, cov_0_100, failcode=1, cencode=0,
#                        na.action=na.omit, variance=TRUE)
# summary(multi_trm_0_100)
# 
# trm_0_100<-fx_model_surv(multi_trm_0_100,cov_0_100)



