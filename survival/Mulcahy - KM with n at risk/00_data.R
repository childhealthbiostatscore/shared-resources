#########Mulcahy LGG manuscript#########
library(Hmisc)
library(survival)
library(jskm)
library(gridExtra)
source('S:/Shared Material/Shared Code/R/temp_table1.R')
dat<-read.csv("K:/CCBD/Mulcahy/LGG/Data/Cleaned data_3132020.csv",na.strings=c("","NA",NA))



####change as of 7/23: remove all grade 3/4:
dat<-subset(dat,dat$Grade %in% c(1,2))


###Table 1:
label(dat$Age.at.Diagnosis)<-"Age at dx, years"
dat$Gender<-as.factor(dat$Gender)
label(dat$Gender)<-"Gender"

dat$Diagnosis[dat$Diagnosis=="ASTROCYTOMA "]<-"ASTROCYTOMA"
dat$Diagnosis<-factor(dat$Diagnosis)
dat$Diagnosis<-as.factor(dat$Diagnosis)
label(dat$Diagnosis)<-"Diagnosis"

dat$Grade<-as.factor(dat$Grade)
label(dat$Grade)<-"Grade"

dat$Location.of.Tumor<-as.factor(dat$Location.of.Tumor)
label(dat$Location.of.Tumor)<-"Region of Spine"

dat$Surgical.Outcome<-as.factor(dat$Surgical.Outcome)
label(dat$Surgical.Outcome)<-"Treatment (Surgery)"

dat$Chemo.Diag.Y.N.<-as.factor(dat$Chemo.Diag.Y.N.)
label(dat$Chemo.Diag.Y.N.)<-"Treatment (Chemotherapy)"

dat$Radiation.Diag..Y.N.<-as.factor(dat$Radiation.Diag..Y.N.)
label(dat$Radiation.Diag..Y.N.)<-"Treatment (Radiation)"

dat$Relapse.if.Applicable.Y.or.N.[dat$Relapse.if.Applicable.Y.or.N.=="NO "]<-"NO"
dat$Relapse.if.Applicable.Y.or.N.[dat$Relapse.if.Applicable.Y.or.N.=="YES   "]<-"YES"
dat$Relapse.if.Applicable.Y.or.N.<-factor(dat$Relapse.if.Applicable.Y.or.N.)
label(dat$Relapse.if.Applicable.Y.or.N.)<-"Relapse"

dat$Death<-as.factor(dat$Death)
label(dat$Death)<-"Overall Survival"

dat$Months.to.Relapse[dat$Months.to.Relapse=="NO"]<-NA
dat$Months.to.Relapse<-as.numeric(as.character(dat$Months.to.Relapse))
label(dat$Months.to.Relapse)<-"Months to Relapse"

dat$Months.to.Death.Followup<-as.numeric(as.character(dat$Months.to.Death.Followup))
label(dat$Months.to.Death.Followup)<-"Months to Death/Followup"

dat$BRAF.fusion<-as.factor(dat$BRAF.fusion)
label(dat$BRAF.fusion)<-"BRAF Fusion"

##new, added on 3/9:
dat$time_to_death<-NA
dat$time_to_death[dat$Death=="YES"]<-dat$Months.to.Death.Followup[dat$Death=="YES"]
dat$time_to_death[dat$Death=="NO"]<-dat$Months.to.Death.Followup[dat$Death=="NO"]
label(dat$time_to_death)<-"Months to Death/Followup"

dat$time_to_death_yes<-NA
dat$time_to_death_yes[dat$Death=="YES"]<-dat$Months.to.Death.Followup[dat$Death=="YES"]
label(dat$time_to_death_yes)<-"Months to Death, if died"

dat$time_to_rel<-NA
dat$time_to_rel[dat$Relapse.if.Applicable.Y.or.N.=="YES"]<-dat$Months.to.Relapse[dat$Relapse.if.Applicable.Y.or.N.=="YES"]
dat$time_to_rel[dat$Relapse.if.Applicable.Y.or.N.=="NO"]<-dat$Months.to.Death.Followup[dat$Relapse.if.Applicable.Y.or.N.=="NO"]
label(dat$time_to_rel)<-"Months to Relapse/Followup"

dat$time_to_rel_yes<-NA
dat$time_to_rel_yes[dat$Relapse.if.Applicable.Y.or.N.=="YES"]<-dat$Months.to.Relapse[dat$Relapse.if.Applicable.Y.or.N.=="YES"]
label(dat$time_to_rel_yes)<-"Months to Relapse, if relapsed"

dat$progression<-NA
dat$progression[dat$Relapse.if.Applicable.Y.or.N.=="YES"]<-1
dat$progression[dat$Relapse.if.Applicable.Y.or.N.=="NO" & dat$Death=="NO"]<-0
dat$progression[dat$Relapse.if.Applicable.Y.or.N.=="NO" & dat$Death=="YES"]<-1
label(dat$progression)<-"Progression"
dat$time_to_prog<-NA
dat$time_to_prog[dat$progression==1]<-dat$Months.to.Relapse[dat$progression==1]
dat$time_to_prog[dat$progression==0]<-dat$Months.to.Death.Followup[dat$progression==0]
dat$time_to_prog[dat$Relapse.if.Applicable.Y.or.N.=="NO" & dat$Death=="YES"]<-dat$Months.to.Death.Followup[dat$Relapse.if.Applicable.Y.or.N.=="NO" & dat$Death=="YES"]
label(dat$time_to_prog)<-"Time to progression"
dat$time_to_prog_yes<-NA
dat$time_to_prog_yes[dat$progression==1]<-dat$Months.to.Relapse[dat$progression==1]
dat$time_to_prog_yes[dat$Relapse.if.Applicable.Y.or.N.=="NO" & dat$Death=="YES"]<-dat$Months.to.Death.Followup[dat$Relapse.if.Applicable.Y.or.N.=="NO" & dat$Death=="YES"]
label(dat$time_to_prog_yes)<-"Time to progression, if progressed"

# dat$progression_free<-as.factor(dat$progression_free)
# label(dat$progression_free)<-"Progression Free Survival"


####TABLE 1:

#abstract:
table(dat$death_yn[dat$BRAF.fusion=="YES"])
quantile(dat$time_to_death[dat$BRAF.fusion=="YES"],na.rm=T)
table(dat$death_yn[dat$BRAF.fusion=="NO"])
quantile(dat$time_to_death[dat$BRAF.fusion=="NO"],na.rm=T)

#methods:
quantile(dat$Age.at.Diagnosis[dat$grade_cat=="I/II"])

tab.1<-final_table(dat,c('Age.at.Diagnosis','Gender','Diagnosis','Location.of.Tumor',
                         'Surgical.Outcome','Chemo.Diag.Y.N.','Radiation.Diag..Y.N.',
                         'Relapse.if.Applicable.Y.or.N.','time_to_rel_yes',
                         'Death','time_to_death_yes','BRAF.fusion'),
                        dat$Grade,margin=2,single=F,1,col.names=T, summary.stat='median')

###KM Curves:
tab.outcomes<-final_table(dat,c('Death','time_to_death','time_to_death_yes',
                                'Relapse.if.Applicable.Y.or.N.','Months.to.Relapse',
                                'progression','time_to_prog','time_to_prog_yes'),
                   dat$Grade,margin=2,single=F,1,col.names=T, summary.stat='median')


quantile(dat$Months.to.Death.Followup[dat$death_yn==0])
dat$Months.to.Death.Followup<-as.numeric(dat$Months.to.Death.Followup)
dat$death_yn<-0
dat$death_yn[dat$Death=="YES"]<-1
 
quantile(dat$Months.to.Death.Followup[dat$BRAF.fusion=="YES"],na.rm=T)
quantile(dat$Months.to.Death.Followup[dat$BRAF.fusion=="NO"],na.rm=T)


####Overall Survival:
#by BRAF:

#Fig 1:
os_grade <- survfit(Surv(dat$time_to_death, dat$death_yn) ~ Grade, data=dat)
summary(os_grade)
fig_1a<-jskm(os_grade, timeby=50, table=T,legend=F,
     surv.scale="percent",ystrataname = "",
     pval=T,pval.coord = c(100,0.4),pval.size=4,
     main="1A) Time-to-Death by Grade",xlab="Time (months)",ylab="Percent Survival")

pfs_grade <- survfit(Surv(dat$time_to_prog, dat$progression) ~ Grade, data=dat)
summary(pfs_grade)
fig_1b<-jskm(pfs_grade, timeby=50, ystratalabs=c("Grade I", "Grade II"),
             table=T,
     surv.scale="percent",legendposition=c(0.85,0.8),ystrataname = "",
     pval=T,pval.coord = c(100,0.4),pval.size=4,
     main="1B) Time-to-Progression by Grade",
     xlab="Time (months)",ylab="Percent Progression Free")

jpeg("K:/CCBD/Mulcahy/LGG/Results/fig1_bygrade_07232020.jpeg",width=8, height=6, units='in',res=300)
        grid.arrange(fig_1a, fig_1b, ncol = 2)
dev.off()

os_grade <- survfit(Surv(dat$time_to_death, dat$death_yn) ~ 1, data=dat)
summary(os_grade)
fig_1a<-jskm(os_grade, timeby=50, table=T,legend=F,
             surv.scale="percent",ystrataname = "",
             pval=T,pval.coord = c(100,0.5),pval.size=4,
             main="1A) Time-to-Death by Grade",xlab="Time (months)",ylab="Percent Survival")

pfs_grade <- survfit(Surv(dat$time_to_prog, dat$progression) ~ 1, data=dat)
summary(pfs_grade)
fig_1b<-jskm(pfs_grade, timeby=50, ystratalabs="",
             table=T,
             surv.scale="percent",legend=F,ystrataname = "",
             pval=T,pval.coord = c(100,0.5),pval.size=4,
             main="1B) Time-to-Progression",
             xlab="Time (months)",ylab="Percent Progression Free")

jpeg("K:/CCBD/Mulcahy/LGG/Results/fig1_overall_07232020.jpeg",width=8, height=6, units='in',res=300)
grid.arrange(fig_1a, fig_1b, ncol = 2)
dev.off()

dat.res<-subset(dat,dat$Surgery.Y_N=="YES" & dat$Radiation.Diag..Y.N.=="NO" &
                  dat$Chemo.Diag.Y.N.=="NO")
dat.res.rad<-subset(dat,dat$Surgery.Y_N=="YES" & dat$Radiation.Diag..Y.N.=="YES" &
                      dat$Chemo.Diag.Y.N.=="NO")
dat.res.chem<-subset(dat,dat$Surgery.Y_N=="YES" & dat$Radiation.Diag..Y.N.=="NO" &
                       dat$Chemo.Diag.Y.N.=="YES")
dat.all<-subset(dat,dat$Surgery.Y_N=="YES" & dat$Radiation.Diag..Y.N.=="YES" &
                       dat$Chemo.Diag.Y.N.=="YES")
# fig2a <- survfit(Surv(dat.res$Months.to.Death.Followup, dat.res$death_yn) ~ Grade, data=dat.res)
# summary(fig2a)
# 
# fig_2a<-jskm(fig2a, timeby=50, table=T,
#      surv.scale="percent",legendposition=c(0.85,0.2),ystrataname = "",
#      pval=T,pval.coord = c(100,0.5),pval.size=4,
#      main="2A) Resection Only",xlab="Time (months)",ylab="Percent Survival",
#      linecols="black")
# 
# fig2b <- survfit(Surv(dat.res.chem$Months.to.Death.Followup, dat.res.chem$death_yn) ~
#                          Grade, data=dat.res.chem)
# summary(fig2b)
# fig_2b<-jskm(fig2b, timeby=50, table=T,ystratalabs=c("Grade I", "Grade II"),
#      surv.scale="percent",legendposition=c(0.85,0.2),ystrataname = "",
#      pval=T,pval.coord = c(100,0.5),pval.size=4,
#      main="2B) Resection and Chemotherapy",xlab="Time (months)",
#      ylab="Percent Survival")
# 
# fig2c <- survfit(Surv(dat.res.rad$Months.to.Death.Followup, dat.res.rad$death_yn) ~
#                    Grade, data=dat.res.rad)
# fig_2c<-jskm(fig2c, timeby=50, table=T,ystratalabs=c("Grade I", "Grade II"),
#      surv.scale="percent",legendposition=c(0.85,0.8),ystrataname = "",
#      pval=T,pval.coord = c(100,0.5),pval.size=4,
#      main="2C) Resection and Radiation",xlab="Time (months)",ylab="Percent Survival")
# 
# fig2d <- survfit(Surv(dat.all$Months.to.Death.Followup, dat.all$death_yn) ~
#                    Grade, data=dat.all)
# fig_2d<-jskm(fig2d, timeby=50, table=T,ystratalabs=c("Grade I/II", "Grade III/IV"),
#      surv.scale="percent",legendposition=c(0.85,0.8),ystrataname = "",
#      pval=T,pval.coord = c(100,0.5),pval.size=4,
#      main="2D)Resection, Radiation,and \n  Chemotherapy",xlab="Time (months)",
#      ylab="Percent Survival")
# 
# jpeg("K:/CCBD/Mulcahy/LGG/Results/fig2_07232020.jpeg",width=8, height=10, units='in',res=300)
#         grid.arrange(fig_2a, fig_2b, fig_2c, fig_2d, ncol = 2)
# dev.off()

fig2a <- survfit(Surv(dat.res$Months.to.Death.Followup, dat.res$death_yn) ~ 1, data=dat.res)
summary(fig2a)

fig_2a<-jskm(fig2a, timeby=50, table=T,
             surv.scale="percent",legend=F,ystrataname = "",
             pval=T,pval.coord = c(100,0.5),pval.size=4,
             main="2A) Resection Only",xlab="Time (months)",ylab="Percent Survival")

fig2b <- survfit(Surv(dat.res.chem$Months.to.Death.Followup, dat.res.chem$death_yn) ~
                         1, data=dat.res.chem)
summary(fig2b)
fig_2b<-jskm(fig2b, timeby=50, table=T,ystratalabs="",
             surv.scale="percent",legend=F,ystrataname = "",
             pval=T,pval.coord = c(100,0.5),pval.size=4,
             main="2B) Resection and Chemotherapy",xlab="Time (months)",
             ylab="Percent Survival")

fig2c <- survfit(Surv(dat.res.rad$Months.to.Death.Followup, dat.res.rad$death_yn) ~
                         1, data=dat.res.rad)
summary(fig2c)
fig_2c<-jskm(fig2c, timeby=50, table=T,ystratalabs="",
             surv.scale="percent",legend=F,ystrataname = "",
             pval=T,pval.coord = c(100,0.5),pval.size=4,
             main="2C) Resection and Radiation",xlab="Time (months)",ylab="Percent Survival")

fig2d <- survfit(Surv(dat.all$Months.to.Death.Followup, dat.all$death_yn) ~
                         1, data=dat.all)
summary(fig2d)
fig_2d<-jskm(fig2d, timeby=50, table=T,ystratalabs="",
             surv.scale="percent",legend=F,ystrataname = "",
             pval=T,pval.coord = c(100,0.5),pval.size=4,
             main="2D)Resection, Radiation,and \n  Chemotherapy",xlab="Time (months)",
             ylab="Percent Survival")

jpeg("K:/CCBD/Mulcahy/LGG/Results/fig2_overall_07232020.jpeg",width=8, height=10, units='in',res=300)
grid.arrange(fig_2a, fig_2b, fig_2c, fig_2d, ncol = 2)
dev.off()



os_braf <- survfit(Surv(dat$time_to_death, dat$death_yn) ~ BRAF.fusion, data=dat)
summary(os_braf)

jpeg("K:/CCBD/Mulcahy/LGG/Results/fig3_07232020.jpeg",width=6, height=6, units='in',res=300)
jskm(os_braf, timeby=50, ystratalabs=c("Negative", "Positive"),table=T,
     surv.scale="percent",legendposition=c(0.85,0.2),ystrataname = "Groups:",
     pval=T,pval.coord = c(100,0.5),pval.size=4,
     main="Time-to-Death by BRAF-KIAA1549",xlab="Time (months)",
     ylab="Percent Survival",linecols="Set2")
dev.off()

