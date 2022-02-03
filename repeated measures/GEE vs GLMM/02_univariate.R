###Univariate regression:
for (i in c(1:ncol(dat))){
  label(dat[,i])<-paste0(colnames(dat[i]))
}  
dat<-dat[order(dat$Subject.ID,dat$Age..years.),]
dat.1<-dat[!duplicated(dat$Subject.ID),]
covariates_all<-c("Age..years.","diagnosis_cat","Lactate...initial.within.8.hours","lactate_cut_CDE","lactate_cut_INF","lactate_2","lactate_2_4",
              "maxtemp_ED","Hypotension_ED","tachy_goldstein",
              "tachypneic_goldstein","Chills.or.rigors","URI.symptoms",
              "Stat_activation_final","WBC...initial","ANC.initial",
              "Neutropenic","AMC.initial","ALC.initial",
              "Plt...initial","Hemoglobin...initial",
              #"CRP_within8","PCT_within8","ESR_within8",
               "D.dimer_within8",#"D.dimer_within8_ISTH",
              # "Fibrinogen_within8",#"Fibrinogen_within8_ISTH",
               "PT_within8",#"PT_within8_ISTH",
              "intensity_grp")
# Would like to git a GLMM, but most of the models don't converge. Just fit logistic? Sensitity analysis of first obs per patient?

dat$CDE.within.48.hrs_01<-ifelse(dat$CDE.within.48.hrs=="yes",1,0)
dat$inv_bacteria_yn48hr_01<-ifelse(dat$inv_bacteria_yn48hr=="yes",1,0)

#secondary outcomes:
# dat$Documented.infection.within.48.hours._01<-ifelse(dat$Documented.infection.within.48.hours.=="yes",1,0)
# dat$virus_yn48hr_01<-ifelse(dat$virus_yn48hr=="yes",1,0)
#dat$fungus_yn48hr_01<-ifelse(dat$fungus_yn48hr=="yes",1,0)
#dat$parasite_yn48hr_01<-ifelse(dat$parasite_yn48hr=="yes",1,0)


CDE_uni_gee<-data.frame(do.call(rbind,lapply(dat[covariates_all],logit_gee,outcome=dat$CDE.within.48.hrs_01)),row.names=NULL)
BACT_uni_gee<-data.frame(do.call(rbind,lapply(dat[covariates_all],logit_gee,outcome=dat$inv_bacteria_yn48hr_01)),row.names=NULL)

colnames(CDE_uni_gee)<-c("Variable","Reference","OR (95% CI) - CDE","pval - CDE")
colnames(BACT_uni_gee)<-c("Variable","Reference","OR (95% CI) - Bacteria","pval - Bacteria")
CDE_uni_gee$row<-rep(1:nrow(CDE_uni_gee))

BACT_uni_gee<-BACT_uni_gee[,-2]
tab_uni<-merge(CDE_uni_gee,BACT_uni_gee,by="Variable",all=T)
tab_uni<-tab_uni[order(tab_uni$row),]
tab_uni<-tab_uni[,-5]

#overall f-test results:
######DIAGNOSIS#####
mod<-gee(CDE.within.48.hrs_01 ~ diagnosis_cat,data = dat, id = id_n, family = binomial,
         corstr = "exchangeable")
coef<-mod$coefficients
var<-mod$robust.variance #D_11, D_12, D_22 are covariance of random effects

#(count):
c.w <- matrix(c(0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,0,
                0,0,0,0,1), byrow=T, 4, 5)

c.wa <- c.w%*%coef
## chisquare test value
chisq.v <- t(c.wa)%*%solve(c.w%*% var %*%t(c.w)) %*% c.wa
pval.w <- pchisq(chisq.v, nrow(c.w), ncp=0, lower.tail=FALSE, log.p=FALSE)
p<-round(pval.w,4)
p

mod<-gee(inv_bacteria_yn48hr_01 ~ diagnosis_cat,data = dat, id = id_n, family = binomial,
         corstr = "exchangeable")
coef<-mod$coefficients
var<-mod$robust.variance #D_11, D_12, D_22 are covariance of random effects

#(count):
c.w <- matrix(c(0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,0,
                0,0,0,0,1), byrow=T, 4, 5)

c.wa <- c.w%*%coef
## chisquare test value
chisq.v <- t(c.wa)%*%solve(c.w%*% var %*%t(c.w)) %*% c.wa
pval.w <- pchisq(chisq.v, nrow(c.w), ncp=0, lower.tail=FALSE, log.p=FALSE)
p<-round(pval.w,4)
p

######LACTATE GRP#####
mod<-gee(CDE.within.48.hrs_01 ~ lactate_2_4,data = dat, id = id_n, family = binomial,
         corstr = "exchangeable")
coef<-mod$coefficients
var<-mod$robust.variance #D_11, D_12, D_22 are covariance of random effects

#(count):
c.w <- matrix(c(0,1,0,
                0,0,1), byrow=T, 2, 3)

c.wa <- c.w%*%coef
## chisquare test value
chisq.v <- t(c.wa)%*%solve(c.w%*% var %*%t(c.w)) %*% c.wa
pval.w <- pchisq(chisq.v, nrow(c.w), ncp=0, lower.tail=FALSE, log.p=FALSE)
p<-round(pval.w,4)
p

mod<-gee(inv_bacteria_yn48hr_01 ~ lactate_2_4,data = dat, id = id_n, family = binomial,
         corstr = "exchangeable")
coef<-mod$coefficients
var<-mod$robust.variance #D_11, D_12, D_22 are covariance of random effects

#(count):
c.w <- matrix(c(0,1,0,
                0,0,1), byrow=T, 2, 3)

c.wa <- c.w%*%coef
## chisquare test value
chisq.v <- t(c.wa)%*%solve(c.w%*% var %*%t(c.w)) %*% c.wa
pval.w <- pchisq(chisq.v, nrow(c.w), ncp=0, lower.tail=FALSE, log.p=FALSE)
p<-round(pval.w,4)
p

######INTENSITY GROUP#####
mod<-gee(CDE.within.48.hrs_01 ~ intensity_grp,data = dat, id = id_n, family = binomial,
         corstr = "exchangeable")
coef<-mod$coefficients
var<-mod$robust.variance #D_11, D_12, D_22 are covariance of random effects

#(count):
c.w <- matrix(c(0,1,0,
                0,0,1), byrow=T, 2, 3)

c.wa <- c.w%*%coef
## chisquare test value
chisq.v <- t(c.wa)%*%solve(c.w%*% var %*%t(c.w)) %*% c.wa
pval.w <- pchisq(chisq.v, nrow(c.w), ncp=0, lower.tail=FALSE, log.p=FALSE)
p<-round(pval.w,4)
p

mod<-gee(inv_bacteria_yn48hr_01 ~ intensity_grp,data = dat, id = id_n, family = binomial,
         corstr = "exchangeable")
coef<-mod$coefficients
var<-mod$robust.variance #D_11, D_12, D_22 are covariance of random effects

#(count):
c.w <- matrix(c(0,1,0,
                0,0,1), byrow=T, 2, 3)

c.wa <- c.w%*%coef
## chisquare test value
chisq.v <- t(c.wa)%*%solve(c.w%*% var %*%t(c.w)) %*% c.wa
pval.w <- pchisq(chisq.v, nrow(c.w), ncp=0, lower.tail=FALSE, log.p=FALSE)
p<-round(pval.w,4)
p
