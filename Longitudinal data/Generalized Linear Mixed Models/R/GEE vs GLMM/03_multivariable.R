### multivariable for CDE:
library(gee)
library(MuMIn)
#try GEEs:

dat$CDE.within.48.hrs_01<-ifelse(dat$CDE.within.48.hrs=="yes",1,0)

##########CDE MODEL#########
#full model:
dat<-subset(dat,!is.na(dat$intensity_grp))
CDE_full<-gee(CDE.within.48.hrs_01 ~
                #forced in:
                lactate_2+
                Neutropenic+
                intensity_grp+
                #options:
                Age..years.+
                Hypotension_ED+
                tachy_goldstein+
                Chills.or.rigors+
                Stat_activation_final,
                  data = dat, id = id_n, family = binomial,corstr = "exchangeable",na.action="na.fail")
summary(CDE_full)
QIC(CDE_full)
model.sel(CDE_full,rank="QIC")
CDE_select<-dredge(CDE_full,beta="none",rank="QIC",subset=~lactate_2&Neutropenic&intensity_grp)

CDE_final<-gee(CDE.within.48.hrs_01 ~
                #forced in:
                lactate_2+
                Neutropenic+
                intensity_grp+
                #options:
                Age..years.+
                Hypotension_ED+
                tachy_goldstein+
                Stat_activation_final,
              data = dat, id = id_n, family = binomial,corstr = "exchangeable",na.action="na.fail")
QIC(CDE_final)

cde_sum<-data.frame(summary(CDE_final)$coefficients)
cde_sum$p<-2 * pnorm(abs(cde_sum$Robust.z), lower.tail = FALSE)
cde_sum$Var<-row.names(cde_sum)
cde_sum<-cde_sum[-1,-c(2,3,5)]
cde_sum$p<-round(cde_sum$p,4)
cde_sum$p[cde_sum$p<0.001]<-"<0.001"
cde_sum$OR<-paste0(round(exp(cde_sum[,1]),2)," 95% CI: (",round(exp(cde_sum[,1]-1.96*cde_sum[,2]),2),
                   ", ",round(exp(cde_sum[,1]+1.96*cde_sum[,2]),2),")")

CDE_sum<-cde_sum[,c(4,5,3)]

#intensity grp overall p-val:
coef<-CDE_final$coefficients
var<-CDE_final$robust.variance #D_11, D_12, D_22 are covariance of random effects

#(count):
c.w <- matrix(c(0,0,0,1,0,0,0,0,0,
                0,0,0,0,1,0,0,0,0), byrow=T, 2, 9)

c.wa <- c.w%*%coef
## chisquare test value
chisq.v <- t(c.wa)%*%solve(c.w%*% var %*%t(c.w)) %*% c.wa
pval.w <- pchisq(chisq.v, nrow(c.w), ncp=0, lower.tail=FALSE, log.p=FALSE)
p<-round(pval.w,4)
p


#######Bacterial infection
INF_full<-gee(inv_bacteria_yn48hr_01 ~
                #forced in:
                lactate_2+
                Neutropenic+
                intensity_grp+
                #options:
                Hypotension_ED+
                tachy_goldstein+
                Chills.or.rigors+
                Stat_activation_final,
              data = dat, id = id_n, family = binomial,corstr = "exchangeable",na.action="na.fail")
summary(INF_full)
QIC(INF_full)
INF_select<-dredge(INF_full,beta="none",rank="QIC",subset=~lactate_2&Neutropenic&intensity_grp)


INF_final<-gee(inv_bacteria_yn48hr_01 ~
                #forced in:
                lactate_2+
                Neutropenic+
                intensity_grp+
                #options:
                Hypotension_ED+
                Chills.or.rigors+
                Stat_activation_final,
              data = dat, id = id_n, family = binomial,corstr = "exchangeable",na.action="na.fail")
QIC(INF_final)
inf_sum<-data.frame(summary(INF_final)$coefficients)
inf_sum$p<-2 * pnorm(abs(inf_sum$Robust.z), lower.tail = FALSE)
inf_sum$Var<-row.names(inf_sum)
inf_sum<-inf_sum[-1,-c(2,3,5)]
inf_sum$p<-round(inf_sum$p,4)
inf_sum$p[inf_sum$p<0.001]<-"<0.001"
inf_sum$OR<-paste0(round(exp(inf_sum[,1]),2)," 95% CI: (",round(exp(inf_sum[,1]-1.96*inf_sum[,2]),2),
                   ", ",round(exp(inf_sum[,1]+1.96*inf_sum[,2]),2),")")

INF_sum<-inf_sum[,c(4,5,3)]
