####LYMPHOCYTES AS OUTCOME####
library(pscl)
library(lmtest)

#univarite zero-inflated negative binomial (ZINB) models:

######LYMPHS######
lymphs_perc<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                                   "vent_2",
                                                   "location_2",
                                                   "death_yn",
                                                   "age_at_diagnosis_use",
                                                   "OnDexamethasone",
                                                   "ATRT",
                                                   "Ependymoma",
                                                   "Germinoma",
                                                   "HGG",
                                                   "LGG",
                                                   "Medullo",
                                                   "OET")],zinb_mod,y=dat$CSF_LYMPHS..)),row.names=NULL)

lymphs_perc$Independent_variable<-c(lymphs_perc$Independent_variable[c(1:6)],"ATRT",
                       "Ependymoma",
                       "Germinoma",
                       "HGG",
                       "LGG",
                       "Medullo",
                       "OET")

dat.complete<-subset(dat,!is.na(dat$vent_2))
dat.complete<-subset(dat.complete,!is.na(dat.complete$location_2))
l1_full<-zeroinfl(dat.complete$CSF_LYMPHS.. ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                  OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                    Medullo,
                  dist = 'negbin',
                  data = dat.complete)

l1_forced<-zeroinfl(dat.complete$CSF_LYMPHS.. ~ 1,
                    dist = 'negbin',
                    data = dat.complete)

#step(l1_forced, scope=list(lower=l1_forced, upper=l1_full), direction="both")


#started with all significatn vars, then did backward selection:
lymphs_perc_multi <- zinb_mod_multi(zeroinfl(CSF_LYMPHS.. ~ CSF_Malignant+death_yn+ATRT+LGG+age_at_diagnosis_use,
                              dist = 'negbin',
                              data = dat))

lymphs_abs<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                                  "vent_2",
                                                  "location_2",
                                                  "death_yn",
                                                  "age_at_diagnosis_use",
                                                  "OnDexamethasone",
                                                  "ATRT",
                                                  "Ependymoma",
                                                  "Germinoma",
                                                  "HGG",
                                                  "LGG",
                                                  "Medullo",
                                                  "OET")],zinb_mod,y=dat$Lymphs.absolute)),row.names=NULL)
lymphs_abs$Independent_variable<-c(lymphs_abs$Independent_variable[c(1:6)],"ATRT",
                                    "Ependymoma",
                                    "Germinoma",
                                    "HGG",
                                    "LGG",
                                    "Medullo",
                                    "OET")

l2_full<-zeroinfl(dat.complete$Lymphs.absolute ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                    OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                    Medullo,
                  dist = 'negbin',
                  data = dat.complete)

l2_forced<-zeroinfl(dat.complete$Lymphs.absolute ~ 1,
                    dist = 'negbin',
                    data = dat.complete)

#step(l2_forced, scope=list(lower=l2_forced, upper=l2_full), direction="both")

lymphs_abs_multi <- zinb_mod_multi(zeroinfl(Lymphs.absolute ~ LGG+age_at_diagnosis_use+location_2+vent_2+ATRT+CSF_Malignant,
                                             dist = 'negbin',
                                             data = dat))

######MONOS######

monos_perc<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                                  "vent_2",
                                                  "location_2",
                                                  "death_yn",
                                                  "age_at_diagnosis_use",
                                                  "OnDexamethasone",
                                                  "ATRT",
                                                  "Ependymoma",
                                                  "Germinoma",
                                                  "HGG",
                                                  "LGG",
                                                  "Medullo",
                                                  "OET")],zinb_mod,y=dat$CSF_MONOCYTES..)),row.names=NULL)
monos_perc$Independent_variable<-c(monos_perc$Independent_variable[c(1:6)],"ATRT",
                                   "Ependymoma",
                                   "Germinoma",
                                   "HGG",
                                   "LGG",
                                   "Medullo",
                                   "OET")


mo1_full<-zeroinfl(dat.complete$CSF_MONOCYTES.. ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                    OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                     Medullo,
                  dist = 'negbin',
                  data = dat.complete)

mo1_forced<-zeroinfl(dat.complete$CSF_MONOCYTES.. ~ 1,
                    dist = 'negbin',
                    data = dat.complete)

#step(mo1_forced, scope=list(lower=mo1_forced, upper=mo1_full), direction="both")

monos_perc_multi <- zinb_mod_multi(zeroinfl(CSF_MONOCYTES.. ~ LGG+age_at_diagnosis_use+CSF_Malignant+HGG+Ependymoma,
                                            dist = 'negbin',
                                            data = dat))

monos_abs<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                                 "vent_2",
                                                 "location_2",
                                                 "death_yn",
                                                 "age_at_diagnosis_use",
                                                 "OnDexamethasone",
                                                 "ATRT",
                                                 "Ependymoma",
                                                 "Germinoma",
                                                 "HGG",
                                                 "LGG",
                                                 "Medullo",
                                                 "OET")],zinb_mod,y=dat$Monos.absolute)),row.names=NULL)
monos_abs$Independent_variable<-c(monos_abs$Independent_variable[c(1:6)],"ATRT",
                                   "Ependymoma",
                                   "Germinoma",
                                   "HGG",
                                   "LGG",
                                   "Medullo",
                                   "OET")


mo2_full<-zeroinfl(dat.complete$Monos.absolute ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                     OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                     Medullo,
                   dist = 'negbin',
                   data = dat.complete)

mo2_forced<-zeroinfl(dat.complete$Monos.absolute ~ 1,
                     dist = 'negbin',
                     data = dat.complete)

#step(mo2_forced, scope=list(lower=mo2_forced, upper=mo2_full), direction="both")

monos_abs_multi <- zinb_mod_multi(zeroinfl(Monos.absolute ~ LGG+location_2+vent_2+OnDexamethasone+
                                             death_yn+CSF_Malignant,
                                            dist = 'negbin',
                                            data = dat))



######SEGS - WON'T CONVERGE######
#logistic regression instead:
segs_perc<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                                "vent_2",
                                                "location_2",
                                                "death_yn",
                                                "age_at_diagnosis_use",
                                                "OnDexamethasone",
                                                "ATRT",
                                                "Ependymoma",
                                                "Germinoma",
                                                "HGG",
                                                "LGG",
                                                "Medullo",
                                                "OET")],log_mod,y=dat$SEGS.absolute_zero)),row.names=NULL)

segs_perc$Parameter<-c(segs_perc$Parameter[c(1:6)],"ATRT",
                      "Ependymoma",
                      "Germinoma",
                      "HGG",
                      "LGG",
                      "Medullo",
                      "OET")


segs1_full<-glm(dat.complete$SEGS.absolute_zero ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                  OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                  Medullo,
                family = 'binomial',
                data = dat.complete)

segs1_forced<-glm(dat.complete$SEGS.absolute_zero ~ 1,
                  family = 'binomial',
                  data = dat.complete)

#step(segs1_forced, scope=list(lower=segs1_forced, upper=segs1_full), direction="both")

segs_perc_multi <- fx_model_logistic(y="SEGS.absolute_zero",x=c("Medullo","LGG"),data=
                                      dat.complete)

segs_abs<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                               "vent_2",
                                               "location_2",
                                               "death_yn",
                                               "age_at_diagnosis_use",
                                               "OnDexamethasone",
                                               "ATRT",
                                               "Ependymoma",
                                               "Germinoma",
                                               "HGG",
                                               "LGG",
                                               "Medullo",
                                               "OET")],log_mod,y=dat$CSF_SEGS_zero)),row.names=NULL)

segs_abs$Parameter<-c(segs_abs$Parameter[c(1:6)],"ATRT",
                     "Ependymoma",
                     "Germinoma",
                     "HGG",
                     "LGG",
                     "Medullo",
                     "OET")


segs2_full<-glm(dat.complete$CSF_SEGS_zero ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                 OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                 Medullo,
               family = 'binomial',
               data = dat.complete)

segs2_forced<-glm(dat.complete$CSF_SEGS_zero ~ 1,
                 family = 'binomial',
                 data = dat.complete)

#step(segs2_forced, scope=list(lower=segs2_forced, upper=segs2_full), direction="both")

segs_abs_multi <- fx_model_logistic(y="CSF_SEGS_zero",x=c("Medullo","OnDexamethasone","ATRT"),data=
                                     dat.complete)
######MACROS######

macros_perc<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                                  "vent_2",
                                                  "location_2",
                                                  "death_yn",
                                                  "age_at_diagnosis_use",
                                                  "OnDexamethasone",
                                                  "ATRT",
                                                  "Ependymoma",
                                                  "Germinoma",
                                                  "HGG",
                                                  "LGG",
                                                  "Medullo",
                                                  "OET")],zinb_mod,y=dat$CSF_MACROPHAGES..)),row.names=NULL)

macros_perc$Independent_variable<-c(macros_perc$Independent_variable[c(1:6)],"ATRT",
                                   "Ependymoma",
                                   "Germinoma",
                                   "HGG",
                                   "LGG",
                                   "Medullo",
                                   "OET")


ma1_full<-zeroinfl(dat.complete$CSF_MACROPHAGES.. ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                     OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                     Medullo,
                   dist = 'negbin',
                   data = dat.complete)

ma1_forced<-zeroinfl(dat.complete$CSF_MACROPHAGES.. ~ 1,
                     dist = 'negbin',
                     data = dat.complete)

#step(ma1_forced, scope=list(lower=ma1_forced, upper=ma1_full), direction="both")

macros_perc_multi <- zinb_mod_multi(zeroinfl(CSF_MACROPHAGES.. ~ age_at_diagnosis_use+LGG+Ependymoma+HGG,
                                            dist = 'negbin',
                                            data = dat))

macros_abs<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                                 "vent_2",
                                                 "location_2",
                                                 "death_yn",
                                                 "age_at_diagnosis_use",
                                                 "OnDexamethasone",
                                                 "ATRT",
                                                 "Ependymoma",
                                                 "Germinoma",
                                                 "HGG",
                                                 "LGG",
                                                 "Medullo",
                                                 "OET")],zinb_mod,y=dat$Macros.absolute)),row.names=NULL)
macros_abs$Independent_variable<-c(macros_abs$Independent_variable[c(1:6)],"ATRT",
                                  "Ependymoma",
                                  "Germinoma",
                                  "HGG",
                                  "LGG",
                                  "Medullo",
                                  "OET")


ma2_full<-zeroinfl(dat.complete$Macros.absolute ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                     OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                     Medullo,
                   dist = 'negbin',
                   data = dat.complete)

ma2_forced<-zeroinfl(dat.complete$Macros.absolute ~ 1,
                     dist = 'negbin',
                     data = dat.complete)

#step(ma2_forced, scope=list(lower=ma2_forced, upper=ma2_full), direction="both")

macros_abs_multi <- zinb_mod_multi(zeroinfl(Macros.absolute ~ Germinoma+LGG+HGG+OnDexamethasone+vent_2,
                                           dist = 'negbin',
                                           data = dat))


######TOTAL - WON'T CONVERGE######
#do logisitic regression instead?
wbc_abs<-data.frame(do.call(rbind,lapply(dat[c("CSF_Malignant",
                                                  "vent_2",
                                                  "location_2",
                                                  "death_yn",
                                                  "age_at_diagnosis_use",
                                                  "OnDexamethasone",
                                                  "ATRT",
                                                  "Ependymoma",
                                                  "Germinoma",
                                                  "HGG",
                                                  "LGG",
                                                  "Medullo",
                                                  "OET")],log_mod,y=dat$CSF_NUCLEATED_CELLS_zero)),row.names=NULL)

wbc_abs$Parameter<-c(wbc_abs$Parameter[c(1:6)],"ATRT",
                                   "Ependymoma",
                                   "Germinoma",
                                   "HGG",
                                   "LGG",
                                   "Medullo",
                                   "OET")


wbc2_full<-glm(dat.complete$CSF_NUCLEATED_CELLS_zero ~ CSF_Malignant+vent_2+location_2+death_yn+age_at_diagnosis_use+
                     OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                     Medullo,
                   family = 'binomial',
                   data = dat.complete)

wbc2_forced<-glm(dat.complete$CSF_NUCLEATED_CELLS_zero ~ 1,
                 family = 'binomial',
                     data = dat.complete)

#step(wbc2_forced, scope=list(lower=wbc2_forced, upper=wbc2_full), direction="both")

wbc_abs_multi <- fx_model_logistic(y="CSF_NUCLEATED_CELLS_zero",x=c("LGG","vent_2","ATRT","location_2"),data=
                                     dat.complete)
 
#### MALIGNANT CELLS #### 
mal_abs<-data.frame(do.call(rbind,lapply(dat[c("vent_2",
                                               "location_2",
                                               "death_yn",
                                               "age_at_diagnosis_use",
                                               "OnDexamethasone",
                                               "ATRT",
                                               "Ependymoma",
                                               "Germinoma",
                                               "HGG",
                                               "LGG",
                                               "Medullo",
                                               "OET")],log_mod,y=dat$CSF_Malignant)),row.names=NULL)

mal_abs$Parameter<-c(mal_abs$Parameter[c(1:5)],"ATRT",
                     "Ependymoma",
                     "Germinoma",
                     "HGG",
                     "LGG",
                     "Medullo",
                     "OET")


mal2_full<-glm(dat.complete$CSF_Malignant ~ vent_2+location_2+death_yn+age_at_diagnosis_use+
                 OnDexamethasone+ATRT+Ependymoma+Germinoma+HGG+LGG+
                 OET,
               family = 'binomial',
               data = dat.complete)

mal2_forced<-glm(dat.complete$CSF_Malignant ~ 1,
                 family = 'binomial',
                 data = dat.complete)

#step(mal2_forced, scope=list(lower=mal2_forced, upper=mal2_full), direction="both")

mal_abs_multi <- fx_model_logistic(y="CSF_Malignant",x=c("location_2","Ependymoma","vent_2",
                                                         "death_yn","Germinoma","LGG"),data=
                                     dat.complete)
#dat.complete$csf_mal<-ifelse(dat.complete$CSF_Malignant=="Malignant ",1,0)

#zero inflated poisson model:
#check percentage of zeros
# 100*sum(dat$CSF_SEGS.. == 0)/nrow(dat) #11%
# 
# #is mean close to variance? then can use ZIP, otherwise should use ZINB
# mean(dat$CSF_SEGS..)
# var(dat$CSF_SEGS..) #MUCH larger, so should use ZINB, but here we compare both choices:
# 
# M3<-zeroinfl(dat$CSF_SEGS..~dat$CSF_Malignant+dat$vent_2+dat$Dx.Category+dat$death_yn
#              ,data=dat)
# # Dispersion statistic
# E2 <- resid(M3, type = "pearson")
# N  <- nrow(dat)
# p  <- length(coef(M3))
# sum(E2^2) / (N - p)
# 
# #try a zero inflated negative binomial model:
# M4 <- zeroinfl(CSF_NUCLEATED_CELLS ~ dat$CSF_Malignant+dat$vent_2+dat$Dx.Category+dat$death_yn,
#                dist = 'negbin',
#                data = dat)
# summary(M4)
# # Dispersion Statistic
# E2 <- resid(M4, type = "pearson")
# N  <- nrow(dat)
# p  <- length(coef(M4)) + 1 # '+1' is due to theta
# sum(E2^2) / (N - p)
# lrtest(M3,M4)
# 
# AIC(lymph_perc)
# summary(lymph_perc)

# ####build final lymph model:
# model_lymph_1 <- zeroinfl(CSF_LYMPHS.. ~ 1|dat$CSF_Malignant,
#                           dist = 'negbin',
#                           data = dat)
# summary(model_lymph_1)
# AIC(model_lymph_1)
# 
# model_lymph_2 <- zeroinfl(CSF_LYMPHS.. ~ dat$vent_2|1,
#                           dist = 'negbin',
#                           data = dat)
# summary(model_lymph_2)
# AIC(model_lymph_2)
# 
# model_lymph_3 <- zeroinfl(CSF_LYMPHS.. ~ dat$death_yn|1,
#                           dist = 'negbin',
#                           data = dat)
# summary(model_lymph_3)
# AIC(model_lymph_3)
# 
# model_lymph_4 <- zeroinfl(CSF_LYMPHS.. ~ dat$Dx.Category,
#                           dist = 'negbin',
#                           data = dat)
# summary(model_lymph_4)
# AIC(model_lymph_4)
# 
# #########FINAL MODEL FOR LYMPHOCYTE PERCENTAGE#########
# model_lymph_full <- zeroinfl(CSF_LYMPHS.. ~ dat$Dx.Category+dat$death_yn+dat$vent_2|dat$Dx.Category+dat$CSF_Malignant,
#                              dist = 'negbin',
#                              data = dat)
# AIC(model_lymph_full)
# summary(model_lymph_full)
# anova(model_lymph_full)
# 
# model_lymph_sub <- zeroinfl(CSF_LYMPHS.. ~ dat$Dx.Category+dat$death_yn+dat$vent_2|dat$Dx.Category+dat$CSF_Malignant,
#                             dist = 'negbin',
#                             data = dat)
# lrtest(model_lymph_sub,model_lymph_full) #LRtest<0.05, so need the more complex model
# 
# ####MONOCYTES AS OUTCOME####
# 
# #zero inflated poisson model:
# #check percentage of zeros
# 100*sum(dat$CSF_MONOCYTES.. == 0)/nrow(dat) #15%
# 
# #is mean close to variance? then can use ZIP, otherwise should use ZINB
# mean(dat$CSF_MONOCYTES..)
# var(dat$CSF_MONOCYTES..) #MUCH larger, so should use ZINB, but here we compare both choices:
# 
# M3<-zeroinfl(dat$CSF_MONOCYTES..~dat$CSF_Malignant
#              ,data=dat)
# # Dispersion statistic
# E2 <- resid(M3, type = "pearson")
# N  <- nrow(dat)
# p  <- length(coef(M3))  
# sum(E2^2) / (N - p)
# 
# #try a zero inflated negative binomial model:
# M4 <- zeroinfl(CSF_MONOCYTES.. ~ dat$CSF_Malignant,
#                dist = 'negbin',
#                data = dat)
# summary(M4)
# # Dispersion Statistic
# E2 <- resid(M4, type = "pearson")
# N  <- nrow(dat)
# p  <- length(coef(M4)) + 1 # '+1' is due to theta
# sum(E2^2) / (N - p)
# lrtest(M3,M4)
# 
# ####build final lymph model:
# model_mono_1 <- zeroinfl(CSF_MONOCYTES.. ~ 1|dat$CSF_Malignant,
#                          dist = 'negbin',
#                          data = dat)
# summary(model_mono_1)
# AIC(model_mono_1)
# 
# model_mono_2 <- zeroinfl(CSF_MONOCYTES.. ~ 1|vent_2,
#                          dist = 'negbin',
#                          data = dat)
# summary(model_mono_2)
# AIC(model_mono_2)
# 
# model_mono_3 <- zeroinfl(CSF_MONOCYTES.. ~ dat$death_yn|1,
#                          dist = 'negbin',
#                          data = dat)
# summary(model_mono_3)
# AIC(model_mono_3)
# 
# model_mono_4 <- zeroinfl(CSF_MONOCYTES.. ~ dat$Dx.Category,
#                          dist = 'negbin',
#                          data = dat)
# summary(model_mono_4)
# AIC(model_mono_4)
# 
# #########FINAL MODEL FOR MONOCYTE PERCENTAGE#########
# model_mono_full <- zeroinfl(CSF_MONOCYTES.. ~ dat$Dx.Category|dat$Dx.Category+dat$CSF_Malignant,
#                             dist = 'negbin',
#                             data = dat)
# AIC(model_mono_full)
# summary(model_mono_full)
# anova(model_mono_full)
# 
# dat.sub<-subset(dat,!is.na(dat$vent_2))
# model_mono_sub <- zeroinfl(CSF_MONOCYTES.. ~ Dx.Category|Dx.Category+CSF_Malignant,
#                            dist = 'negbin',
#                            data = dat.sub)
# lrtest(model_mono_sub,model_mono_full) #LRtest<0.05, so need the more complex model
# 
# ####SEGS AS OUTCOME####
# 
# #zero inflated poisson model:
# #check percentage of zeros
# 100*sum(dat$CSF_SEGS.. == 0)/nrow(dat) #44%
# 
# #is mean close to variance? then can use ZIP, otherwise should use ZINB
# mean(dat$CSF_SEGS..)
# var(dat$CSF_SEGS..) #MUCH larger, so should use ZINB, but here we compare both choices:
# 
# M3<-zeroinfl(dat$CSF_SEGS..~dat$CSF_Malignant
#              ,data=dat)
# # Dispersion statistic
# E2 <- resid(M3, type = "pearson")
# N  <- nrow(dat)
# p  <- length(coef(M3))  
# sum(E2^2) / (N - p)
# 
# #try a zero inflated negative binomial model:
# M4 <- zeroinfl(CSF_SEGS.. ~ dat$CSF_Malignant,
#                dist = 'negbin',
#                data = dat)
# summary(M4)
# # Dispersion Statistic
# E2 <- resid(M4, type = "pearson")
# N  <- nrow(dat)
# p  <- length(coef(M4)) + 1 # '+1' is due to theta
# sum(E2^2) / (N - p)
# lrtest(M3,M4)
# 
# ####build final segs model:
# model_mono_1 <- zeroinfl(CSF_SEGS.. ~ dat$CSF_Malignant|1,
#                          dist = 'negbin',
#                          data = dat)
# summary(model_mono_1)
# AIC(model_mono_1)
# 
# model_mono_2 <- zeroinfl(CSF_SEGS.. ~ 1|vent_2,
#                          dist = 'negbin',
#                          data = dat)
# summary(model_mono_2)
# AIC(model_mono_2)
# 
# model_mono_3 <- zeroinfl(CSF_SEGS.. ~ dat$death_yn|1,
#                          dist = 'negbin',
#                          data = dat)
# summary(model_mono_3)
# AIC(model_mono_3)
# 
# model_mono_4 <- zeroinfl(CSF_SEGS.. ~dat$Dx.Category|1,
#                          dist = 'negbin',
#                          data = dat)
# summary(model_mono_4)
# AIC(model_mono_4)
# 
# #########FINAL MODEL FOR SEGS PERCENTAGE#########
# model_mono_full <- zeroinfl(CSF_SEGS.. ~ dat$Dx.Category+dat$death_yn+dat$CSF_Malignant|dat$vent_2,
#                             dist = 'negbin',
#                             data = dat)
# AIC(model_mono_full)
# summary(model_mono_full)
# anova(model_mono_full)
# 
# dat.sub<-subset(dat,!is.na(dat$vent_2))
# model_mono_sub <- zeroinfl(CSF_SEGS.. ~ 1,
#                            dist = 'negbin',
#                            data = dat.sub)
# lrtest(model_mono_sub,model_mono_full) #LRtest>0.05, so no need for model