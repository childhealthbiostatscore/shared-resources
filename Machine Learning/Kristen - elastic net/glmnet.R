###########LASSO/ELASTIC NET REGRESSION - FINAL VERSION###########
library(ensr)

####CLUSTER NUMBER####
lasso.clust <- dat.apcd[, which(colnames(dat.apcd) %in% 
                           c("clustnum",
                             "age_at_mv","Gender","race_cat",
                             "Score","pim2_conv","ccc_final",
                             "admt_source","cpr_prehosp","severity_of_ards_3",
                             "ards_mod_sev_total","corr_mv_time_days",
                             "PICULOS.days.","ECMO_yn","HFOV_yn","cpr_inhosp_yn",
                             "trac_yn","complicated_icu",
                             "Prim.Respiratory","Prim.Injury","Prim.Neurologic",
                             "Prim.Infectious","Prim.Gastrointestinal",
                             "Prim.Cardiovascular","prim.poinsoning","Prim.Oncologic"))]

#LASSO
x_clust <- model.matrix(clustnum~.,lasso.clust)
y_clust <- as.matrix(lasso.clust$clustnum)

lasso_clust <- cv.glmnet(x_clust,y_clust,
                    alpha=1,
                    family="multinomial",
                    type.multinomial="ungrouped",
                    type.measure ="mse" )
plot(lasso_clust)
#min value of lambda
lambda_min <- lasso_clust$lambda.min
coef(lasso_clust,s=lambda_min)
#best value of lambda
lambda_1se <- lasso_clust$lambda.1se
#regression coefficients
coef(lasso_clust,s=lambda_1se) 

#ELASTIC NET

en_clust <-   ensr(y = y_clust, x = x_clust,
                       standardize = FALSE, 
                       standardize.response = FALSE,
                       #foldid = foldid,
                       family = "multinomial",
                       type.multinomial = "ungrouped")

plot(en_clust)
summary(en_clust)[cvm == min(cvm)]
str(preferable(en_clust), max.level = 1L)
multi_best <- coef(en_clust[[17]], s = 0.01028492)

summary(en_clust)[, .SD[cvm == min(cvm)], by = nzero]
plot(en_clust, type = 2)  

summary(en_clust)[nzero %in% c(4)] [, .SD[cvm == min(cvm)], by = nzero]
summary(en_clust)[nzero %in% c(5)] [, .SD[cvm == min(cvm)], by = nzero]
summary(en_clust)[nzero %in% c(6)] [, .SD[cvm == min(cvm)], by = nzero]

multi_4 <- coef(en_clust[[17]], s = 0.01507097)
multi_5 <- coef(en_clust[[17]], s = 0.01096247)
multi_6 <- coef(en_clust[[2]], s = 1.827517)


####CLAIM TRAJECTORY####
dat.claim$claim_trend_3<-NA
dat.claim$claim_trend_3[dat.claim$claim_trend=="Low"]<-0
dat.claim$claim_trend_3[dat.claim$claim_trend=="Decreasing"]<-1
dat.claim$claim_trend_3[dat.claim$claim_trend=="High"]<-2
dat.claim$claim_trend_3[dat.claim$claim_trend=="Increasing"]<-2
lasso.claim <- dat.claim[, which(colnames(dat.claim) %in% 
                                  c("claim_trend_3",
                                    "age_at_mv","Gender","race_cat",
                                    "Score","pim2_conv","ccc_final",
                                    "admt_source","cpr_prehosp","severity_of_ards_3",
                                    "ards_mod_sev_total","corr_mv_time_days",
                                    "PICULOS.days.","ECMO_yn","HFOV_yn","cpr_inhosp_yn",
                                    "trac_yn","complicated_icu",
                                    "Prim.Respiratory","Prim.Injury","Prim.Neurologic",
                                    "Prim.Infectious","Prim.Gastrointestinal",
                                    "Prim.Cardiovascular","prim.poinsoning","Prim.Oncologic"
                                  ))]

#LASSO
x_claim <- model.matrix(claim_trend_3~.,lasso.claim)
y_claim <- lasso.claim$claim_trend_3

lasso_claim <- cv.glmnet(x_claim,y_claim,
                         alpha=1,
                         family="multinomial",
                         type.multinomial="ungrouped",
                         type.measure ="mse" )
plot(lasso_claim)
#min value of lambda
lambda_min <- lasso_claim$lambda.min
coef(lasso_claim,s=lambda_min) #
#best value of lambda
lambda_1se <- lasso_claim$lambda.1se
#regression coefficients
coef(lasso_claim,s=lambda_1se) 

#ELASTIC NET
en_claim <-   ensr(y = y_claim, x = x_claim,
                   standardize = FALSE, standardize.response = FALSE,
                   #foldid = foldid,
                   family = "multinomial",
                   type.multinomial = "ungrouped")

plot(en_claim)
summary(en_claim)[cvm == min(cvm)]
claim_best <- coef(en_claim[[9]], s = 0.01362418)


str(preferable(en_claim), max.level = 1L)

summary(en_claim)[, .SD[cvm == min(cvm)], by = nzero]
plot(en_claim, type = 2)  

multi_6 <- coef(en_claim[[15]], s = 0.015048915)

