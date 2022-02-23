#FINAL CLUSTERING:

library(cluster)
library(factoextra)
library(Rtsne)
library(Hmisc)
library(dplyr)
library(clValid)
library(fpc)
library(NbClust)

##variable names for cluster:
vars_num <- c("MRN","number_elig_months","death_1st_yr",
              #resource use
              "Hospitalization","ED","Pulmonary","Rehabilitation","Primarycare","Therapy",
              #medications
              "meds_gastrointestinal_num","meds_pulmonary_num","meds_anti_infectives_num","meds_pain_num",
              "meds_endocrine_num","meds_cardiac_num","meds_diuretics_num",
              "meds_steroids_num","meds_infl_num","meds_hema_num","meds_n.muscle_num","meds_n.psych_num",
              "meds_n.epi_num","meds_n.stim_num","meds_n.oth_num","meds_uro_num",
              #CPT
              #"speech_num","pt_ot_num", #removed on 11/10/2020
              #DME
              "dme_feeding_num","dme_mobility_num","dme_respiratory_num","dme_infusion_num",
              "dme_therapy_num","dme_foley_num", "dme_tracheostomy_num","dme_hear_comm_num" #hearing/communication added on 10/23/2020
)
#subset to only relavent variables for cluster:
dat.sub<-dat.apcd[vars_num]
N<-ncol(dat.sub)


#adjust each variable by number of eligible months:
divide<-function(x){
  #x<-dat$speech_yn
  x<-as.numeric(as.character(x))
  x<-x/dat.apcd$number_elig_months
  return(as.numeric(x))
}
standardize<-function(x){
  x<-x-mean(x)/sd(x)
  return(x)
}


dat.sub[,c(4:N)]<-lapply(dat.sub[,c(4:N)],divide)
dat.use<-dat.sub
dat.use[,c(4:N)]<-lapply(dat.sub[,c(4:N)],standardize)
dat.use<-dat.use[,-c(1,2)]
N<-ncol(dat.use)

#####CARTER'S CODE#######
dissmat <- cluster::daisy(dat.use,metric='gower', 
                          
                          #  type=list(ordratio=c(1:21) 
                          type=list(asymm=c(1)
                          ))
sil_width <- c(NA)

for(i in 2:10){
  
  pam_fit <- cluster::pam(dissmat,
                          diss = TRUE,
                          k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}


gowerClust <- cluster::pam(x=dissmat, diss=T, k = 3)
#table(gowerClust$clustering)
dat.apcd["clustnum"] <- gowerClust$clustering


