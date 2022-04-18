# Function for calculating "GFR_Schwartz","GFR_FAS","GFR_Zappitelli","GFR_CKIDU25"
# based on Laura's code for Petter Bjornstad's eGFR vs mGFR study

# male and female arguments specify how sex is coded
calculate_egfr = function(age,sex,male = "Male",female = "Female"){
  age = floor(age)
  sex[sex == male] = "m"
  sex[sex == female] = "f"
}


# calculate QCR values for FAS equations
alldata$age <- floor(alldata$age)
alldata$qcr[alldata$age==8] <- 0.46
alldata$qcr[alldata$age==9] <- 0.49
alldata$qcr[alldata$age==10] <- 0.51 
alldata$qcr[alldata$age==11] <- 0.53 
alldata$qcr[alldata$age==12] <- 0.57
alldata$qcr[alldata$age==13] <- 0.59
alldata$qcr[alldata$age==14] <- 0.61
# females
alldata$qcr[alldata$age==15 & alldata$sex_MF=="F"] <- 0.64
# males
alldata$qcr[alldata$age==15 & alldata$sex_MF=="M"] <- 0.72
# females
alldata$qcr[alldata$age==16 & alldata$sex_MF=="F"] <- 0.67
# males
alldata$qcr[alldata$age==16 & alldata$sex_MF=="M"] <- 0.78
# females
alldata$qcr[alldata$age==17 & alldata$sex_MF=="F"] <- 0.69
# males
alldata$qcr[alldata$age==17 & alldata$sex_MF=="M"] <- 0.82
# females
alldata$qcr[alldata$age==18 & alldata$sex_MF=="F"] <- 0.69
# males
alldata$qcr[alldata$age==18 & alldata$sex_MF=="M"] <- 0.85
# females
alldata$qcr[alldata$age==19 & alldata$sex_MF=="F"] <- 0.70
# males
alldata$qcr[alldata$age==19 & alldata$sex_MF=="M"] <- 0.88
# females
alldata$qcr[alldata$age>19 & alldata$sex_MF=="F"] <- 0.70
# males
alldata$qcr[alldata$age>19 & alldata$sex_MF=="M"] <- 0.90

# eGFR FAS creatinine
alldata$eGFR.fas_cr <-107.3/(alldata$serum_creatinine/alldata$qcr)

# eGFR FAS combined creatinine and cystatin-C
alldata$f1 <- alldata$serum_creatinine/alldata$qcr
alldata$f2 <- 1-0.5
alldata$f3 <- alldata$cystatin_c/0.82
alldata$eGFR.fas_cr_cysc <- 107.3 / ((0.5*alldata$f1) + (alldata$f2*alldata$f3))

# eGFR Zapatelli
alldata$eGFR.Zap <- (507.76*exp(0.003*(alldata$clamp_height)))/((alldata$cystatin_c^0.635)*((alldata$serum_creatinine*88.4)^0.547))

# eGFR bedside Schwartz
alldata$eGFR.bedside_Schwartz <- (41.3*(alldata$clamp_height/100))/alldata$serum_creatinine

# eGFR CKiD U25 age and sex-dependent sCR
alldata$k_CKiD_scr = apply(alldata,1,function(r){
  # Sex and age
  sex = as.character(r["sex_MF"])
  age = as.numeric(r["age"])
  age_cat = cut(age,c(-Inf,12,15,18,Inf),right = F)
  # Different equations for age and sex
  if(age_cat == "[-Inf,12)"){
    if (sex == "M"){
      k = 39*(1.008^(age-12))
    } else if (sex == "F"){
      k = 36.1*(1.008^(age-12))
    }
  } else if (age_cat == "[12,15)") {
    if (sex == "M"){
      k = 39*(1.045^(age-12))
    } else if (sex == "F"){
      k = 36.1*(1.023^(age-12))
    }
  } else if (age_cat == "[15,18)") {
    if (sex == "M"){
      k = 39*(1.045^(age-12))
    } else if (sex == "F"){
      k = 36.1*(1.023^(age-12))
    }
  } else if (age_cat == "[18, Inf)") {
    if (sex == "M"){
      k = 50.8
    } else if (sex == "F"){
      k = 41.4
    }
  }
  return(k)
})
alldata$eGFR.CKiD_scr <- alldata$k_CKiD_scr * ((alldata$clamp_height/100)/alldata$serum_creatinine)

# eGFR CKiD U25 age and sex-dependent cystatin-C
alldata$k_CKiD_cysc = apply(alldata,1,function(r){
  # Sex and age
  sex = as.character(r["sex_MF"])
  age = as.numeric(r["age"])
  age_cat = cut(age,c(-Inf,12,15,18,Inf),right = F)
  # Different equations for age and sex
  if(age_cat == "[-Inf,12)"){
    if (sex == "M"){
      k = 87.2*(1.011^(age-15))
    } else if (sex == "F"){
      k = 79.9*(1.004^(age-12))
    }
  } else if (age_cat == "[12,15)") {
    if (sex == "M"){
      k = 87.2*(1.011^(age-15))
    } else if (sex == "F"){
      k = 79.9*(0.974^(age-12))
    }
  } else if (age_cat == "[15,18)") {
    if (sex == "M"){
      k = 87.2*(0.960^(age-15))
    } else if (sex == "F"){
      k = 79.9*(0.974^(age-12))
    }
  } else if (age_cat == "[18, Inf)") {
    if (sex == "M"){
      k = 77.1
    } else if (sex == "F"){
      k = 68.3
    }
  }
  return(k)
})
alldata$eGFR.CKiD_cysc <- alldata$k_CKiD_cysc * (1/alldata$cystatin_c)