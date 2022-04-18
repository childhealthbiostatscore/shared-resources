# Function for calculating "GFR_Schwartz","GFR_FAS","GFR_Zappitelli","GFR_CKIDU25"
# based on Laura's code for Petter Bjornstad's eGFR vs mGFR study

# male and female arguments specify how sex is coded
# Returns a dataframe
calculate_egfr = function(age,serum_creatinine,cystatin_c,height,
                          sex,male = "Male",female = "Female",alpha = 0.5){
  # Format input
  age = floor(age)
  sex = as.character(sex)
  sex[sex == male] = "M"
  sex[sex == female] = "F"
  sex[sex != "M" & sex != "F"] = NA
  # Get qcr
  qcr = age
  # Younger participants
  qcr[age==8] <- 0.46
  qcr[age==9] <- 0.49
  qcr[age==10] <- 0.51 
  qcr[age==11] <- 0.53 
  qcr[age==12] <- 0.57
  qcr[age==13] <- 0.59
  qcr[age==14] <- 0.61
  # Females
  qcr[age==15 & sex=="F"] <- 0.64
  qcr[age==16 & sex=="F"] <- 0.67
  qcr[age==17 & sex=="F"] <- 0.69
  qcr[age==18 & sex=="F"] <- 0.69
  qcr[age==19 & sex=="F"] <- 0.70
  qcr[age>19 & sex=="F"] <- 0.70
  # Males
  qcr[age==15 & sex=="M"] <- 0.72
  qcr[age==16 & sex=="M"] <- 0.78
  qcr[age==17 & sex=="M"] <- 0.82
  qcr[age==18 & sex=="M"] <- 0.85
  qcr[age==19 & sex=="M"] <- 0.88
  qcr[age>19 & sex=="M"] <- 0.90
  # Calculate final metrics
  eGFR_fas_cr = 107.3/(serum_creatinine/qcr)
  # eGFR FAS combined creatinine and cystatin-C
  f1 <- serum_creatinine/qcr
  f2 <- 1-0.5
  f3 <- cystatin_c/0.82
  eGFR_fas_cr_cysc <- 107.3 / ((0.5*f1) + (f2*f3))
  # eGFR Zapatelli
  eGFR_Zap <- (507.76*exp(0.003*(height)))/
    ((cystatin_c^0.635)*((serum_creatinine*88.4)^0.547))
  # eGFR bedside Schwartz
  eGFR_bedside_Schwartz <- (41.3*(height/100))/serum_creatinine
  # Return dataframe
  return(data.frame(list("eGFR_fas_cr" = eGFR_fas_cr,"eGFR_fas_cr_cysc" = eGFR_fas_cr_cysc,
                         "eGFR_Zap" = eGFR_Zap)))
}



# females

# males








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