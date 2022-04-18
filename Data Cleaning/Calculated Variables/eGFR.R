# Function for calculating "GFR_Schwartz","GFR_FAS","GFR_Zappitelli","GFR_CKIDU25"
# based on Laura's code for Petter Bjornstad's eGFR vs mGFR study

# male and female arguments specify how sex is coded
# Returns a dataframe
egfr_calc = function(age,serum_creatinine,cystatin_c,bun = NULL,height,
                sex,male = "Male",female = "Female",alpha = 0.5){
  # Format input
  sex = as.character(sex)
  sex[sex == male] = "M"
  sex[sex == female] = "F"
  sex[sex != "M" & sex != "F"] = NA
  # BUN to NA
  if(is.null(bun)){
    bun = NA
    length(bun) = length(age)
  }
  # Get qcr
  qcr = floor(age)
  # Younger participants
  qcr[qcr==8] <- 0.46
  qcr[qcr==9] <- 0.49
  qcr[qcr==10] <- 0.51 
  qcr[qcr==11] <- 0.53 
  qcr[qcr==12] <- 0.57
  qcr[qcr==13] <- 0.59
  qcr[qcr==14] <- 0.61
  # Females
  qcr[qcr==15 & sex=="F"] <- 0.64
  qcr[qcr==16 & sex=="F"] <- 0.67
  qcr[qcr==17 & sex=="F"] <- 0.69
  qcr[qcr==18 & sex=="F"] <- 0.69
  qcr[qcr==19 & sex=="F"] <- 0.70
  qcr[qcr>19 & sex=="F"] <- 0.70
  # Males
  qcr[qcr==15 & sex=="M"] <- 0.72
  qcr[qcr==16 & sex=="M"] <- 0.78
  qcr[qcr==17 & sex=="M"] <- 0.82
  qcr[qcr==18 & sex=="M"] <- 0.85
  qcr[qcr==19 & sex=="M"] <- 0.88
  qcr[qcr>19 & sex=="M"] <- 0.90
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
  # eGFR Schwartz
  m = sex
  m[m == "M"] = 1
  m[m == "F"] = 0
  m = suppressWarnings(as.numeric(m))
  eGFR_Schwartz = 39.1*((height/serum_creatinine)^0.516) * ((1.8/cystatin_c)^0.294) * 
    ((30/bun)^0.169) * (1.099^m) * ((height/1.4)^0.188)
  # eGFR bedside Schwartz
  eGFR_bedside_Schwartz <- (41.3*(height/100))/serum_creatinine
  # eGFR CKiD U25 age and sex-dependent sCR
  df = data.frame(cbind(age,sex))
  k_CKiD_scr = apply(df,1,function(r){
    # Sex and age
    s = as.character(r["sex"])
    age = as.numeric(r["age"])
    age_cat = cut(age,c(-Inf,12,15,18,Inf),right = F)
    # Different equations for age and sex
    if(is.na(age_cat) | is.na(s)){
      k = NA
    } else if (age_cat == "[-Inf,12)"){
      if (s == "M"){
        k = 39*(1.008^(age-12))
      } else if (s == "F"){
        k = 36.1*(1.008^(age-12))
      }
    } else if (age_cat == "[12,15)") {
      if (s == "M"){
        k = 39*(1.045^(age-12))
      } else if (s == "F"){
        k = 36.1*(1.023^(age-12))
      }
    } else if (age_cat == "[15,18)") {
      if (s == "M"){
        k = 39*(1.045^(age-12))
      } else if (s == "F"){
        k = 36.1*(1.023^(age-12))
      }
    } else if (age_cat == "[18, Inf)") {
      if (s == "M"){
        k = 50.8
      } else if (s == "F"){
        k = 41.4
      }
    }
    return(k)
  })
  eGFR_CKiD_scr <- k_CKiD_scr * ((height/100)/serum_creatinine)
  # eGFR CKiD U25 age and sex-dependent cystatin-C
  k_CKiD_cysc = apply(df,1,function(r){
    # Sex and age
    s = as.character(r["sex"])
    age = as.numeric(r["age"])
    age_cat = cut(age,c(-Inf,12,15,18,Inf),right = F)
    # Different equations for age and sex
    if(is.na(age_cat) | is.na(s)){
      k = NA
    } else if(age_cat == "[-Inf,12)"){
      if (s == "M"){
        k = 87.2*(1.011^(age-15))
      } else if (s == "F"){
        k = 79.9*(1.004^(age-12))
      }
    } else if (age_cat == "[12,15)") {
      if (s == "M"){
        k = 87.2*(1.011^(age-15))
      } else if (s == "F"){
        k = 79.9*(0.974^(age-12))
      }
    } else if (age_cat == "[15,18)") {
      if (s == "M"){
        k = 87.2*(0.960^(age-15))
      } else if (s == "F"){
        k = 79.9*(0.974^(age-12))
      }
    } else if (age_cat == "[18, Inf)") {
      if (s == "M"){
        k = 77.1
      } else if (s == "F"){
        k = 68.3
      }
    }
    return(k)
  })
  eGFR_CKiD_cysc <- k_CKiD_cysc * (1/cystatin_c)
  # Return dataframe
  return(data.frame(list("eGFR_fas_cr" = eGFR_fas_cr,"eGFR_fas_cr_cysc" = eGFR_fas_cr_cysc,
                         "eGFR_Zap" = eGFR_Zap,"eGFR_Schwartz" = eGFR_Schwartz,
                         "eGFR_bedside_Schwartz" = eGFR_bedside_Schwartz,
                         "eGFR_CKiD_scr" = eGFR_CKiD_scr,"eGFR_CKiD_cysc" = eGFR_CKiD_cysc)))
}
