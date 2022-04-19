# Function for calculating various hemodynamics for Petter Bjornstad's 
# studies

# SAS code from Petter:
# data SUA.r01;
# set SUA.r01;
# 
# 
# 
# TotProtgdL = total_protein;
# HCt = mean (hematocrit_minus_5, hematocrit_90, hematocrit_120);
# MAP_clamp = clamp_map;
# HCt_clamp = HCt / 100;
# 
# 
# 
# * Absolute;
# 
# GFRabsmLs= GFR/60;
# ERPFabsmLs= RPF/60;
# FFabs= GFR/RPF;
# RBFabs= RPF /(1-HCt_clamp);
# RBFGomabs= (ERPFabsmLs/(1-HCt_clamp));
# 
# 
# 
# run;
# 
# 
# 
# data SUA.r01;
# set SUA.r01;
# 
# 
# 
# RVRAbs= clamp_map/RBFabs ;
# 
# if group = 2 then Kfg=0.1012;
# if group = 3 then Kfg=0.1733;
# if group = 4 then Kfg=0.1733; *** for non-DM if DM=0 then Kfg=0.1733;*;
# 
# 
# 
# 
# deltaPFabs= GFRabsmLs/Kfg;
# 
# Cm_abs= (TotProtgdL/FFabs)*log(1/(1-FFabs));
# Pg_abs= 5*(Cm_abs-2);
# Pglo_abs= Pg_abs+deltaPFabs+10;
# Ra_abs= ((MAP_clamp-Pglo_abs)/RBFGomabs)*1328;
# Re_abs= (GFRabsmLs/(Kfg*(RBFGomabs-GFRabsmLs)))*1328;
# 
# 
# 
# run;

# In general studies with have either a minus 10 or minus 5 value for hematocrit,
# but double check that it's okay to include both in the average if available.
hemodynamics_calc = function(total_protein,gfr,rpf,map,
                             hematocrit_minus_10,hematocrit_minus_5,
                             hematocrit_90,hematocrit_120,group){
  # Format input
  TotProtgdL = total_protein
  HCt = rowMeans(cbind(hematocrit_minus_10,hematocrit_minus_5,
                       hematocrit_90,hematocrit_120),na.rm = T)
  MAP_clamp = map;
  HCt_clamp = HCt / 100;
  # Calculate variables
  GFRabsmLs= gfr/60;
  ERPFabsmLs= rpf/60;
  FFabs= gfr/rpf;
  RBFabs= rpf /(1-HCt_clamp);
  RBFGomabs= (ERPFabsmLs/(1-HCt_clamp));
  RVRAbs = map/RBFabs ;

  
  if group = 2 then Kfg=0.1012;
  if group = 3 then Kfg=0.1733;
  if group = 4 then Kfg=0.1733; *** for non-DM if DM=0 then Kfg=0.1733;*;




  deltaPFabs= GFRabsmLs/Kfg;

  Cm_abs = (TotProtgdL/FFabs)*log(1/(1-FFabs));
  Pg_abs = 5*(Cm_abs-2);
  Pglo_abs = Pg_abs+deltaPFabs+10;
  Ra_abs = ((MAP_clamp-Pglo_abs)/RBFGomabs)*1328;
  Re_abs = (GFRabsmLs/(Kfg*(RBFGomabs-GFRabsmLs)))*1328;
  
  
  # Return dataframe
  return(data.frame(list()))
}
