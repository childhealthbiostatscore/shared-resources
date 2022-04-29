# Function for calculating various hemodynamics for Petter Bjornstad's 
# studies. 
# Groups are T1D, T2D, Obese Controls, Lean Lontrols.
hemodynamics_calc = function(total_protein,gfr,rpf,map,
                             hematocrit_minus_10,hematocrit_minus_5,
                             hematocrit_90,hematocrit_120,group){
  # Format input
  TotProtgdL = total_protein
  HCt = rowMeans(cbind(hematocrit_minus_10,hematocrit_minus_5,
                       hematocrit_90,hematocrit_120),na.rm = T)
  HCt[is.nan(HCt)] = NA
  MAP_clamp = map
  HCt_clamp = HCt / 100
  # Calculate variables
  GFRabsmLs = gfr/60
  ERPFabsmLs = rpf/60
  FFabs = gfr/rpf
  RBFabs = rpf /(1-HCt_clamp)
  RBFGomabs = (ERPFabsmLs/(1-HCt_clamp))
  RVRAbs = map/RBFabs 
  # Kfg by group
  Kfg = NA
  length(Kfg) = length(group)
  Kfg[group %in% c("T1D","T2D")] = 0.1012
  Kfg[group %in% c("Obese Control","Lean Control")] = 0.1733
  # More calculations
  deltaPFabs = GFRabsmLs/Kfg;
  Cm_abs = (TotProtgdL/FFabs)*log(1/(1-FFabs))
  Pg_abs = 5*(Cm_abs-2)
  Pglo_abs = Pg_abs+deltaPFabs+10
  Ra_abs = ((MAP_clamp-Pglo_abs)/RBFGomabs)*1328
  Re_abs = (GFRabsmLs/(Kfg*(RBFGomabs-GFRabsmLs)))*1328
  # Return dataframe
  return(data.frame(list("GFRabsmLs" = GFRabsmLs,
                         "ERPFabsmLs" = ERPFabsmLs,"FFabs" = FFabs,
                         "RBFabs" = RBFabs,"RBFGomabs" = RBFGomabs,
                         "RVRAbs" = RVRAbs,"deltaPFabs" = deltaPFabs,
                         "Cm_abs" = Cm_abs,"Pg_abs" = Pg_abs,
                         "Pglo_abs" = Pglo_abs,"Ra_abs" = Ra_abs,
                         "Re_abs" = Re_abs)))
}
