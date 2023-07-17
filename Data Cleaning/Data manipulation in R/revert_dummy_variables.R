# This function takes a dataframe, column name, and prefix as 
# arguments. This assumes that the dummy variables are named following the 
# format "prefix{value}" (e.g. "race1", "race2", etc.). There is an option to 
# keep the dummy variables if desired. The levels argument links numeric values
# factor names, but is not required. If provided, the levels need to be listed 
# in numeric order, not the order that the columns appear in the dataframe.

revert_dummies = function(df,new_col_name,prefix,levels,keep_cols = F){
  reg = paste0(prefix,"\\d")
  cols = sort(colnames(df)[grep(reg,colnames(df))])
  df[,new_col_name] = apply(df[,cols],1,function(r){
    w = which(r==1)
    if(!is.null(levels)){
      w = levels[w]
    }
    return(paste(w,collapse = ", "))
  })
  df[,new_col_name] = factor(df[,new_col_name])
  if(keep_cols==F){
    df = df[,-match(cols,colnames(df))]
  }
  return(df)
}
levels = c("Health Insurance Policy (e.g. Private Insurance)",
                   "Medicare",
                   "Medicaid",
                   "State special needs program, e.g., BCMH, CCS, CRS, GHPP, etc.",
                   "TriCare or other military health plan",
                   "Indian Health Service",
                   "Other")
