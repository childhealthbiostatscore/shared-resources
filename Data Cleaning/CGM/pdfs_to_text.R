library(pdftools)
library(tidyverse)
# Function
# Skip: argument must be numeric and refer to indices of files list
indir = "/Users/timvigers/Dropbox/Mac/Documents/PDFs"
pdf_cgm_data = function(indir,skip = NULL){
  # List files 
  files = list.files(indir,pattern = "*pdf",full.names = T)
  # Summary data frame
  pdf_summary = data.frame()
  # Iterate through files
  for (f in 1:length(files)) {
    # Skip if relevant
    #if (f %in% skip){next}
    # Read PDF into list
    pdf = pdf_data(files[f])
    # First page
    # Page as a dataframe, sort by x and y values
    df = as.data.frame(pdf[[1]])
    df = df %>% arrange(x,y)
    # Get name
    name = paste(gsub('[[:punct:]]','',df$text[which(df$x == 46 & df$y ==3)]),
                 gsub('[[:punct:]]','',df$text[which(df$x > 50 & df$x < 200 & df$y ==3)]))
    # Average and readings
    avg_bg = as.numeric(df$text[which(df$x == 107 & df$y == 186)])
    readings = as.numeric(df$text[which(df$x == 107 & df$y == 206)])
    # Insulin and diet
    basal = as.numeric(df$text[which(df$x == 374 & df$y == 116)])
    bolus = as.numeric(df$text[which(df$x == 469 & df$y == 116)])
    overrides = df$text[which(df$x == 466 & df$y == 161)]
    bolus_day = as.numeric(df$text[which(df$x == 466 & df$y == 171)])
    carbs_day = as.numeric(df$text[which(df$x == 467 & df$y == 227)])
    entries_day = as.numeric(df$text[which(df$x == 467 & df$y == 238)])
    # Add to summary df
    pdf_summary[f,"first_name"] = strsplit(name," ")[[1]][1]
    pdf_summary[f,"last_name"] = strsplit(name," ")[[1]][2]
    pdf_summary[f,"average_bg"] = avg_bg
    pdf_summary[f,"readings_per_day"] = readings
    pdf_summary[f,"basal_units"] = basal
    pdf_summary[f,"bolus_units"] = bolus
    pdf_summary[f,"overrides"] = overrides
    pdf_summary[f,"bolus_per_day"] = bolus_day
    pdf_summary[f,"carbs_per_day"] = carbs_day
    pdf_summary[f,"carb_entries_per_day"] = entries_day
    # Second page
    ## Page as a dataframe, sort by x and y values
    df = as.data.frame(pdf[[2]])
    df = df %>% arrange(x,y)
    # CGM data
    very_high = df$text[which(df$x == 67 & df$y == 95)]
    high = df$text[which(df$x == 67 & df$y == 111)]
    target = df$text[which(df$x == 67 & df$y == 127)]
    low = df$text[which(df$x == 67 & df$y == 143)]
    very_low = df$text[which(df$x == 67 & df$y == 159)]
    avg_sg = as.numeric(df$text[which(df$x == 317 & df$y == 130)])
    sd_sg = as.numeric(sub("mg/dL","",df$text[which(df$x == 548 & df$y == 94)]))
    percent_time = df$text[which(df$x == 317 & df$y == 157)]
    # Add to summary df
    pdf_summary[f,"very_high"] = very_high
    pdf_summary[f,"high"] = high
    pdf_summary[f,"tir"] = target
    pdf_summary[f,"low"] = low
    pdf_summary[f,"very_low"] = very_low
    pdf_summary[f,"avg_sg"] = avg_sg
    pdf_summary[f,"sd_sg"] = sd_sg
    pdf_summary[f,"percent_time"] = percent_time
  }
  return(pdf_summary)
}

pdf_summary = pdf_cgm_data("/Users/timvigers/Dropbox/Mac/Documents/PDFs")
write.csv(pdf_summary,"/Users/timvigers/Dropbox/Mac/Documents/pdf_summary.csv",
          row.names = F,na = "")
