# This script provides exploratory data analysis (EDA) functions that work
# with longitudinal data.


# Basic EDA. By default it assumes that data is in long format
# and requires the user to specify the name of the column that indexes time.

easy_eda <- function(data, output_dir = tempdir(), report_name = "report.html",
                     time_variable = NULL, long_format = T, groups = NULL) {
  # Check longitudinal formatting
  if (is.null(time_variable)) {
    stop("Please provide a value for the time_variable argument or set long_format = FALSE")
  }
  # Summarize variable types
  
}
