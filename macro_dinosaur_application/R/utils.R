# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project: Early SARS-CoV-2 Europe Phylodynamics
#        V\ Y /V    Util functions for analysis workflow
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------



read_trace <- function(traceFile, burninFrac){
  df <- read_table(traceFile, comment = "#")
  
  if (burninFrac > 0) {
    n <- dim(df)[1]
    df <- df[-(1:ceiling(burninFrac * n)), ]
  }
  
  return(df)
}
