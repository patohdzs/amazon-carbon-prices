# Start timer
tictoc::tic(msg = "_masterfile.R script", log = TRUE)

source("rsrc/cleaning/_masterfile.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

print("Data cleaning is done")

source("rsrc/processing/_masterfile.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

print("Data processing is done")

source("rsrc/calibration/_masterfile.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

print("Calibration is done")

# End timer
tictoc::toc(log = TRUE)
