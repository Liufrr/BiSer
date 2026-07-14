################################################################################
# BiSer0608 low-noise subset (12 of the 24 submitted simulation conditions)
################################################################################

source("simulations/config_submission.R")
configs <- Filter(function(x) identical(x$noise_level, "low"), configs)
