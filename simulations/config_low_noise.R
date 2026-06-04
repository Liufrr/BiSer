################################################################################
# Simulation Config: Low Noise Scenarios
#
# Normal, Normal-exclusive, NBP, and Poisson distributions.
# Two matrix sizes: small (250x150) and large (1000x650).
################################################################################

source("R/data_generation.R")

configs <- list(

  # --- Normal, overlapping, small ---
  list(
    name = "norm_low", generator = generate_norm, n_iter = 150,
    bicrnum = c(20, 50, 50, 30, 100), biccnum = c(20, 20, 20, 40, 50),
    overr = c(0, 0, 8, 3, 20), overc = c(0, 0, 10, 8, 15),
    prominence = 0.002, distance = 5,
    gen_args = list(bicmean = c(2,1,4,2,5)+10, bicsd = rep(0.5,5),
                    noisemean = 0, noisesd = 0.2)
  ),

  # --- Normal, overlapping, large ---
  list(
    name = "norm_low_1000", generator = generate_norm, n_iter = 60,
    bicrnum = c(100, 250, 200, 150, 300), biccnum = c(100, 150, 100, 200, 100),
    overr = c(0, 10, 20, 30, 20), overc = c(0, 20, 10, 20, 15),
    prominence = 0.02, distance = 10,
    gen_args = list(bicmean = c(2,1,4,2,5)+10, bicsd = rep(0.5,5),
                    noisemean = 0, noisesd = 0.2)
  ),

  # --- Normal, exclusive, small ---
  list(
    name = "norm_exc_low", generator = generate_norm, n_iter = 150,
    bicrnum = c(20, 50, 50, 30, 100), biccnum = c(20, 20, 20, 40, 50),
    overr = c(0,0,0,0,0), overc = c(0,0,0,0,0),
    prominence = 0.01, distance = 10,
    gen_args = list(bicmean = c(2,1,4,2,5)+10, bicsd = rep(0.5,5),
                    noisemean = 0, noisesd = 0.2)
  ),

  # --- Normal, exclusive, large ---
  list(
    name = "norm_exc_low_1000", generator = generate_norm, n_iter = 60,
    bicrnum = c(100, 250, 200, 150, 300), biccnum = c(100, 150, 100, 200, 100),
    overr = c(0,0,0,0,0), overc = c(0,0,0,0,0),
    prominence = 0.1, distance = 5,
    gen_args = list(bicmean = c(2,1,4,2,5)+10, bicsd = rep(0.5,5),
                    noisemean = 0, noisesd = 0.2)
  ),

  # --- NBP, overlapping, small ---
  list(
    name = "NBP_low", generator = generate_NBP, n_iter = 150,
    bicrnum = c(20, 50, 50, 30, 100), biccnum = c(20, 20, 20, 40, 50),
    overr = c(0, 0, 8, 3, 20), overc = c(0, 0, 10, 8, 15),
    prominence = 0.07, distance = 5,
    gen_args = list(bicmean = c(10,18,15,12,16)+70, bicsize = rep(15,5),
                    noisemean = 0, noisesize = 1)
  ),

  # --- NBP, exclusive, small ---
  list(
    name = "NBP_exc_low", generator = generate_NBP, n_iter = 150,
    bicrnum = c(20, 50, 50, 30, 100), biccnum = c(20, 20, 20, 40, 50),
    overr = c(0,0,0,0,0), overc = c(0,0,0,0,0),
    prominence = 0.1, distance = 10,
    gen_args = list(bicmean = c(10,18,15,12,16)+70, bicsize = rep(15,5),
                    noisemean = 0, noisesize = 1)
  )
)
