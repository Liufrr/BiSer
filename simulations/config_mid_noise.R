################################################################################
# Simulation Config: Mid Noise Scenarios
################################################################################

source("R/data_generation.R")

configs <- list(

  list(
    name = "norm_mid", generator = generate_norm, n_iter = 150,
    bicrnum = c(20,50,50,30,100), biccnum = c(20,20,20,40,50),
    overr = c(0,0,8,3,20), overc = c(0,0,10,8,15),
    prominence = 0.005, distance = 6,
    gen_args = list(bicmean = c(2,1,4,2,5), bicsd = rep(0.5,5),
                    noisemean = 0, noisesd = 0.2)
  ),

  list(
    name = "norm_mid_1000", generator = generate_norm, n_iter = 20,
    bicrnum = c(100,250,200,150,300), biccnum = c(100,150,100,200,100),
    overr = c(0,10,20,30,20), overc = c(0,20,10,20,15),
    prominence = 0.02, distance = 30,
    gen_args = list(bicmean = c(2,1,4,2,5), bicsd = rep(0.5,5),
                    noisemean = 0, noisesd = 0.2)
  ),

  list(
    name = "norm_exc_mid", generator = generate_norm, n_iter = 150,
    bicrnum = c(20,50,50,30,100), biccnum = c(20,20,20,40,50),
    overr = c(0,0,0,0,0), overc = c(0,0,0,0,0),
    prominence = 0.05, distance = 5,
    gen_args = list(bicmean = c(2,1,4,2,5), bicsd = rep(0.5,5),
                    noisemean = 0, noisesd = 0.2)
  ),

  list(
    name = "norm_exc_mid_1000", generator = generate_norm, n_iter = 20,
    bicrnum = c(100,250,200,150,300), biccnum = c(100,150,100,200,100),
    overr = c(0,0,0,0,0), overc = c(0,0,0,0,0),
    prominence = 0.03, distance = 10,
    gen_args = list(bicmean = c(2,1,4,2,5), bicsd = rep(0.5,5),
                    noisemean = 0, noisesd = 0.2)
  ),

  list(
    name = "NBP_mid", generator = generate_NBP, n_iter = 150,
    bicrnum = c(20,50,50,30,100), biccnum = c(20,20,20,40,50),
    overr = c(0,0,8,3,20), overc = c(0,0,10,8,15),
    prominence = 0.05, distance = 5,
    gen_args = list(bicmean = c(10,18,15,12,16)+70, bicsize = rep(15,5),
                    noisemean = 1, noisesize = 0.1)
  ),

  list(
    name = "NBP_exc_mid", generator = generate_NBP, n_iter = 150,
    bicrnum = c(20,50,50,30,100), biccnum = c(20,20,20,40,50),
    overr = c(0,0,0,0,0), overc = c(0,0,0,0,0),
    prominence = 0.1, distance = 5,
    gen_args = list(bicmean = c(10,18,15,12,16)+70, bicsize = rep(15,5),
                    noisemean = 1, noisesize = 0.1)
  )
)
