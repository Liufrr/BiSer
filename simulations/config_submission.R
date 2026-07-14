################################################################################
# BiSer0608 submission simulation configuration
#
# Full factorial design reported in the submitted manuscript:
#   3 distributions x 2 structures x 2 noise levels x 2 matrix sizes
# Each of the 24 conditions is evaluated in 100 independent replicates.
################################################################################

source("R/data_generation.R")

submission_boundary <- list(
  window = 10L,
  smooth = "gaussian",
  sigma = 3,
  prominence = 0.01,
  distance = 5L
)

submission_sizes <- list(
  small = list(
    bicrnum = c(20, 50, 50, 30, 100),
    biccnum = c(20, 20, 20, 40, 50),
    overr = c(0, 0, 8, 3, 20),
    overc = c(0, 0, 10, 20, 15)
  ),
  large = list(
    bicrnum = c(100, 250, 200, 150, 300),
    biccnum = c(100, 150, 100, 200, 100),
    overr = c(0, 10, 20, 30, 20),
    overc = c(0, 20, 10, 20, 15)
  )
)

submission_distributions <- list(
  gaussian = list(
    prefix = "norm",
    generator = generate_norm,
    low = list(
      bicmean = c(12, 11, 14, 12, 15),
      bicsd = rep(0.5, 5), noisemean = 0, noisesd = 0.2
    ),
    high = list(
      bicmean = c(1, 1.5, 1, 1.5, 1.5),
      bicsd = rep(0.3, 5), noisemean = 0.2, noisesd = 0.3
    )
  ),
  poisson = list(
    prefix = "poisson",
    generator = generate_poisson,
    low = list(bicmean = c(11, 12, 14, 13, 12), noisemean = 0),
    high = list(bicmean = c(10, 12, 14, 9, 11), noisemean = 1)
  ),
  negative_binomial = list(
    prefix = "NBP",
    generator = generate_NBP,
    low = list(
      bicmean = c(80, 88, 85, 82, 86),
      bicsize = rep(15, 5), noisemean = 0, noisesize = 1
    ),
    high = list(
      bicmean = c(40, 48, 45, 42, 46),
      bicsize = rep(15, 5), noisemean = 1, noisesize = 0.1
    )
  )
)

make_submission_config <- function(distribution, structure, noise_level, scale,
                                   condition_id) {
  dist_spec <- submission_distributions[[distribution]]
  size_spec <- submission_sizes[[scale]]
  is_exclusive <- identical(structure, "exclusive")

  name <- paste0(
    dist_spec$prefix,
    if (is_exclusive) "_exc" else "",
    "_", noise_level,
    if (identical(scale, "large")) "_1000" else ""
  )

  c(
    list(
      name = name,
      distribution = distribution,
      structure = structure,
      noise_level = noise_level,
      scale = scale,
      generator = dist_spec$generator,
      n_iter = 100L,
      n_starts = 100L,
      nmf_nrun = 10L,
      candidate_k = 2:10,
      seed = 20260608L + condition_id * 1000L,
      bicrnum = size_spec$bicrnum,
      biccnum = size_spec$biccnum,
      overr = if (is_exclusive) rep(0, 5) else size_spec$overr,
      overc = if (is_exclusive) rep(0, 5) else size_spec$overc,
      gen_args = dist_spec[[noise_level]]
    ),
    submission_boundary
  )
}

design <- expand.grid(
  distribution = names(submission_distributions),
  structure = c("exclusive", "overlapping"),
  noise_level = c("low", "high"),
  scale = c("small", "large"),
  stringsAsFactors = FALSE
)

configs <- lapply(seq_len(nrow(design)), function(i) {
  make_submission_config(
    distribution = design$distribution[i],
    structure = design$structure[i],
    noise_level = design$noise_level[i],
    scale = design$scale[i],
    condition_id = i
  )
})

names(configs) <- vapply(configs, `[[`, character(1), "name")

stopifnot(
  length(configs) == 24L,
  length(unique(names(configs))) == 24L,
  all(vapply(configs, `[[`, integer(1), "n_iter") == 100L)
)
