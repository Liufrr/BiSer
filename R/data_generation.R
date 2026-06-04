################################################################################
# BiSer Benchmarking: Simulation Data Generation
#
# Functions to generate block-structured matrices under different distributional
# assumptions: Normal, Negative Binomial (NB), and Poisson.
# Both exclusive (non-overlapping) and overlapping bicluster configurations
# are supported via the overr/overc parameters.
################################################################################

# ==============================================================================
# Normal distribution block matrix
# ==============================================================================
generate_norm <- function(bicrnum = c(20, 50, 50, 30, 100),
                          biccnum = c(20, 20, 20, 40, 50),
                          overr   = c(0, 0, 5, 3, 15),
                          overc   = c(0, 2, 4, 0, 8),
                          bicmean = c(2, 1, 4, 2, 5),
                          bicsd   = rep(0.5, 5),
                          noisemean = 0,
                          noisesd   = 0.2,
                          seed = NA) {
  if (!is.na(seed)) set.seed(seed)
  m <- sum(bicrnum)
  p <- sum(biccnum)
  mat <- matrix(rnorm(m * p, mean = noisemean, sd = noisesd), nrow = m)

  bicnum  <- length(bicrnum)
  centerr <- c(1, cumsum(bicrnum[-bicnum]) + 1)
  centerc <- c(1, cumsum(biccnum[-bicnum]) + 1)

  for (i in 1:bicnum) {
    row_start <- max(1, centerr[i] - overr[i])
    row_end   <- min(m, centerr[i] + bicrnum[i] + overr[i] - 1)
    col_start <- max(1, centerc[i] - overc[i])
    col_end   <- min(p, centerc[i] + biccnum[i] + overc[i] - 1)

    block_rows <- row_start:row_end
    block_cols <- col_start:col_end

    block_signal <- matrix(rnorm(length(block_rows) * length(block_cols),
                                 mean = bicmean[i], sd = bicsd[i]),
                           nrow = length(block_rows))
    mat[block_rows, block_cols] <- mat[block_rows, block_cols] + block_signal
  }

  return(mat)
}

# ==============================================================================
# Negative Binomial distribution block matrix
# ==============================================================================
generate_NBP <- function(bicrnum = c(20, 50, 50, 30, 100),
                         biccnum = c(20, 20, 20, 40, 50),
                         overr   = c(0, 0, 5, 3, 15),
                         overc   = c(0, 2, 4, 0, 8),
                         bicmean = c(2, 1, 4, 2, 5),
                         bicsize = rep(0.5, 5),
                         noisemean = 0,
                         noisesize = 0.2,
                         seed = NA) {
  if (!is.na(seed)) set.seed(seed)
  m <- sum(bicrnum)
  p <- sum(biccnum)
  mat <- matrix(rnbinom(m * p, mu = noisemean, size = noisesize), nrow = m)

  bicnum  <- length(bicrnum)
  centerr <- c(1, cumsum(bicrnum[-bicnum]) + 1)
  centerc <- c(1, cumsum(biccnum[-bicnum]) + 1)

  # Clear block regions before adding signal
  for (i in 1:bicnum) {
    row_start <- max(1, centerr[i] - overr[i])
    row_end   <- min(m, centerr[i] + bicrnum[i] + overr[i] - 1)
    col_start <- max(1, centerc[i] - overc[i])
    col_end   <- min(p, centerc[i] + biccnum[i] + overc[i] - 1)
    block_rows <- row_start:row_end
    block_cols <- col_start:col_end
    mat[block_rows, block_cols] <- matrix(0, nrow = length(block_rows),
                                          ncol = length(block_cols))
  }

  # Add block signal
  for (i in 1:bicnum) {
    row_start <- max(1, centerr[i] - overr[i])
    row_end   <- min(m, centerr[i] + bicrnum[i] + overr[i] - 1)
    col_start <- max(1, centerc[i] - overc[i])
    col_end   <- min(p, centerc[i] + biccnum[i] + overc[i] - 1)
    block_rows <- row_start:row_end
    block_cols <- col_start:col_end

    block_signal <- matrix(rnbinom(length(block_rows) * length(block_cols),
                                   mu = bicmean[i], size = bicsize[i]),
                           nrow = length(block_rows))
    mat[block_rows, block_cols] <- mat[block_rows, block_cols] + block_signal
  }

  return(mat)
}

# ==============================================================================
# Poisson distribution block matrix
# ==============================================================================
generate_poisson <- function(bicrnum = c(20, 50, 50, 30, 100),
                             biccnum = c(20, 20, 20, 40, 50),
                             overr   = c(0, 0, 5, 3, 15),
                             overc   = c(0, 2, 4, 0, 8),
                             bicmean = c(2, 1, 4, 2, 5),
                             noisemean = 0.5,
                             seed = NA) {
  if (!is.na(seed)) set.seed(seed)
  m <- sum(bicrnum)
  p <- sum(biccnum)
  mat <- matrix(rpois(m * p, lambda = noisemean), nrow = m)

  bicnum  <- length(bicrnum)
  centerr <- c(1, cumsum(bicrnum[-bicnum]) + 1)
  centerc <- c(1, cumsum(biccnum[-bicnum]) + 1)

  for (i in 1:bicnum) {
    row_start <- max(1, centerr[i] - overr[i])
    row_end   <- min(m, centerr[i] + bicrnum[i] + overr[i] - 1)
    col_start <- max(1, centerc[i] - overc[i])
    col_end   <- min(p, centerc[i] + biccnum[i] + overc[i] - 1)
    block_rows <- row_start:row_end
    block_cols <- col_start:col_end

    block_signal <- matrix(rpois(length(block_rows) * length(block_cols),
                                 lambda = bicmean[i]),
                           nrow = length(block_rows))
    mat[block_rows, block_cols] <- mat[block_rows, block_cols] + block_signal
  }

  return(mat)
}

# ==============================================================================
# Simple block matrix for scalability testing
# ==============================================================================
generate_block_matrix <- function(m, p, n_blocks = 3, signal = 2, seed = 42) {
  set.seed(seed)
  mat <- matrix(rnorm(m * p), nrow = m, ncol = p)
  row_blocks <- cut(1:m, breaks = n_blocks, labels = FALSE)
  col_blocks <- cut(1:p, breaks = n_blocks, labels = FALSE)
  for (b in 1:n_blocks) {
    ri <- which(row_blocks == b)
    ci <- which(col_blocks == b)
    mat[ri, ci] <- mat[ri, ci] + signal
  }
  mat <- mat[sample(m), sample(p)]
  rownames(mat) <- paste0("R", 1:m)
  colnames(mat) <- paste0("C", 1:p)
  return(mat)
}

# ==============================================================================
# Ground truth label generation
# ==============================================================================

#' Label non-zero overlap regions distinctly from core blocks
label_nonzero_regions <- function(x, x0, start_label = 5) {
  rle_x  <- rle(x != 0)
  labels <- rep(0, length(rle_x$lengths))

  label <- start_label
  for (i in seq_along(rle_x$lengths)) {
    if (rle_x$values[i]) {
      labels[i] <- label
      label <- label + 1
    }
  }

  tag_vec <- inverse.rle(list(lengths = rle_x$lengths, values = labels))
  x0[x != 0] <- tag_vec[x != 0]
  return(x0)
}

#' Generate ground truth labels for overlapping bicluster configurations
generate_true_labels <- function(bicrnum, biccnum, overr, overc, indr, indc) {
  bicnum  <- length(bicrnum)
  m <- sum(bicrnum)
  p <- sum(biccnum)

  centerr <- c(1, cumsum(bicrnum[-bicnum]) + 1)
  centerc <- c(1, cumsum(biccnum[-bicnum]) + 1)

  tmpr1 <- tmpr2 <- tmpc1 <- tmpc2 <- numeric()
  for (i in 1:bicnum) {
    tmpr1[centerr[i]:(centerr[i] + bicrnum[i] - 1)] <- i
    tmpr2[max(1, centerr[i] - overr[i]):min(m, centerr[i] + bicrnum[i] + overr[i] - 1)] <- i
    tmpc1[centerc[i]:(centerc[i] + biccnum[i] - 1)] <- i
    tmpc2[max(1, centerc[i] - overc[i]):min(p, centerc[i] + biccnum[i] + overc[i] - 1)] <- i
  }

  truelabel <- list(
    row = label_nonzero_regions(tmpr1 - tmpr2, tmpr1, bicnum + 1)[indr],
    col = label_nonzero_regions(tmpc1 - tmpc2, tmpc1, bicnum + 1)[indc]
  )

  return(truelabel)
}
