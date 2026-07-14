################################################################################
# BiSer0608: evaluation metrics and discretization-penalty calculation
################################################################################

library(mclust)
library(aricode)
library(clue)

submission_metric_names <- c(
  "NMI_row", "NMI_col", "NMI",
  "purity_row", "purity_col", "purity",
  "ARI_row", "ARI_col", "ARI",
  "Accuracy_row", "Accuracy_col", "Accuracy",
  "precision_row", "precision_col", "precision",
  "recall_row", "recall_col", "recall", "sensitive",
  "F1_row", "F1_col", "F1",
  "recovery_row", "recovery_col", "recovery",
  "relevance_row", "relevance_col", "relevance",
  "ARI_pre", "ARI_post", "discretization_penalty"
)

mean_row_col <- function(row_score, col_score) {
  mean(c(unname(row_score), unname(col_score)), na.rm = TRUE)
}

purity_score <- function(truth, pred) {
  tab <- table(truth, pred)
  sum(apply(tab, 2, max)) / length(truth)
}

set_recovery_relevance <- function(truth, pred) {
  truth_groups <- split(seq_along(truth), truth)
  pred_groups <- split(seq_along(pred), pred)
  jac <- vapply(pred_groups, function(pred_set) {
    vapply(truth_groups, function(true_set) {
      length(intersect(true_set, pred_set)) / length(union(true_set, pred_set))
    }, numeric(1))
  }, numeric(length(truth_groups)))
  if (is.null(dim(jac))) jac <- matrix(jac, nrow = length(truth_groups))
  list(
    recovery = mean(apply(jac, 1, max)),
    relevance = mean(apply(jac, 2, max))
  )
}

# Match predicted cluster identifiers to ground-truth identifiers with a
# Hungarian assignment, then calculate macro precision/recall/F1 and accuracy.
classification_scores <- function(truth, pred) {
  truth <- as.character(truth)
  pred <- as.character(pred)
  truth_levels <- unique(truth)
  pred_levels <- unique(pred)
  k <- max(length(truth_levels), length(pred_levels))
  score <- matrix(0, nrow = k, ncol = k)
  tab <- table(
    factor(pred, levels = pred_levels),
    factor(truth, levels = truth_levels)
  )
  score[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  assignment <- as.integer(solve_LSAP(score, maximum = TRUE))

  mapped <- rep(NA_character_, length(pred))
  pred_index <- match(pred, pred_levels)
  assigned_truth <- assignment[pred_index]
  valid <- assigned_truth <= length(truth_levels)
  mapped[valid] <- truth_levels[assigned_truth[valid]]
  mapped[!valid] <- paste0("unmatched_", pred[!valid])

  accuracy <- mean(mapped == truth)
  per_class <- vapply(truth_levels, function(g) {
    tp <- sum(mapped == g & truth == g)
    fp <- sum(mapped == g & truth != g)
    fn <- sum(mapped != g & truth == g)
    precision <- if (tp + fp == 0) 0 else tp / (tp + fp)
    recall <- if (tp + fn == 0) 0 else tp / (tp + fn)
    f1 <- if (precision + recall == 0) 0 else {
      2 * precision * recall / (precision + recall)
    }
    c(precision = precision, recall = recall, F1 = f1)
  }, numeric(3))

  c(
    Accuracy = accuracy,
    precision = mean(per_class["precision", ]),
    recall = mean(per_class["recall", ]),
    F1 = mean(per_class["F1", ])
  )
}

validate_labels <- function(labels, truth, axis) {
  if (length(labels) != length(truth)) {
    stop(axis, " predicted and true labels have different lengths.")
  }
  if (anyNA(labels) || any(labels == "")) stop(axis, " labels are incomplete.")
}

# The submission states that separate row and column scores are averaged. This
# function therefore uses an unweighted arithmetic mean, not an (m+p)-weighted
# mean as in the earlier development script.
overlapbic_metric <- function(biclabel, truelabel, rsim = NULL, csim = NULL,
                              m = NULL, p = NULL) {
  validate_labels(biclabel$row, truelabel$row, "Row")
  validate_labels(biclabel$col, truelabel$col, "Column")

  out <- setNames(rep(NA_real_, length(submission_metric_names)),
                  submission_metric_names)

  out["NMI_row"] <- NMI(truelabel$row, biclabel$row)
  out["NMI_col"] <- NMI(truelabel$col, biclabel$col)
  out["NMI"] <- mean_row_col(out["NMI_row"], out["NMI_col"])

  out["purity_row"] <- purity_score(truelabel$row, biclabel$row)
  out["purity_col"] <- purity_score(truelabel$col, biclabel$col)
  out["purity"] <- mean_row_col(out["purity_row"], out["purity_col"])

  out["ARI_row"] <- adjustedRandIndex(truelabel$row, biclabel$row)
  out["ARI_col"] <- adjustedRandIndex(truelabel$col, biclabel$col)
  out["ARI"] <- mean_row_col(out["ARI_row"], out["ARI_col"])

  class_row <- classification_scores(truelabel$row, biclabel$row)
  class_col <- classification_scores(truelabel$col, biclabel$col)
  for (metric in c("Accuracy", "precision", "recall", "F1")) {
    out[paste0(metric, "_row")] <- class_row[metric]
    out[paste0(metric, "_col")] <- class_col[metric]
    out[metric] <- mean_row_col(class_row[metric], class_col[metric])
  }
  out["sensitive"] <- out["recall"]

  set_row <- set_recovery_relevance(truelabel$row, biclabel$row)
  set_col <- set_recovery_relevance(truelabel$col, biclabel$col)
  for (metric in c("recovery", "relevance")) {
    out[paste0(metric, "_row")] <- set_row[[metric]]
    out[paste0(metric, "_col")] <- set_col[[metric]]
    out[metric] <- mean_row_col(set_row[[metric]], set_col[[metric]])
  }

  out["ARI_post"] <- out["ARI"]
  out
}

# Compatibility wrapper used by legacy supplementary scripts.
metric_micro <- function(allclus, label, m = NULL, p = NULL) {
  row <- classification_scores(label$row, allclus$row)
  col <- classification_scores(label$col, allclus$col)
  list(
    accurate = mean_row_col(row["Accuracy"], col["Accuracy"]),
    precision = mean_row_col(row["precision"], col["precision"]),
    sensitive = mean_row_col(row["recall"], col["recall"]),
    F1 = mean_row_col(row["F1"], col["F1"])
  )
}

excbic_metric <- function(truebic, cluslist, m, p, rowend = NULL) {
  truth <- list(row = integer(m), col = integer(p))
  for (i in seq_along(truebic$bicr)) truth$row[truebic$bicr[[i]]] <- i
  for (i in seq_along(truebic$bicc)) truth$col[truebic$bicc[[i]]] <- i
  metrics <- overlapbic_metric(cluslist, truth)
  metrics[c("recovery", "relevance", "recovery_row", "recovery_col",
            "relevance_row", "relevance_col", "precision", "sensitive",
            "F1", "Accuracy")]
}

sanitize_boundaries <- function(boundaries, n) {
  boundaries <- sort(unique(as.integer(boundaries)))
  boundaries[boundaries > 0L & boundaries < n]
}

# Convert shared-order boundaries into row and column labels. Segment IDs are
# retained on both node types so the labels genuinely derive from one partition.
labelinput <- function(simorder, matname, boundaries) {
  n <- length(simorder)
  cuts <- c(0L, sanitize_boundaries(boundaries, n), n)
  row_label <- rep(NA_integer_, length(matname$row))
  col_label <- rep(NA_integer_, length(matname$col))

  for (segment in seq_len(length(cuts) - 1L)) {
    members <- simorder[(cuts[segment] + 1L):cuts[segment + 1L]]
    row_index <- match(members, matname$row, nomatch = 0L)
    col_index <- match(members, matname$col, nomatch = 0L)
    row_label[row_index[row_index > 0L]] <- segment
    col_label[col_index[col_index > 0L]] <- segment
  }
  if (anyNA(row_label) || anyNA(col_label)) {
    stop("Joint order does not contain every row and column name.")
  }
  list(row = row_label, col = col_label)
}

labelinput2 <- function(simorder, matname, boundaries) {
  n <- length(simorder)
  cuts <- c(0L, sanitize_boundaries(boundaries, n), n)
  label <- rep(NA_integer_, length(matname))
  for (segment in seq_len(length(cuts) - 1L)) {
    members <- simorder[(cuts[segment] + 1L):cuts[segment + 1L]]
    index <- match(members, matname)
    label[index] <- segment
  }
  if (anyNA(label)) stop("Order does not contain every entity name.")
  label
}

# Oracle segmentation used for ARI_pre. The number and exact sizes of the true
# groups are known. Dynamic programming chooses the sequence of those group
# sizes that maximizes agreement with the labels along the recovered order.
oracle_partition <- function(order_names, truth, item_names) {
  truth <- as.character(truth)
  names(truth) <- item_names
  truth_ordered <- truth[order_names]
  if (anyNA(truth_ordered)) stop("Oracle order contains unknown entity names.")

  groups <- unique(truth)
  sizes <- vapply(groups, function(g) sum(truth == g), integer(1))
  k <- length(groups)
  n_states <- bitwShiftL(1L, k)
  score <- rep(-Inf, n_states)
  sequence_by_state <- vector("list", n_states)
  score[1L] <- 0
  sequence_by_state[[1L]] <- integer()

  for (mask in 0:(n_states - 2L)) {
    state_index <- mask + 1L
    if (!is.finite(score[state_index])) next
    used <- vapply(seq_len(k), function(g) bitwAnd(mask, bitwShiftL(1L, g - 1L)) != 0L,
                   logical(1))
    start <- sum(sizes[used]) + 1L
    for (g in which(!used)) {
      end <- start + sizes[g] - 1L
      gain <- sum(truth_ordered[start:end] == groups[g])
      next_mask <- bitwOr(mask, bitwShiftL(1L, g - 1L))
      next_index <- next_mask + 1L
      candidate <- score[state_index] + gain
      if (candidate > score[next_index]) {
        score[next_index] <- candidate
        sequence_by_state[[next_index]] <- c(sequence_by_state[[state_index]], g)
      }
    }
  }

  group_sequence <- sequence_by_state[[n_states]]
  ordered_prediction <- character(length(order_names))
  start <- 1L
  for (g in group_sequence) {
    end <- start + sizes[g] - 1L
    ordered_prediction[start:end] <- groups[g]
    start <- end + 1L
  }
  prediction <- ordered_prediction[match(item_names, order_names)]
  if (anyNA(prediction)) stop("Oracle segmentation failed to label all entities.")
  prediction
}

oracle_joint_partition <- function(joint_order, truelabel, matname) {
  item_names <- c(matname$row, matname$col)
  joint_truth <- c(truelabel$row, truelabel$col)
  joint_prediction <- oracle_partition(joint_order, joint_truth, item_names)
  list(
    row = joint_prediction[seq_along(matname$row)],
    col = joint_prediction[length(matname$row) + seq_along(matname$col)]
  )
}

add_discretization_metrics <- function(post_metric, pre_labels, truelabel) {
  pre_metric <- overlapbic_metric(pre_labels, truelabel)
  post_metric["ARI_pre"] <- pre_metric["ARI"]
  post_metric["ARI_post"] <- post_metric["ARI"]
  post_metric["discretization_penalty"] <-
    post_metric["ARI_pre"] - post_metric["ARI_post"]
  post_metric
}

seriout <- function(out, truelabel, rsim, csim, m, p, mat1,
                    r_to_py, find_auto_boundaries,
                    prominence = 0.01, distance = 5L, window = 10L,
                    smooth = "gaussian", sigma = 3) {
  rsim2 <- rsim[out$row_order, out$row_order, drop = FALSE]
  csim2 <- csim[out$col_order, out$col_order, drop = FALSE]

  boundaries_row <- find_auto_boundaries(
    r_to_py(rsim2), valley = "find_peaks", smooth = smooth,
    window = as.integer(window), sigma = sigma,
    prominence = prominence, distance = as.integer(distance)
  )
  boundaries_col <- find_auto_boundaries(
    r_to_py(csim2), valley = "find_peaks", smooth = smooth,
    window = as.integer(window), sigma = sigma,
    prominence = prominence, distance = as.integer(distance)
  )

  labels <- list(
    row = labelinput2(rownames(rsim2), rownames(mat1), boundaries_row),
    col = labelinput2(rownames(csim2), colnames(mat1), boundaries_col)
  )
  pre_labels <- list(
    row = oracle_partition(rownames(rsim2), truelabel$row, rownames(mat1)),
    col = oracle_partition(rownames(csim2), truelabel$col, colnames(mat1))
  )
  metric <- add_discretization_metrics(
    overlapbic_metric(labels, truelabel), pre_labels, truelabel
  )

  list(
    clus = labels, pre_clus = pre_labels, metric = metric,
    boundaries = list(row = as.integer(boundaries_row),
                      col = as.integer(boundaries_col))
  )
}
