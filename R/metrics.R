################################################################################
# BiSer Benchmarking: Evaluation Metrics
#
# Clustering and seriation quality metrics:
#   NMI, ARI, Purity, Path length, Smoothness, Bandwidth, Block contrast,
#   Recovery, Relevance, Precision, Sensitivity, F1, Accuracy
#
# Dependencies: mclust, aricode, clue
################################################################################

library(mclust)
library(aricode)
library(clue)

# ==============================================================================
# Overlap-aware biclustering metric (row + column combined)
# ==============================================================================
overlapbic_metric <- function(biclabel, truelabel, rsim, csim, m, p) {
  rsim <- rsim[biclabel$row, biclabel$row]
  csim <- csim[biclabel$col, biclabel$col]
  metric <- numeric()

  metric["NMI_row"] <- NMI(truelabel$row, biclabel$row)
  metric["NMI_col"] <- NMI(truelabel$col, biclabel$col)
  metric["NMI"]     <- (metric["NMI_row"]*m + metric["NMI_col"]*p) / (m+p)

  purity_fn <- function(truth, pred) {
    tab <- table(truth, pred)
    sum(apply(tab, 2, max)) / length(truth)
  }
  metric["purity_row"] <- purity_fn(truelabel$row, biclabel$row)
  metric["purity_col"] <- purity_fn(truelabel$col, biclabel$col)
  metric["purity"]     <- (metric["purity_row"]*m + metric["purity_col"]*p) / (m+p)

  path_length <- function(dist_mat, ord) {
    total <- 0
    for (i in 1:(length(ord)-1)) total <- total + dist_mat[ord[i], ord[i+1]]
    total
  }
  metric["pathlen_row"] <- path_length(1 - rsim, 1:m)
  metric["pathlen_col"] <- path_length(1 - csim, 1:p)
  metric["pathlen"]     <- (metric["pathlen_row"]*m + metric["pathlen_col"]*p) / (m+p)

  metric["ARI_row"] <- adjustedRandIndex(truelabel$row, biclabel$row)
  metric["ARI_col"] <- adjustedRandIndex(truelabel$col, biclabel$col)
  metric["ARI"]     <- (m*metric["ARI_row"] + p*metric["ARI_col"]) / (m+p)

  profile_smooth <- function(mat, window = 10) {
    n <- nrow(mat); profile <- numeric(n)
    for (i in 1:n) {
      lo <- max(1, i-window); hi <- min(n, i+window)
      profile[i] <- mean(mat[i, lo:hi])
    }
    var(diff(profile))
  }
  metric["smooth_row"] <- profile_smooth(rsim, 10)
  metric["smooth_col"] <- profile_smooth(csim, 10)
  metric["smooth"]     <- (m*metric["smooth_row"] + p*metric["smooth_col"]) / (m+p)

  bw_fn <- function(mat) {
    n <- nrow(mat)
    sum(mat * abs(row(mat) - col(mat))) / (n^2)
  }
  metric["bandwidth_row"] <- bw_fn(rsim)
  metric["bandwidth_col"] <- bw_fn(csim)
  metric["bandwidth"]     <- (m*metric["bandwidth_row"] + p*metric["bandwidth_col"]) / (m+p)

  bc_fn <- function(mat, band_width = 10) {
    diag_band <- abs(row(mat) - col(mat)) <= band_width
    off_band  <- abs(row(mat) - col(mat)) >  band_width
    mean(mat[diag_band]) / mean(mat[off_band])
  }
  metric["block_contrast_row"] <- bc_fn(rsim)
  metric["block_contrast_col"] <- bc_fn(csim)
  metric["block_contrast"]     <- (m*metric["block_contrast_row"] +
                                    p*metric["block_contrast_col"]) / (m+p)
  return(metric)
}

# ==============================================================================
# Micro-averaged precision / sensitivity / F1 / accuracy
# ==============================================================================
metric_micro <- function(allclus, label, m, p) {
  match_label <- function(clus1, label1) {
    label1 <- as.factor(label1); label2 <- levels(label1)
    ktrue <- length(label2); k <- length(unique(clus1))
    ratio <- matrix(nrow=k, ncol=ktrue,
                    dimnames=list(paste0("clus",1:k), label2))
    for (j in 1:k) for (t in 1:ktrue)
      ratio[j,t] <- length(intersect(which(clus1==j),which(label1==label2[t])))/
                     length(label1==label2[t])
    if (k<=ktrue) return(solve_LSAP(ratio, maximum=TRUE))
    matchout <- numeric()
    tmp <- solve_LSAP(t(ratio), maximum=TRUE)
    for (i in 1:ktrue) matchout[tmp[i]] <- i
    res <- (1:k)[-tmp]
    for (i in seq_along(res)) matchout[res[i]] <- ktrue + i
    matchout
  }
  matchrow <- match_label(allclus$row, label$row)
  matchcol <- match_label(allclus$col, label$col)
  ktrue <- length(unique(label$row))
  k <- min(length(unique(allclus$row)), length(unique(allclus$col)))
  truemat <- clusmat <- matrix(0, nrow=m, ncol=p)
  for (i in 1:ktrue) truemat[which(label$row==i), which(label$col==i)] <- i
  for (i in 1:k)
    clusmat[which(allclus$row==which(matchrow==i)),
            which(allclus$col==which(matchcol==i))] <- i
  accurate <- sum(as.numeric(truemat)==as.numeric(clusmat)) / (m*p)
  tp <- tn <- fp <- fn <- numeric()
  for (i in 1:ktrue) {
    clus1 <- true1 <- rep(0, m*p)
    clus1[which(as.numeric(clusmat)==i)] <- i
    true1[which(as.numeric(truemat)==i)] <- i
    tab <- table(true1, clus1)
    if (sum(as.numeric(clusmat)==i)==0) tab <- cbind(tab, rep(0,2))
    tp[i]<-tab[2,2]; tn[i]<-tab[1,1]; fp[i]<-tab[1,2]; fn[i]<-tab[2,1]
  }
  tp<-sum(tp); tn<-sum(tn); fp<-sum(fp); fn<-sum(fn)
  prec <- tp/(tp+fp); sens <- tp/(tp+fn)
  list(accurate=accurate, precision=prec, sensitive=sens,
       F1=2*(prec*sens)/(prec+sens))
}

# ==============================================================================
# Exclusive bicluster metrics (recovery, relevance, precision, F1)
# ==============================================================================
excbic_metric <- function(truebic, cluslist, m, p, rowend) {
  metric2 <- c("recovery","relevance","recovery_row","recovery_col",
               "relevance_row","relevance_col","precision","sensitive","F1","Accuracy")
  grouptmp <- list(row=numeric(), col=numeric())
  for (i in 1:(length(rowend)-1)) {
    grouptmp$row[truebic$bicr[[paste0("bic",i)]]] <- i
    grouptmp$col[truebic$bicc[[paste0("bic",i)]]] <- i
  }
  metricout <- setNames(rep(NA, length(metric2)), metric2)

  # Jaccard-based recovery/relevance for rows
  ktrue <- length(truebic$bicr); clustmp <- cluslist$row
  kmeth <- length(unique(clustmp))
  jac <- matrix(nrow=ktrue, ncol=kmeth)
  for (j in 1:kmeth) { ci <- which(clustmp==j)
    for (t in 1:ktrue) { bi <- truebic$bicr[[paste0("bic",t)]]
      jac[t,j] <- length(intersect(ci,bi))/length(union(ci,bi)) }}
  metricout["recovery_row"]  <- mean(apply(jac,1,max))
  metricout["relevance_row"] <- mean(apply(jac,2,max))

  # Jaccard-based recovery/relevance for cols
  ktrue <- length(truebic$bicc); clustmp <- cluslist$col
  kmeth <- length(unique(clustmp))
  jac <- matrix(nrow=ktrue, ncol=kmeth)
  for (j in 1:kmeth) { ci <- which(clustmp==j)
    for (t in 1:ktrue) { bi <- truebic$bicc[[paste0("bic",t)]]
      jac[t,j] <- length(intersect(ci,bi))/length(union(ci,bi)) }}
  metricout["recovery_col"]  <- mean(apply(jac,1,max))
  metricout["relevance_col"] <- mean(apply(jac,2,max))
  metricout["recovery"]  <- (m*metricout["recovery_row"]+p*metricout["recovery_col"])/(m+p)
  metricout["relevance"] <- (m*metricout["relevance_row"]+p*metricout["relevance_col"])/(m+p)

  out <- metric_micro(cluslist, grouptmp, m, p)
  metricout["precision"]<-out$precision; metricout["sensitive"]<-out$sensitive
  metricout["F1"]<-out$F1; metricout["Accuracy"]<-out$accurate
  return(metricout)
}

# ==============================================================================
# Label extraction utilities
# ==============================================================================

#' BiSer output: convert boundary + ordering to row/col cluster labels
labelinput <- function(simorder, matname, boundaries) {
  boundaries <- c(0, boundaries, length(simorder))
  group <- list(row=list(), col=list()); numr <- numc <- 1
  posr <- posc <- numeric()
  for (i in 1:(length(boundaries)-1)) {
    tmp <- simorder[(boundaries[i]+1):boundaries[i+1]]
    tmpr <- tmp %in% matname$row; tmpc <- tmp %in% matname$col
    if (sum(tmpr)>0) {
      group$row[[numr]] <- match(tmp[tmpr], matname$row)
      posr[numr] <- mean(group$row[[numr]]); numr <- numr+1 }
    if (sum(tmpc)>0) {
      group$col[[numc]] <- match(tmp[tmpc], matname$col)
      posc[numc] <- mean(group$col[[numc]]); numc <- numc+1 }
  }
  group$row <- group$row[order(posr)]; group$col <- group$col[order(posc)]
  label <- list(row=numeric(), col=numeric())
  for (i in seq_along(group$row)) label$row[group$row[[i]]] <- i
  for (i in seq_along(group$col)) label$col[group$col[[i]]] <- i
  return(label)
}

#' Single-dimension label extraction from boundaries
labelinput2 <- function(simorder, matname, boundaries) {
  boundaries <- c(0, boundaries, length(simorder))
  group <- list(); pos <- numeric()
  for (i in 1:(length(boundaries)-1)) {
    tmp <- simorder[(boundaries[i]+1):boundaries[i+1]]
    group[[i]] <- match(tmp, matname); pos[i] <- mean(group[[i]])
  }
  group <- group[order(pos)]; label <- numeric()
  for (i in seq_along(group)) label[group[[i]]] <- i
  return(label)
}

#' Evaluate seriation methods using boundary detection
seriout <- function(out, truelabel, rsim, csim, m, p, mat1,
                    r_to_py, find_auto_boundaries, prominence, distance) {
  rsim2 <- rsim[out$row_order, out$row_order]
  csim2 <- csim[out$col_order, out$col_order]
  py_rsim <- r_to_py(rsim2); py_csim <- r_to_py(csim2)
  ab_r <- find_auto_boundaries(py_rsim, valley="find_peaks",
                                prominence=prominence, distance=distance)
  labelr <- labelinput2(rownames(rsim2), rownames(mat1), ab_r)
  ab_c <- find_auto_boundaries(py_csim, valley="find_peaks",
                                prominence=prominence, distance=distance)
  labelc <- labelinput2(rownames(csim2), colnames(mat1), ab_c)
  clus <- list(row=labelr, col=labelc)
  metric <- overlapbic_metric(clus, truelabel, rsim, csim, m, p)
  list(clus=clus, metric=metric)
}
