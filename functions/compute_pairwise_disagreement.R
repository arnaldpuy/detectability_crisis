
# FUNCTION TO COMPUTE PAIRWISE DISAGREEMENT ####################################

compute_pairwise_disagreement <- function(dtp, dataset_names) {

  pairs <- as.data.table(t(combn(dataset_names, 2)))
  setnames(pairs, c("dataset_i", "dataset_j"))

  out_list <- vector("list", nrow(pairs))

  for (k in seq_len(nrow(pairs))) {
    di <- pairs$dataset_i[k]
    dj <- pairs$dataset_j[k]

    tmp <- dtp[, .(xi = get(di), xj = get(dj))]

    # keep comparable non-missing cells--------------

    tmp <- tmp[!is.na(xi) & !is.na(xj)]

    n_cells <- nrow(tmp)

    # basic pairwise metrics ------------------------ ---

    n_00 <- tmp[, sum(xi == 0 & xj == 0)]
    n_11 <- tmp[, sum(xi == 1 & xj == 1)]
    n_10 <- tmp[, sum(xi == 1 & xj == 0)]
    n_01 <- tmp[, sum(xi == 0 & xj == 1)]

    n_disagree <- n_10 + n_01
    share_disagree_all <- n_disagree / n_cells

    # conditional on at least one dataset detecting irrigation

    n_union <- n_11 + n_10 + n_01
    share_disagree_union <- ifelse(n_union > 0, n_disagree / n_union, NA_real_)

    # Jaccard similarity / dissimilarity for presence masks

    jaccard_similarity <- ifelse(n_union > 0, n_11 / n_union, NA_real_)
    jaccard_dissimilarity <- ifelse(n_union > 0, 1 - jaccard_similarity, NA_real_)

    out_list[[k]] <- data.table(
      dataset_i = di,
      dataset_j = dj,
      n_cells = n_cells,
      n_00 = n_00,
      n_11 = n_11,
      n_10 = n_10,
      n_01 = n_01,
      n_disagree = n_disagree,
      share_disagree_all = share_disagree_all,
      share_disagree_union = share_disagree_union,
      jaccard_similarity = jaccard_similarity,
      jaccard_dissimilarity = jaccard_dissimilarity
    )
  }

  rbindlist(out_list)
}
