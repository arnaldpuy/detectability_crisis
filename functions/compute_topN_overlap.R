
# FUNCTION TO COMPUTE PAIRWISE OVERLAP OF TOP-N IRRIGATION COUNTRIES ###########

compute_topN_overlap <- function(country_totals, dataset_names, N = 20) {

  topN_list <- lapply(dataset_names, function(d) {
    country_totals[dataset == d][order(-irrigated_area_mha)][1:N, country]
  })

  names(topN_list) <- dataset_names

  # Jaccard & intersection size for each pair
  pairs <- combn(dataset_names, 2, simplify = FALSE)

  rbindlist(lapply(pairs, function(p) {
    d1 <- p[1]
    d2 <- p[2]
    c1 <- topN_list[[d1]]
    c2 <- topN_list[[d2]]

    inter <- intersect(c1, c2)
    union <- union(c1, c2)

    data.table(dataset_1 = d1,
               dataset_2 = d2,
               N = N,
               intersection_N = length(inter),
               jaccard_topN = length(inter) / length(union))
  }))
}
