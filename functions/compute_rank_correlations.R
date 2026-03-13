
# FUNCTION TO COMPUTE PAIRWISE RANK CORRELATIONS BETWEEN DATASETS #############@

compute_rank_correlations <- function(country_totals) {
  pairs <- combn(unique(country_totals$dataset), 2, simplify = FALSE)

  rbindlist(lapply(pairs, function(p) {
    d1 <- p[1]
    d2 <- p[2]

    x <- country_totals[dataset == d1, .(country, rank1 = rank)]
    y <- country_totals[dataset == d2, .(country, rank2 = rank)]

    m <- merge(x, y, by = "country", all = FALSE)

    data.table(dataset_1 = d1,
               dataset_2 = d2,
               kendall = cor(m$rank1, m$rank2, method = "kendall"),
               spearman = cor(m$rank1, m$rank2, method = "spearman"))
  }))
}
