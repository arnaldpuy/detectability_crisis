
# Label pairs, run Wilcoxon + permutation for one attribute ####################

run_test <- function(a, pairs_dt, res_dt, dt_meta, n_perm = 10000) {
  lookup <- dt_meta[, .(dataset, val = get(a))]
  mc <- paste0(a, "_match")
  
  # Label pooled pairs---------------------------
  
  p <- copy(pairs_dt)
  p[lookup, on = .(ds1 = dataset), v1:= i.val]
  p[lookup, on = .(ds2 = dataset), v2:= i.val]
  p[, (mc) := fifelse(v1 == v2, "same", "different")]
  
  # Pooled Wilcoxon------------------------------
  
  s <- p[get(mc) == "same", disagreement]
  d <- p[get(mc) == "different", disagreement]
  wt <- wilcox.test(d, s, alternative = "greater")
  
  # Label resolution-stratified pairs-----------
  
  r <- copy(res_dt)
  r[lookup, on = .(ds1 = dataset), v1:= i.val]
  r[lookup, on = .(ds2 = dataset), v2:= i.val]
  r[, (mc):= fifelse(v1 == v2, "same", "different")]
  
  # Resolution-stratified Wilcoxon-------------
  
  res_w <- r[, {
    s <- .SD[get(mc) == "same", disagreement]
    d <- .SD[get(mc) == "different", disagreement]
    
    if (length(s) < 2 || length(d) < 2) {
      
      data.table(n_same = length(s), n_diff = length(d),
                 median_same = median(s), median_diff = median(d),
                 W = NA_real_, p_value = NA_real_)
    } else {
      
      wt <- wilcox.test(d, s, alternative = "greater")
      data.table(n_same = length(s), n_diff = length(d),
                 median_same = median(s), median_diff = median(d),
                 W = wt$statistic, p_value = wt$p.value)
    }
  }, by = resolution]
  
  # Permutation test (resolution-blocked)-------
  
  obs_delta <- r[, .(
    delta = mean(.SD[get(mc) == "different", disagreement]) -
      mean(.SD[get(mc) == "same", disagreement])
  ), by = resolution][, mean(delta, na.rm = TRUE)]
  
  null_dist <- replicate(n_perm, {
    pl <- data.table(dataset = dt_meta$dataset, val = sample(dt_meta[[a]]))
    tmp <- copy(r)
    tmp[pl, on = .(ds1 = dataset), v1:= i.val]
    tmp[pl, on = .(ds2 = dataset), v2:= i.val]
    tmp[, pm := fifelse(v1 == v2, "same", "different")]
    tmp[, .(delta = mean(.SD[pm == "different", disagreement]) -
              mean(.SD[pm == "same", disagreement])
    ), by = resolution][, mean(delta, na.rm = TRUE)]
  })
  
  list(
    attribute = a,
    pairs = p,
    res_pairs = r,
    res_wilcox = res_w,
    pooled_W = wt$statistic,
    pooled_p = wt$p.value,
    obs_delta = obs_delta,
    null_dist = null_dist,
    perm_p = mean(null_dist >= obs_delta),
    match_col = mc
  )
}