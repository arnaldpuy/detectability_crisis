
# FUNCTION TO COMPUTE HOW MUCH IRRIGATION LIES IN CELLS NOT IN THE K-CORE ######

compute_loss_for_k <- function(k) {
  grid_long[present == 1,  # cells where this dataset says "irrigated"
            {
              core_k <- n_pos >= k
              total_cells <- .N
              core_cells  <- sum(core_k, na.rm = TRUE)
              .(total_cells = total_cells,
                core_cells  = core_cells)
            },
            by = .(country, dataset)
  ][
    total_cells > 0,
    `:=`(share_lost = 1 - core_cells / total_cells, k = k)
  ]
}
