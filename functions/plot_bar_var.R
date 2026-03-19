
# FUNCTION TO CREATE BARPLOTS ##################################################

plot_bar_var <- function(freq_dt, varname) {
  tmp <- copy(freq_dt[variable == varname])
  ggplot(tmp, aes(x = fct_reorder(category, N), y = N)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "Number of datasets", title = varname) +
    theme_AP()
}
