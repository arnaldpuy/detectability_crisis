
plot_attr <- function(r) {
  
  mc <- r$match_col
  
  p1 <- ggplot(r$pairs, aes(.data[[mc]], disagreement, fill = .data[[mc]])) +
    geom_boxplot(outlier.size = 0.8) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
    scale_fill_manual(values = c("different" = "#F8D7DA", "same" = "#D5E8D4")) +
    labs(title = sprintf("Pairwise disagreement: %s", r$attribute),
         subtitle = sprintf("Wilcoxon p = %.3f", r$pooled_p),
         x = NULL, y = "Disagreement rate") +
    theme_AP() + theme(plot.subtitle = element_text(size = 7)) +
    guides(fill = "none")
  
  p2 <- ggplot(r$res_pairs, aes(.data[[mc]], disagreement, fill = .data[[mc]])) +
    geom_boxplot(outlier.size = 0.8) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
    scale_fill_manual(values = c("different" = "#F8D7DA", "same" = "#D5E8D4")) +
    facet_wrap(~resolution) +
    labs(title = sprintf("%s × resolution", r$attribute),
         x = NULL, y = "Disagreement rate") +
    theme_AP() + guides(fill = "none")
  
  p3 <- ggplot(r$res_wilcox, aes(x = resolution, y = p_value)) +
    geom_col(fill = "#5B8DB8", width = 0.5, alpha = 0.85) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
    geom_text(aes(label = sprintf("Delta = %.3f", median_diff - median_same)),
              vjust = -0.5, size = 2.2) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), expand = expansion(mult = c(0, 0.15))) +
    labs(title = sprintf("Wilcoxon p: %s", r$attribute),
         x = NULL, y = "p-value") +
    theme_AP()
  
  p4 <- ggplot(data.table(delta = r$null_dist), aes(x = delta)) +
    geom_histogram(bins = 60, fill = "grey70", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = r$obs_delta, colour = "#C0392B", linewidth = 0.8) +
    annotate("text", x = r$obs_delta, y = Inf,
             label = sprintf("Delta = %.4f\np = %.3f", r$obs_delta, r$perm_p),
             vjust = 2, hjust = -0.1, colour = "#C0392B", size = 2.2) +
    labs(title = sprintf("Permutation null: %s", r$attribute),
         x = "Delta disagreement", y = "Count") +
    theme_AP()
  
  plot_grid(p1, p2, p3, p4, ncol = 2, labels = "auto", label_size = 8)
}