
# FUNCTION FOR PROPERLY BREAKING LOGS #############################################

log10_pretty_breaks <- function(n = 3) {
  function(x) {
    10^pretty(log10(x), n = n)
  }
}
