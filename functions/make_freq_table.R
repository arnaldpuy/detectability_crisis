

# FUNCTION TO MAKE FREQUENCY TABLE #############################################

make_freq_table <- function(dt, varname) {
  out <- dt[, .N, by = .(category = get(varname))][order(-N)]
  out[, variable:= varname]
  out[, share:= N / sum(N)]
  out[]
}
