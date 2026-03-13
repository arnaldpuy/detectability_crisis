
# FUNCTION TO EXTRACT CONTINENT FROM COUNTRY ###################################

country_to_continent<- function(dt, country_col = "country",
                                standardize_country_name = TRUE) {
  stopifnot(is.data.table(dt))
  stopifnot(country_col %in% names(dt))

  # work internally with a fixed name
  if (country_col != "country") {
    setnames(dt, country_col, "country")
  }

  # unique countries
  u <- unique(dt[!is.na(country), .(country)])

  # lookup once
  u[, code := countrycode(country,
                          origin = "country.name",
                          destination = "un")]
  u[, continent := countrycode(country,
                               origin = "country.name",
                               destination = "continent")]

  # standardize names if requested
  if (standardize_country_name) {
    u[, country_std := countrycode(code,
                                   origin = "un",
                                   destination = "country.name")]
  }

  # join back (overwrite country in place)
  setkey(u, country)
  dt[u, `:=`(code = i.code,
             continent = i.continent,
             country = if (standardize_country_name) i.country_std else country),
     on = "country"]

  # enforce lowercase column names
  setnames(dt, tolower(names(dt)))

  invisible(dt)
}

