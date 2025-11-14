.onLoad <- function(libname, pkgname) {
  .mz_comp_db <<- memoise::memoise(.mz_comp_db)
}