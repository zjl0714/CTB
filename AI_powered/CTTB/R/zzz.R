# nocov start
.onLoad <- function(libname, pkgname) {
  utils::globalVariables(c(
    "Method", "Estimate", "Type", "CI95_lower", "CI95_upper",
    "CI95_lower_IM", "CI95_upper_IM"
  ))
}
# nocov end
