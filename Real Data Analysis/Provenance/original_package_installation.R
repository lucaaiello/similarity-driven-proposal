library_dir <- Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(library_dir)) stop("R_LIBS_USER must be set explicitly")
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)

archives <- c(
  Rcpp = "https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.12.tar.gz",
  RcppArmadillo = "https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_15.2.3-1.tar.gz",
  RcppEigen = "https://cran.r-project.org/src/contrib/RcppEigen_0.3.4.0.2.tar.gz",
  RcppEnsmallen = "https://cran.r-project.org/src/contrib/Archive/RcppEnsmallen/RcppEnsmallen_0.2.22.1.1.tar.gz",
  huge = "https://cran.r-project.org/src/contrib/Archive/huge/huge_1.3.5.tar.gz"
)
expected <- c(
  Rcpp = "1.0.12",
  RcppArmadillo = "15.2.3.1",
  RcppEigen = "0.3.4.0.2",
  RcppEnsmallen = "0.2.22.1.1",
  huge = "1.3.5"
)

for (package in names(archives)) {
  installed <- tryCatch(
    as.character(packageVersion(package, lib.loc = library_dir)),
    error = function(e) NA_character_
  )
  if (!identical(installed, expected[[package]])) {
    message("Installing pinned ", package, " ", expected[[package]])
    install.packages(
      archives[[package]], lib = library_dir, repos = NULL,
      type = "source", dependencies = FALSE
    )
  }
}

actual <- vapply(
  names(expected),
  function(package) as.character(packageVersion(package, lib.loc = library_dir)),
  character(1)
)
print(data.frame(package = names(expected), expected, actual, row.names = NULL))
if (!identical(unname(actual), unname(expected))) {
  stop("pinned Hoffman2 package versions do not match")
}
