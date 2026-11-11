# Optional setup inside a compute allocation; not called by postprocessing.
lib<-Sys.getenv("R_LIBS_USER");stopifnot(nzchar(lib));dir.create(lib,recursive=TRUE,showWarnings=FALSE)
urls<-c(
 Rcpp="https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.12.tar.gz",
 RcppArmadillo="https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_15.2.3-1.tar.gz",
 RcppEnsmallen="https://cran.r-project.org/src/contrib/Archive/RcppEnsmallen/RcppEnsmallen_0.2.22.1.1.tar.gz")
for(url in urls)install.packages(url,repos=NULL,type="source",lib=lib)
expected<-c(Rcpp="1.0.12",RcppArmadillo="15.2.3.1",RcppEnsmallen="0.2.22.1.1")
actual<-vapply(names(expected),function(p)as.character(packageVersion(p,lib.loc=lib)),"")
stopifnot(identical(actual,expected));print(actual)
