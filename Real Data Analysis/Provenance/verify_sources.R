# Verify historical source receipts against original archived bytes. These
# checks are distinct from the current-source gates used for new sampler runs.
recorded_source_md5 <- function(paths) {
  archive <- "Real Data Analysis/Provenance/original_sources.zip"
  stopifnot(file.exists(archive),is.character(paths),length(paths)>0L,
    !any(grepl("(^|/)\\.\\.(/|$)|^[/\\\\]|^[A-Za-z]:",paths)))
  members <- utils::unzip(archive,list=TRUE)$Name
  if(!all(paths%in%members))stop("Recorded source absent from original archive: ",paste(paths[!paths%in%members],collapse=", "))
  directory <- tempfile("original_sources_");dir.create(directory)
  parent<-normalizePath(tempdir(),winslash="/",mustWork=TRUE)
  resolved<-normalizePath(directory,winslash="/",mustWork=TRUE)
  stopifnot(startsWith(resolved,paste0(parent,"/")))
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  utils::unzip(archive,files=paths,exdir=directory)
  values<-unname(tools::md5sum(file.path(directory,paths)))
  if(anyNA(values))stop("Could not hash extracted original sources")
  setNames(values,paths)
}
verify_recorded_sources <- function(expected) {
  stopifnot(!is.null(names(expected)))
  identical(expected,recorded_source_md5(names(expected)))
}
recorded_file_matches <- function(expected,current) {
  verify_recorded_sources(expected) && identical(unname(expected),unname(tools::md5sum(current)))
}
