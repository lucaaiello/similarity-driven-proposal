# Integrity check only; base R, no fitting/sampling/compilation.
m<-read.csv("reproduce/inputs_md5.csv",stringsAsFactors=FALSE)
if(!all(file.exists(m$file)))stop("Missing saved input: ",paste(m$file[!file.exists(m$file)],collapse=", "))
actual<-unname(tools::md5sum(m$file))
if(!identical(actual,m$md5))stop("Input checksum mismatch: ",paste(m$file[actual!=m$md5|is.na(actual)],collapse=", "))
cat("Saved input integrity checks passed:",nrow(m),"files\n")
