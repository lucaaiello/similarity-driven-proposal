source("Real Data Analysis/Provenance/verify_sources.R")
# Rebuild manuscript assets from the completed downloaded Slurm summary only.
# No MCMC, optimizer, compilation, raw checkpoint reads, or cluster submission.
arguments <- commandArgs(trailingOnly=TRUE)
destination <- grep("^--output-root=",arguments,value=TRUE)
destination <- if(length(destination))sub("^--output-root=","",destination[1])else "reproduced_application"
input <- "Real Data Analysis/Results/summary.rds"
stopifnot(identical(readLines(file.path(dirname(input),"status.txt")),"COMPLETE"))
input_hash <- tools::md5sum(input); s <- readRDS(input)
stopifnot(length(s$records)==12L, s$plan$iterations==50000L,
  s$plan$burn==25000L, s$plan$thin==1L,
  identical(s$source_hashes,recorded_source_md5(names(s$source_hashes))))
x <- readRDS("Real Data Analysis/Data/Prepared/fresh_inputs.rds")
stopifnot(identical(dim(x$X),c(98L,117L)),identical(dim(x$Y),c(98L,30L)))
methods <- s$plan$methods; starts <- s$plan$start_names
method_names <- c("Similarity flip","Similarity + swap","Uniform flip")
start_names <- c("Null","5%","10%","15%")
start_tex <- c("Null","5\\%","10\\%","15\\%")
colors <- c("#0072B2","#009E73","#666666")
start_colors <- c("#0072B2","#D55E00","#009E73","#CC79A7")
images <- file.path(destination,"figures"); generated <- file.path(destination,"tables/application_latex")
out <- file.path(destination,"tables/application")
if(nzchar(destination)) {
  images<-file.path(destination,"figures")
  generated<-file.path(destination,"tables/application_latex")
  out<-file.path(destination,"tables/application")
}
for(d in c(images,generated,out))dir.create(d,recursive=TRUE,showWarnings=FALSE)
rr <- function(m)s$records[vapply(s$records,function(r)r$method==m,TRUE)]
N <- 25000L; it <- 25001:50000
bind <- function(v){d<-do.call(rbind,v);rownames(d)<-NULL;d}
pip_array <- simplify2array(lapply(s$records,function(r)r$pips))
min_pip <- apply(pip_array,c(1,2),min); max_pip <- apply(pip_array,c(1,2),max)
selected <- which(min_pip>=.5,arr.ind=TRUE)
stopifnot(nrow(selected)==7L)
labels <- read.table("Real Data Analysis/Data/ffq.labels2.txt",header=TRUE,stringsAsFactors=FALSE)
associations <- data.frame(genus=sub("^.*\\.","",x$taxa[selected[,2]]),
  category=selected[,2],predictor=selected[,1],code=colnames(x$X)[selected[,1]],
  min_all12=min_pip[selected],max_all12=max_pip[selected])
associations$nutrient<-as.character(labels[associations$code,1])
associations$nutrient[associations$code=="f204"]<-"Arachidonic fatty acid"
associations$nutrient[associations$code=="ucat"]<-"Catechin"
associations<-associations[order(associations$category,associations$predictor),]
pip_ranges<-category_stats<-list()
for(m in 1:3) {
  pa<-simplify2array(lapply(rr(methods[m]),function(r)r$pips))
  for(j in 1:30) {
    span<-apply(pa[,j,],1,function(v)diff(range(v)))
    category_stats[[length(category_stats)+1L]]<-data.frame(method=methods[m],
      genus=sub("^.*\\.","",x$taxa[j]),category=j,max_spread=max(span),
      spread_ge_half=sum(span>=.5),spread_ge_quarter=sum(span>=.25))
  }
  for(i in seq_len(nrow(associations))) {
    a<-associations[i,]; v<-pa[a$predictor,a$category,]
    pip_ranges[[length(pip_ranges)+1L]]<-data.frame(method=methods[m],genus=a$genus,
      nutrient=a$nutrient,code=a$code,association=i,min_pip=min(v),mean_pip=mean(v),max_pip=max(v))
  }
}
pip_ranges<-bind(pip_ranges); category_stats<-bind(category_stats)
stopifnot(requireNamespace("posterior",quietly=TRUE))
series<-function(r,field)if(field=="other_model_size")r$model_size-r$prev_size else r[[field]]
comparison<-bind(lapply(1:3,function(m) {
  z<-s$tables$convergence[s$tables$convergence$method==methods[m],]
  getr<-function(v)z$rhat[z$variable==v]
  d<-category_stats[category_stats$method==methods[m]&category_stats$genus!="Prevotella",]
  data.frame(method=methods[m],label=method_names[m],
    median_hours=median(s$tables$chain_metrics$seconds_including_startup[s$tables$chain_metrics$method==methods[m]])/3600,
    model_size_rhat=getr("model_size"),prev_size_rhat=getr("prev_size"),
    other_model_size_rhat=posterior::rhat(vapply(rr(methods[m]),function(r)series(r,"other_model_size"),numeric(N))),
    subject4014_rhat=getr("subject4014"),other_genera_median_max_spread=median(d$max_spread),
    all_genera_median_max_spread=median(category_stats$max_spread[category_stats$method==methods[m]]),
    prev_spread_ge_half=category_stats$spread_ge_half[category_stats$method==methods[m]&category_stats$genus=="Prevotella"])
}))
stopifnot(all(category_stats$spread_ge_half[category_stats$genus!="Prevotella"]==0))
write.csv(comparison,file.path(out,"method_comparison.csv"),row.names=FALSE)
write.csv(associations,file.path(out,"shared_empirical_support.csv"),row.names=FALSE)
write.csv(pip_ranges,file.path(out,"association_method_ranges.csv"),row.names=FALSE)
write.csv(category_stats,file.path(out,"category_pip_spreads.csv"),row.names=FALSE)
fmt<-function(v,d=3)sprintf(paste0("%.",d,"f"),v)
row_end<-strrep(intToUtf8(92),2)
tex<-c("% Generated from the saved application summary; do not edit numbers manually.",
  "\\begin{table}[t]","\\centering","\\small","\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Model-size $\\widehat R$ and median genus-wise maximum PIP spread, with and without \\textit{Prevotella}.}",
  "\\label{tab:real_data_comparison}","\\begin{tabular}{lrrrr}","\\toprule",
  paste0(" & \\multicolumn{2}{c}{Model-size $\\widehat R$} & \\multicolumn{2}{c}{PIP spread} ",row_end),
  paste0("Method & All genera & Other 29 & All genera & Other 29 ",row_end),"\\midrule")
for(i in 1:3)tex<-c(tex,paste0(comparison$label[i]," & ",fmt(comparison$model_size_rhat[i]),
  " & ",fmt(comparison$other_model_size_rhat[i]),
  " & ",fmt(comparison$all_genera_median_max_spread[i]),
  " & ",fmt(comparison$other_genera_median_max_spread[i])," ",row_end))
tex<-c(tex,"\\bottomrule","\\end{tabular}","\\end{table}")
writeLines(tex,file.path(generated,"real_data_comparison.tex"))
tex<-c("% Generated empirical support ranges; no FDR or convergence guarantee.",
  "\\begin{table}[t]","\\centering","\\small",
  "\\caption{PIP ranges across four starts for seven genus--nutrient associations with PIP at least 0.5 in all twelve chains.}",
  "\\label{tab:vs_taxa_covariate}","\\resizebox{\\textwidth}{!}{%","\\begin{tabular}{llrrr}","\\toprule",
  paste0("Genus & Nutrient & Similarity flip & Similarity + swap & Uniform flip ",row_end),"\\midrule")
for(i in seq_len(nrow(associations))) {
  a<-associations[i,]; ranges<-vapply(methods,function(m) {
    d<-pip_ranges[pip_ranges$method==m&pip_ranges$association==i,]
    paste0(fmt(d$min_pip),"--",fmt(d$max_pip))
  },"")
  tex<-c(tex,paste0("\\textit{",a$genus,"} & ",a$nutrient," & ",paste(ranges,collapse=" & ")," \\\\"))
}
tex<-c(tex,"\\bottomrule","\\end{tabular}","}","\\end{table}")
writeLines(tex,file.path(generated,"real_data_associations.tex"))
tex<-c("% Generated chain-specific retained summaries.","\\begin{table}[t]","\\centering","\\small",
  "\\caption{Ladder wall time, retained mean model size, and exchange acceptance between the two coldest replicas ($A_{\\rm cold}$).}",
  "\\label{tab:real_data_chain_diagnostics}","\\begin{tabular}{llrrr}","\\toprule",
  paste0("Method & Start & Hours & Model size & $A_{\\rm cold}$ ",row_end),"\\midrule")
for(m in 1:3)for(k in 1:4) {
  r<-rr(methods[m])[[k]]; d<-s$tables$chain_metrics
  d<-d[d$method==methods[m]&d$start==starts[k],]
  e<-s$tables$exchange_edges
  e<-e[e$method==methods[m]&e$start==starts[k]&e$edge==1,]
  tex<-c(tex,paste0(method_names[m]," & ",start_tex[k]," & ",fmt(d$seconds_including_startup/3600,2),
    " & ",fmt(mean(r$model_size),2),
    " & ",fmt(e$acceptance_post)," \\\\"))
}
tex<-c(tex,"\\bottomrule","\\end{tabular}","\\end{table}")
writeLines(tex,file.path(generated,"real_data_chain_diagnostics.tex"))

stopifnot(requireNamespace("posterior",quietly=TRUE))
scalar_fields<-c("model_size","other_model_size","prev_size","prev_alpha","log_posterior","positive_subject_mean")
scalar_labels<-c("Overall model size","Other 29 genera model size","Prevotella model size","Prevotella intercept","Log posterior",
  "Mean Prevotella proportion")
scalar<-list()
for(m in 1:3)for(v in seq_along(scalar_fields)) {
  mat<-vapply(rr(methods[m]),function(r)series(r,scalar_fields[v]),numeric(N))
  h<-posterior::rhat(mat)
  stored<-s$tables$convergence
  stored<-stored$rhat[stored$method==methods[m]&stored$variable==scalar_fields[v]]
  if(scalar_fields[v]!="other_model_size")stopifnot(length(stored)==1L,abs(h-stored)<1e-10)
  scalar[[length(scalar)+1L]]<-data.frame(method=methods[m],variable=scalar_fields[v],rhat=h,
    bulk_ess=posterior::ess_bulk(mat),tail_ess=posterior::ess_tail(mat))
}
scalar<-bind(scalar)
write.csv(scalar,file.path(out,"scalar_diagnostics.csv"),row.names=FALSE)
tex<-c("% Generated diagnostics from all four retained cold chains per method.",
  "\\begin{table}[t]","\\centering","\\small",
  "\\caption{Rank-normalized folded split $\\widehat R$ and bulk/tail ESS across four cold chains per method.}",
  "\\label{tab:real_data_scalar_diagnostics}","\\renewcommand{\\arraystretch}{0.9}",
  "\\fontsize{10}{11}\\selectfont","\\begin{tabular}{llrrr}","\\toprule",
  paste0("Method & Quantity & $\\widehat R$ & Bulk ESS & Tail ESS ",row_end),"\\midrule")
for(m in 1:3)for(v in which(scalar_fields!="other_model_size")) {
  d<-scalar[scalar$method==methods[m]&scalar$variable==scalar_fields[v],]
  tex<-c(tex,paste0(if(v==1)method_names[m]else""," & ",scalar_labels[v]," & ",fmt(d$rhat),
    " & ",fmt(d$bulk_ess,1)," & ",fmt(d$tail_ess,1)," ",row_end))
  if(v==length(scalar_fields)&&m<3)tex<-c(tex,"\\addlinespace")
}
tex<-c(tex,"\\bottomrule","\\end{tabular}","\\end{table}")
writeLines(tex,file.path(generated,"real_data_scalar_diagnostics.tex"))

source("Real Data Analysis/Figures_and_Tables/plot_figures.R")
stopifnot(identical(input_hash,tools::md5sum(input)),
  identical(s$source_hashes,recorded_source_md5(names(s$source_hashes))))
saveRDS(list(input_hash=input_hash,source_hashes=s$source_hashes,plan=s$plan,
  method_comparison=comparison,associations=associations,association_ranges=pip_ranges,
  category_spreads=category_stats,scalar_diagnostics=scalar,
  posterior_version=as.character(packageVersion("posterior"))),file.path(out,"provenance.rds"))
writeLines(c("Generated from completed Slurm summary; no MCMC, fitting or compilation.",
  paste("Input MD5:",unname(input_hash)),"PIP support criterion: empirical PIP >=0.5 in all12 cold chains; no FDR claim.",
  "All methods share adaptive parallel tempering; runtime descriptive, not a controlled efficiency benchmark."),file.path(out,"README.txt"))
cat("REAL-DATA MANUSCRIPT ASSETS COMPLETE\n",out,"\n")
