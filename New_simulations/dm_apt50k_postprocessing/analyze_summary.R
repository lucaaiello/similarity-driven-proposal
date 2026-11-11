# Local read-only analysis of the downloaded completed Slurm summary.
# No MCMC, optimizer, compilation, checkpoint reads, or remote commands.
args <- commandArgs(TRUE)
root <- if (length(args)) args[1] else "New_simulations/outputs/dm_apt50k_fresh_slurm_v1/summary"
path <- file.path(root, "summary.rds")
stopifnot(requireNamespace("posterior", quietly = TRUE),
  identical(readLines(file.path(root, "status.txt")), "COMPLETE"))
before <- tools::md5sum(path)
s <- readRDS(path)
stopifnot(length(s$records) == 12L, length(s$checkpoint_hashes) == 12L,
  identical(s$source_hashes, tools::md5sum(names(s$source_hashes))),
  s$plan$iterations == 50000L, s$plan$burn == 25000L, s$plan$thin == 1L)
x <- readRDS("New_simulations/dm_apt50k_fresh_v1/inputs/fresh.rds")
methods <- s$plan$methods; starts <- s$plan$start_names
short <- c("Similarity flip", "Similarity local", "Uniform flip")
fields <- c("model_size", "prev_size", "prev_alpha", "log_posterior",
  "subject3026", "subject4014", "positive_subject_mean")
N <- 25000L; iteration <- seq.int(25001L, 50000L)
stopifnot(all(vapply(s$records, function(r) all(vapply(r[fields], length, 1L) == N) &&
  all(is.finite(unlist(r[fields]))) && all(dim(r$fitted) == c(98L, N)) &&
  all(is.finite(r$fitted)) && all(r$fitted >= 0 & r$fitted <= 1), TRUE)))
stopifnot(!anyDuplicated(vapply(s$records, function(r) paste(r$method, r$start), "")))
out <- file.path("New_simulations/outputs", paste0("dm_apt50k_assessment_",
  format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
stopifnot(!file.exists(out), dir.create(out, recursive = TRUE))
writeLines("INCOMPLETE", file.path(out, "status.txt"))
cat("LOCAL SUMMARY ANALYSIS; no MCMC, compilation or optimizer.\n", out, "\n")
bind <- function(v) { d <- do.call(rbind, v); rownames(d) <- NULL; d }
safe <- function(f, v) suppressWarnings(tryCatch(as.numeric(f(v)), error = function(e) NA_real_))
getrr <- function(method) s$records[vapply(s$records, function(r) r$method == method, TRUE)]
windows <- list(full = seq_len(N), first_half = 1:12500, last_half = 12501:25000,
  last5000 = 20001:25000)
diagnostics <- chain <- drift <- occupancy <- blocks <- fitted <- list()
for (method in methods) {
  rr <- getrr(method)
  stopifnot(identical(vapply(rr, function(r) as.integer(r$start), 1L), 1:4))
  for (field in fields) {
    mat <- vapply(rr, function(r) r[[field]], numeric(N))
    for (window in names(windows)) {
      v <- mat[windows[[window]], , drop = FALSE]
      diagnostics[[length(diagnostics)+1L]] <- data.frame(method, variable = field, window,
        draws_per_chain = nrow(v), rhat = safe(posterior::rhat, v),
        bulk_ess = safe(posterior::ess_bulk, v), tail_ess = safe(posterior::ess_tail, v),
        mcse_mean = safe(posterior::mcse_mean, v), range_chain_means = diff(range(colMeans(v))))
    }
    for (k in 1:4) {
      v <- mat[, k]; q <- quantile(v, c(.025, .5, .975), names = FALSE)
      chain[[length(chain)+1L]] <- data.frame(method, start = starts[k], variable = field,
        mean = mean(v), sd = sd(v), q025 = q[1], median = q[2], q975 = q[3],
        bulk_ess = safe(posterior::ess_bulk, matrix(v, ncol=1)),
        tail_ess = safe(posterior::ess_tail, matrix(v, ncol=1)),
        mcse_mean = safe(posterior::mcse_mean, matrix(v, ncol=1)))
      drift[[length(drift)+1L]] <- data.frame(method, start=starts[k], variable=field,
        first_half_mean=mean(v[1:12500]), last_half_mean=mean(v[12501:25000]),
        change_last_minus_first=mean(v[12501:25000])-mean(v[1:12500]),
        last5000_mean=mean(tail(v,5000)))
      for (b in 1:4) {
        ii <- seq.int((b-1L)*6250L+1L, b*6250L)
        blocks[[length(blocks)+1L]] <- data.frame(method,start=starts[k],variable=field,
          block=b, iteration_start=iteration[min(ii)], iteration_end=iteration[max(ii)], mean=mean(v[ii]))
      }
    }
  }
  for (k in 1:4) {
    r <- rr[[k]]
    for (subject in c("3026", "4014")) for (window in names(windows)) {
      v <- r[[paste0("subject", subject)]][windows[[window]]]; high <- v >= .05
      runs <- rle(high)
      occupancy[[length(occupancy)+1L]] <- data.frame(method,start=starts[k],subject,window,
        mean_percent=mean(v)*100, fraction_below5=mean(!high),
        end_of_sweep_threshold_crossings=sum(diff(high)!=0),
        longest_low_run=if(any(!runs$values))max(runs$lengths[!runs$values])else 0L,
        longest_high_run=if(any(runs$values))max(runs$lengths[runs$values])else 0L,
        first_percent=v[1]*100,last_percent=tail(v,1)*100)
    }
    for (j in seq_len(nrow(r$fitted))) {
      v <- r$fitted[j, ]
      fitted[[length(fitted)+1L]] <- data.frame(method,start=starts[k],subject=rownames(x$X)[j],
        observed_prevotella_count=x$Y[j, match("Prevotella",sub("^.*\\.","",x$taxa))],
        mean_percent=mean(v)*100,first_half_percent=mean(v[1:12500])*100,
        last_half_percent=mean(v[12501:25000])*100,
        q025_percent=quantile(v,.025,names=FALSE)*100,q975_percent=quantile(v,.975,names=FALSE)*100)
    }
  }
}
diagnostics <- bind(diagnostics); chain <- bind(chain); drift <- bind(drift)
occupancy <- bind(occupancy); blocks <- bind(blocks); fitted <- bind(fitted)
reported <- s$tables$convergence
check <- merge(reported,diagnostics[diagnostics$window=="full",],by=c("method","variable"))
stopifnot(nrow(check)==21L,max(abs(check$rhat.x-check$rhat.y))<1e-10)
cat("Cross-checked all21 reported Rhat values against posterior::rhat.\n")

# All retained draws are available for scalar/fitted summaries, but individual
# retained xi states are NOT in this download: no window-specific PIPs claimed.
pip_ranges <- category <- list()
for (method in methods) {
  rr <- getrr(method); pa <- simplify2array(lapply(rr,function(r)r$pips))
  for (j in 1:ncol(x$Y)) {
    v <- pa[,j,]; span <- apply(v,1,function(z)diff(range(z)))
    category[[length(category)+1L]] <- data.frame(method,genus=sub("^.*\\.","",x$taxa[j]),
      category=j,max_spread=max(span),spreads_ge_point1=sum(span>=.1),
      spreads_ge_quarter=sum(span>=.25),spreads_ge_half=sum(span>=.5),
      all_starts_pip_ge_point9=sum(apply(v,1,min)>=.9),
      all_starts_pip_le_point1=sum(apply(v,1,max)<=.1))
    d <- data.frame(method,genus=sub("^.*\\.","",x$taxa[j]),category=j,
      predictor=seq_len(ncol(x$X)),predictor_code=colnames(x$X),pooled_pip=rowMeans(v),
      min_pip=apply(v,1,min),max_pip=apply(v,1,max),spread=span)
    for(k in 1:4)d[[paste0("pip_",starts[k])]]<-v[,k]
    pip_ranges[[length(pip_ranges)+1L]]<-d
  }
}
pip_ranges<-bind(pip_ranges); category<-bind(category)
shared<-pip_ranges[pip_ranges$method==methods[1],c("genus","category","predictor","predictor_code")]
for(k in 1:3) {
  d<-pip_ranges[pip_ranges$method==methods[k],]
  stopifnot(identical(unname(as.matrix(d[,c("category","predictor")])),
    unname(as.matrix(shared[,c("category","predictor")]))))
  shared[[paste0("min_method",k)]]<-d$min_pip
  shared[[paste0("pooled_method",k)]]<-d$pooled_pip
}
shared$min_all12_pip<-pmin(shared$min_method1,shared$min_method2,shared$min_method3)
robust_high<-shared[shared$min_all12_pip>=.9,]
robust_high<-robust_high[order(-robust_high$min_all12_pip),]
labels<-read.table("Real Data Analysis/Data/ffq.labels2.txt",header=TRUE,stringsAsFactors=FALSE)
shared$nutrient<-as.character(labels[shared$predictor_code,1])
shared$max_all12_pip<-apply(simplify2array(lapply(s$records,function(r)r$pips)),c(1,2),max)[
  cbind(shared$predictor,shared$category)]
robust_half<-shared[shared$min_all12_pip>=.5,]
nonprev<-bind(lapply(methods,function(m) {
  d<-category[category$method==m&category$genus!="Prevotella",]
  data.frame(method=m,genera= nrow(d),median_genus_max_pip_spread=median(d$max_spread),
    predictors_spread_ge_quarter=sum(d$spreads_ge_quarter),
    predictors_spread_ge_half=sum(d$spreads_ge_half),max_pip_spread=max(d$max_spread))
}))
cross_method<-shared[,c("genus","category","predictor","predictor_code","nutrient")]
cross_method$pooled_pip_spread<-apply(shared[,paste0("pooled_method",1:3)],1,function(v)diff(range(v)))
cross_method$min_all12_pip<-shared$min_all12_pip
cross_method$max_all12_pip<-shared$max_all12_pip

metrics<-s$tables$chain_metrics
costs<-bind(lapply(methods,function(m) {
  d<-metrics[metrics$method==m,]
  data.frame(method=m,median_wall_hours=median(d$seconds_including_startup)/3600,
    total_allocated_core_hours=sum(d$allocated_core_hours),
    total_worker_cpu_hours=sum(d$worker_cpu_hours_excluding_startup))
}))
tables<-list(diagnostics_by_window=diagnostics,chain_scalar_summary=chain,
  retained_half_drift=drift,retained_quarter_means=blocks,threshold_occupancy=occupancy,
  fitted_subjects=fitted,pip_ranges=pip_ranges,pip_category_disagreement=category,
  robust_high_pip_all12=robust_high,all12_pip_ge_half=robust_half,
  outside_prevotella_stability=nonprev,cross_method_pip_spread=cross_method,descriptive_costs=costs)
for(name in names(tables))write.csv(tables[[name]],file.path(out,paste0(name,".csv")),row.names=FALSE)

# Standalone scientific plots, with block averaging explicitly labelled.
colors<-c("#0072B2","#D55E00","#009E73","#CC79A7")
plot_starts<-c("Null","5%","10%","15%")
pdf(file.path(out,"retained_trace_assessment.pdf"),width=13,height=10,onefile=TRUE)
for(page in list(c("subject3026","subject4014"),c("prev_alpha","prev_size"),c("model_size","log_posterior"))) {
  par(mfrow=c(3,2),mar=c(3.4,4.1,2.5,.6),oma=c(2,0,2,0))
  for(m in 1:3)for(field in page) {
    rr<-getrr(methods[m]); mat<-vapply(rr,function(r)r[[field]],numeric(N))
    if(grepl("subject",field))mat<-mat*100
    yy<-range(unlist(lapply(s$records,function(r)r[[field]])))
    if(grepl("subject",field))yy<-c(0,100)
    if(diff(yy)==0)yy<-yy+c(-1,1)
    plot(range(iteration),yy,type="n",xlab="Iteration (retained sampling only)",ylab=field,
      main=paste(short[m],field,sep=" / "))
    # Raw values faint; 125-draw block means in strong colors.
    for(k in 1:4)lines(iteration,mat[,k],col=adjustcolor(colors[k],alpha.f=.15),lwd=.3)
    groups<-rep(1:(N/125),each=125)
    for(k in 1:4)lines(tapply(iteration,groups,mean),tapply(mat[,k],groups,mean),col=colors[k],lwd=1.2)
    if(grepl("subject",field))abline(h=5,lty=3,col="grey50")
    legend("topright",legend=plot_starts,col=colors,lty=1,bty="n",cex=.7,ncol=2)
  }
  mtext("Faint lines: raw cold draws; strong lines: 125-draw block means. No pre-burn draws included.",
    outer=TRUE,side=1,line=.4,cex=.85)
}
par(mfrow=c(3,2),mar=c(3.8,4.1,2.5,.5))
for(m in 1:3)for(subject in c("3026","4014")) {
  rr<-getrr(methods[m]); bins<-seq(0,100,length.out=51)
  hh<-lapply(rr,function(r)hist(r[[paste0("subject",subject)]]*100,breaks=bins,plot=FALSE))
  plot(c(0,100),c(0,max(vapply(hh,function(h)max(h$density),0))),type="n",
    xlab="Fitted Prevotella (%)",ylab="Empirical density",main=paste(short[m],subject))
  for(k in 1:4)lines(hh[[k]]$mids,hh[[k]]$density,col=colors[k],lwd=1.2)
  legend("topright",legend=plot_starts,col=colors,lty=1,bty="n",cex=.7)
}
par(mfrow=c(3,2),mar=c(3.6,4.1,2.5,.5))
for(m in 1:3)for(subject in c("3026","4014")) {
  d<-blocks[blocks$method==methods[m]&blocks$variable==paste0("subject",subject),]
  plot(c(1,4),c(0,100),type="n",xaxt="n",xlab="Retained quarter",ylab="Mean fitted Prevotella (%)",
    main=paste(short[m],subject))
  axis(1,at=1:4,labels=c("25-31.25k","31.25-37.5k","37.5-43.75k","43.75-50k"),cex.axis=.7)
  for(k in 1:4)with(d[d$start==starts[k],],lines(block,mean*100,type="b",col=colors[k],pch=k))
  legend("topright",legend=plot_starts,col=colors,lty=1,pch=1:4,bty="n",cex=.7)
}
dev.off()
png(file.path(out,"subject4014_retained_windows.png"),width=1800,height=700,res=140)
par(mfrow=c(1,3),mar=c(5,4,3,1))
for(m in 1:3) {
  d<-blocks[blocks$method==methods[m]&blocks$variable=="subject4014",]
  plot(c(1,4),c(0,70),type="n",xaxt="n",xlab="Retained quarter",ylab="Mean fitted Prevotella (%)",main=short[m])
  axis(1,at=1:4)
  for(k in 1:4)with(d[d$start==starts[k],],lines(block,mean*100,type="b",col=colors[k],pch=k,lwd=1.4))
  legend("topright",legend=plot_starts,col=colors,lty=1,pch=1:4,bty="n",cex=.8)
}
dev.off()

fmt<-function(d)capture.output(print(d,row.names=FALSE,digits=5))
report<-c("COMPLETED FRESH50K SLURM RUN: LOCAL READ-ONLY ASSESSMENT",
  paste("Input:",path),paste("posterior package:",packageVersion("posterior")),
  "12 records x25,000 retained cold draws; source hashes match local audited source graph.",
  "All21 reported full-window Rhat values reproduced using posterior::rhat to1e-10.",
  "No checkpoint was downloaded/read or independently revalidated in this assessment.",
  "No MCMC, fitting, compilation, scheduler submission or alteration of inputs was performed.",
  "A posterior mean/quantile here is an empirical chain summary, not proof of a settled posterior estimate.",
  "Threshold5% is the existing fitted-fraction diagnostic region, not an independently identified posterior mode.",
  "Window-specific PIPs/indicator ESS cannot be obtained from this summary: retained xi states are absent.",
  "Within-chain ESS/MCSE can look favorable inside a confined region; cross-chain evidence remains essential.",
  "Costs are descriptive; hardware/load and the states explored differ. No controlled timing or accuracy win claimed.",
  "", "MULTICHAIN DIAGNOSTICS: FULL RETAINED WINDOW",fmt(diagnostics[diagnostics$window=="full",]),
  "", "MULTICHAIN DIAGNOSTICS: LAST HALF",fmt(diagnostics[diagnostics$window=="last_half",]),
  "", "SUBJECT4014: OCCUPANCY AND END-OF-SWEEP CROSSINGS",fmt(occupancy[occupancy$subject=="4014"&occupancy$window=="full",]),
  "", "SUBJECT4014: RETAINED-HALF CHANGES (fraction units)",fmt(drift[drift$variable=="subject4014",]),
  "", "SUBJECT3026: RETAINED-HALF CHANGES (fraction units)",fmt(drift[drift$variable=="subject3026",]),
  "", "PIP DISAGREEMENT: BACTEROIDES AND PREVOTELLA",fmt(category[category$genus%in%c("Bacteroides","Prevotella"),]),
  "", "OUTSIDE PREVOTELLA: WITHIN-METHOD CROSS-START PIP STABILITY",fmt(nonprev),
  "These are descriptive agreement measures, not independent estimates of true PIP accuracy.",
  "", "PIPS >=0.5 IN ALL12 CHAINS (exploratory support, not convergence or causal evidence)",fmt(robust_half),
  "", "HIGH PIPS >=0.9 IN ALL12 CHAINS (exploratory, no convergence claim)",fmt(robust_high),
  "", "DESCRIPTIVE COMPUTATIONAL COSTS",fmt(costs))
stopifnot(identical(before,tools::md5sum(path)),
  identical(s$source_hashes,tools::md5sum(names(s$source_hashes))))
writeLines(report,file.path(out,"report.txt"))
saveRDS(list(tables=tables,input_hash=before,source_hashes=s$source_hashes,
  posterior_version=as.character(packageVersion("posterior")),session=sessionInfo()),file.path(out,"assessment.rds"))
writeLines("COMPLETE",file.path(out,"status.txt"))
cat("ASSESSMENT COMPLETE\n",out,"\n")
cat(paste(report,collapse="\n"),"\n")
