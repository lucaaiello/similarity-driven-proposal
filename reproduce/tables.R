# Numerical table exports computed from saved experiment summaries/draws.
a<-commandArgs(TRUE)
opt<-function(k,d="") {x<-grep(paste0("^--",k,"="),a,value=TRUE);if(length(x))sub(paste0("^--",k,"="),"",x[1])else d}
out<-opt("output");stopifnot(nzchar(out))
dest<-file.path(out,"tables");dir.create(dest,recursive=TRUE,showWarnings=FALSE)
g<-"New_simulations/outputs/matched_gaussian_20261007_163020"
rd<-function(p)read.csv(p,check.names=FALSE)
avg<-function(d,by,vars)aggregate(d[vars],d[by],mean,na.rm=TRUE)
sm<-function(d,by,vars) {
  pieces<-split(d,interaction(d[by],drop=TRUE,lex.order=TRUE))
  z<-lapply(pieces,function(v) {
    r<-v[1,by,drop=FALSE]
    for(n in vars) {r[[paste0(n,"_mean")]]<-mean(v[[n]],na.rm=TRUE)
      r[[paste0(n,"_se")]]<-sd(v[[n]],na.rm=TRUE)/sqrt(sum(is.finite(v[[n]])))}
    r
  });z<-do.call(rbind,z);rownames(z)<-NULL;z
}
export<-function(id,d)write.csv(d,file.path(dest,paste0(id,".csv")),row.names=FALSE)
sim<-function(f)rd(file.path(dest,"simulation",f))
# S3, S4 and S6: existing matched chain-averaged aggregation.
export("S3_exact_gaussian",sim("table_exact_reference.csv"))
export("S4_scaling",sim("table_highdim_p1000.csv"))
export("S6_correlated_moves",sim("table_correlated_local_moves.csv"))
high<-rd(file.path(g,"scaling/chain_metrics.csv"));high<-high[high$dimension==1000,]
hvars<-c("final_lambda","mean_proposal_effective_neighbors","acceptance_rate")
h<-sm(avg(high,c("dataset","method"),hvars),"method",hvars)
p<-rd(file.path(g,"scaling/pip_estimates.csv"));p<-p[p$dimension==1000,]
range1<-function(x)diff(range(x))
pr<-aggregate(pip~dataset+method+variable,p,range1)
pr<-avg(pr,c("dataset","method"),"pip");pr<-avg(pr,"method","pip")
names(pr)[2]<-"mean_predictor_pip_range"
export("S5_scaling_diagnostics",merge(h,pr,by="method"))
# S7 and S8: keep adaptive score results distinct from the fixed-lambda grid.
e<-rd(file.path(g,"exact_score/chain_metrics.csv"));e<-e[e$regime=="online",]
e$robust_ess_second<-pmin(e$ess_model_size_per_second,e$ess_log_target_per_second)
ev<-c("lambda_final","mean_effective_neighbors","robust_ess_second","mean_pip_error","total_variation")
es<-sm(avg(e,c("dataset","score"),ev),"score",ev)
al<-rd(file.path(g,"exact_score/score_alignment.csv"))
al<-sm(al,"score","posterior_weighted_rank_correlation")
export("S7_exact_score",merge(es,al,by="score"))
corr<-rd(file.path(g,"correlated_score/chain_metrics.csv"))
corr$robust_ess_second<-pmin(corr$ess_model_size_per_second,corr$ess_log_target_per_second)
cv<-c("lambda_final","mean_effective_neighbors","acceptance","robust_ess_second","rhat_model_size","rhat_log_target","mean_active_pip")
cs<-sm(avg(corr,c("dataset","score"),cv),"score",cv)
cp<-rd(file.path(g,"correlated_score/pip_estimates.csv"))
cp<-aggregate(pip~dataset+score+variable,cp,range1)
cp<-avg(avg(cp,c("dataset","score"),"pip"),"score","pip")
names(cp)[2]<-"mean_predictor_pip_range"
export("S8_correlated_score",merge(cs,cp,by="score"))
# Main2/S12/S13: average chains within dataset before MCSE over datasets.
dm<-rd("New_simulations/outputs/dm_gagnon_comparison/combined_chain_metrics.csv")
dm$robust_ess_second<-pmin(dm$ess_model_size_per_second,dm$ess_log_posterior_per_second)
dv<-c("seconds","acceptance_flip","acceptance_local","robust_ess_second","rhat_model_size","rhat_log_posterior",
 "pip_rmse","true_positive_rate_05","false_positive_rate_05","active_coefficient_rmse","active_coefficient_coverage")
ds<-sm(avg(dm,c("dataset","method"),dv),"method",dv)
export("main2_dm_computation",ds[c("method","seconds_mean","seconds_se","robust_ess_second_mean","robust_ess_second_se","rhat_log_posterior_mean","rhat_log_posterior_se")])
export("S12_dm_recovery",ds[c("method",grep("pip_rmse|true_positive|false_positive|active_coefficient",names(ds),value=TRUE))])
export("S13_dm_computation",ds[c("method",grep("seconds_|acceptance_|robust_|rhat_",names(ds),value=TRUE))])
export("S9_dm_dataset_diagnostics",rd(file.path(dest,"dm_dataset_diagnostics/diagnostic_summary.csv")))
# S10: recompute aggregates from parameter-level saved diagnostics.
pd<-rd(opt("parameters"));pd<-pd[pd$quantity!="coefficient"|pd$generating_active,]
safe<-function(x,fun)if(any(is.finite(x)))fun(x[is.finite(x)])else NA_real_
parts<-split(pd,interaction(pd$score_steps,pd$quantity,drop=TRUE))
ps<-do.call(rbind,lapply(parts,function(z)data.frame(score_steps=z$score_steps[1],quantity=z$quantity[1],
 informative=sum(is.finite(z$rhat)),total=nrow(z),rhat_above_1_01=sum(z$rhat>1.01,na.rm=TRUE),
 max_rhat=safe(z$rhat,max),min_bulk_ess=safe(z$bulk_ess,min),max_mcse_mean=safe(z$mcse_mean,max))))
export("S10_parameter_diagnostics",ps[order(ps$score_steps,ps$quantity),])
# S11: chain means within each of the five datasets, then means/MCSEs.
depth<-"New_simulations/outputs/dm_fitting_depth_20261008_104824"
ss<-lapply(c(1,5,10),function(b) {
  z<-rd(file.path(depth,paste0("depth_",b),"chain_metrics.csv"))
  z$robust_ess_second<-pmin(z$ess_model_size_per_second,z$ess_log_posterior_per_second)
  v<-c("seconds","robust_ess_second","pip_rmse","rhat_model_size","rhat_log_posterior")
  d<-avg(z,"dataset",v);r<-sm(transform(d,budget=b),"budget",v)
  r$minutes_per_chain<-r$seconds_mean/60;r$max_rhat_log_posterior<-max(d$rhat_log_posterior)
  r$datasets_above_1_01<-sum(d$rhat_model_size>1.01|d$rhat_log_posterior>1.01)
  r$datasets<-nrow(d);r
});export("S11_fitting_depth",do.call(rbind,ss))
# Main3/Main4/S14/S15: generated directly from the completed cold chains.
app<-file.path(dest,"application")
export("main3_application_comparison",rd(file.path(app,"method_comparison.csv")))
export("main4_application_associations",rd(file.path(app,"association_method_ranges.csv")))
s<-readRDS("New_simulations/outputs/dm_apt50k_fresh_slurm_v1/summary/summary.rds")
chains<-s$tables$chain_metrics
chains$hours<-chains$seconds_including_startup/3600
chains$cold_exchange_acceptance<-chains$cold_swaps_post/25000
chains$mean_model_size<-vapply(seq_len(nrow(chains)),function(i) {
  r<-s$records[vapply(s$records,function(r)r$method==chains$method[i]&&
    r$start==match(chains$start[i],s$plan$start_names),TRUE)]
  stopifnot(length(r)==1L);mean(r[[1]]$model_size)
},1.)
export("S14_application_chains",chains)
export("S15_application_scalars",rd(file.path(app,"scalar_diagnostics.csv")))
# Static design/method tables are definitions, not estimates from chains.
for(doc in c("main_r1","supplements_r1")) {
  text<-paste(readLines(file.path("JCGS_submission",paste0(doc,".tex")),warn=FALSE),collapse="\n")
  blocks<-regmatches(text,gregexpr("(?s)\\\\begin\\{table\\}.*?\\\\end\\{table\\}",text,perl=TRUE))[[1]]
  for(b in blocks) {
    label<-regmatches(b,regexec("\\\\label\\{([^}]+)\\}",b))[[1]]
    if(length(label)>1&&label[2]%in%c("tab:DM_methods","tab:simulation_designs","tab:simulation_methods"))
      writeLines(b,file.path(dest,paste0(gsub(":","_",label[2]),".tex")))
  }
}
cat("All computational table summaries exported; static definition tables extracted.\n")
