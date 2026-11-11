# Compare recomputed numbers with the explicitly typeset tables, at the
# precision printed in each cell. Never edits manuscript content.
a<-commandArgs(TRUE);out<-sub("^--output=","",a[grepl("^--output=",a)])[1]
stopifnot(length(out)==1L,!is.na(out))
tab<-file.path(out,"tables");read<-function(f)read.csv(file.path(tab,f),check.names=FALSE)
methods<-c("uniform_flip","optimized_similarity_fixed","online_similarity_flip","online_similarity_local","gagnon_laplace_rj")
gauss<-c("uniform","f_fixed","f_online","bic_fixed","bic_online","locally_balanced","bic_locally_balanced")
app<-c("online_similarity_flip","online_similarity_local","uniform_flip")
ord<-function(d,key,ids)d[match(ids,d[[key]]),,drop=FALSE]
pair<-function(z,n,m=1)c(z[[paste0(n,"_mean")]],z[[paste0(n,"_se")]])*m
checks<-list()
check<-function(id,doc,label,expected,skip=1L) {
  text<-paste(readLines(file.path("JCGS_submission",paste0(doc,".tex")),warn=FALSE),collapse="\n")
  blocks<-regmatches(text,gregexpr("(?s)\\\\begin\\{table\\}.*?\\\\end\\{table\\}",text,perl=TRUE))[[1]]
  block<-blocks[grepl(paste0("\\label{",label,"}"),blocks,fixed=TRUE)]
  stopifnot(length(block)==1L)
  body<-strsplit(block,"\\midrule",fixed=TRUE)[[1]][2]
  body<-strsplit(body,"\\bottomrule",fixed=TRUE)[[1]][1]
  body<-gsub("(?m)(?<!\\\\)%.*$","",body,perl=TRUE)
  rows<-strsplit(body,"\\\\",fixed=TRUE)[[1]]
  rows<-rows[grepl("&",rows,fixed=TRUE)]
  if(length(rows)!=length(expected))stop("Row count mismatch in ",id)
  count<-0L
  for(i in seq_along(rows)) {
    cells<-strsplit(rows[i],"&",fixed=TRUE)[[1]]
    s<-paste(cells[-seq_len(skip)],collapse=" ")
    s<-gsub("{,}","",s,fixed=TRUE);s<-gsub("\\,","",s,fixed=TRUE)
    tokens<-regmatches(s,gregexpr("[0-9]+(?:[.][0-9]+)?",s,perl=TRUE))[[1]]
    x<-as.numeric(tokens);y<-expected[[i]]
    if(length(x)!=length(y))stop("Cell count mismatch in ",id," row ",i,": ",length(x)," vs ",length(y))
    decimal<-ifelse(grepl(".",tokens,fixed=TRUE),nchar(sub("^.*[.]","",tokens)),0L)
    ok<-abs(x-y)<=.500001*10^(-decimal)+1e-10
    if(!all(ok))stop("Printed-number mismatch in ",id," row ",i,": ",paste(tokens[!ok],collapse=","),
      " vs ",paste(signif(y[!ok],9),collapse=","))
    count<-count+length(x)
  }
  checks[[length(checks)+1L]]<<-data.frame(table=id,rows=length(rows),numeric_entries=count,passed=TRUE)
}
rows<-function(d,f)lapply(seq_len(nrow(d)),function(i)f(d[i,,drop=FALSE]))
d<-ord(read("S3_exact_gaussian.csv"),"method",gauss)
check("S3","supplements_r1","tab:exact_gaussian_results",rows(d,function(z)c(
 pair(z,"mean_abs_pip_error",100),pair(z,"total_variation_distance",100),pair(z,"robust_ess_second"),pair(z,"acceptance_rate"),
 if(z$method%in%c("f_fixed","f_online","bic_fixed","bic_online"))c(z$final_lambda_mean,z$mean_last_quarter_effective_neighbors_mean)else numeric())))
d<-ord(read("S4_scaling.csv"),"method",gauss)
check("S4","supplements_r1","tab:highdim_scaling",rows(d,function(z)c(z$empirical_time_exponent,pair(z,"iterations_per_second"),pair(z,"robust_ess_second"),z$top_signal_recall_mean)))
d<-ord(read("S5_scaling_diagnostics.csv"),"method",gauss)
check("S5","supplements_r1","tab:highdim_proposal_diagnostics",rows(d,function(z)c(
 if(z$method%in%c("f_fixed","bic_fixed"))z$final_lambda_mean else if(z$method%in%c("f_online","bic_online"))pair(z,"final_lambda") else numeric(),
 if(z$method%in%c("f_fixed","f_online","bic_fixed","bic_online"))pair(z,"mean_proposal_effective_neighbors")else numeric(),
 # S5 originally rounded guided-method acceptance summaries to three
 # decimals before adding a trailing zero; preserve that presentation rule.
 if(z$method=="uniform")pair(z,"acceptance_rate")else round(pair(z,"acceptance_rate"),3),
 z$mean_predictor_pip_range)))
d<-ord(read("S6_correlated_moves.csv"),"method",c("f_flip_only","f_uniform_swap","f_informed_swap","bic_flip_only","bic_uniform_swap","bic_informed_swap"))
check("S6","supplements_r1","tab:correlated_local_moves",rows(d,function(z)c(pair(z,"switches_per_1000"),
 if(grepl("flip_only",z$method))z$multi_component_jump_rate_mean*100 else pair(z,"multi_component_jump_rate",100),
 pair(z,"robust_ess_second"),pair(z,"true_pip"),pair(z,"surrogate_pip_sum"))))
d<-ord(read("S7_exact_score.csv"),"score",c("f","bic","marginal"))
check("S7","supplements_r1","tab:exact_score_sensitivity",rows(d,function(z)c(pair(z,"posterior_weighted_rank_correlation"),round(pair(z,"lambda_final"),3),pair(z,"mean_effective_neighbors"),pair(z,"robust_ess_second"),pair(z,"mean_pip_error",100),pair(z,"total_variation",100))))
d<-ord(read("S8_correlated_score.csv"),"score",c("f","bic","marginal"))
check("S8","supplements_r1","tab:correlated_score_sensitivity",rows(d,function(z)c(z$lambda_final_mean,z$mean_effective_neighbors_mean,z$acceptance_mean,pair(z,"robust_ess_second"),z$rhat_model_size_mean,z$rhat_log_target_mean,z$mean_predictor_pip_range,z$mean_active_pip_mean)))
d<-ord(read("main2_dm_computation.csv"),"method",methods)
check("Main2","main_r1","tab:DM_computation",rows(d,function(z)c(pair(z,"seconds"),pair(z,"robust_ess_second"),pair(z,"rhat_log_posterior"))))
d<-ord(read("S9_dm_dataset_diagnostics.csv"),"method",methods)
check("S9","supplements_r1","tab:DM_dataset_diagnostics",rows(d,function(z)unlist(z[c("size_median","size_q75","size_max","log_median","log_q75","log_max","n_above_1_01","n_above_1_05","pip_range_median","pip_range_max")],use.names=FALSE)))
d<-read("S10_parameter_diagnostics.csv")
check("S10","supplements_r1","tab:DM_parameter_diagnostics",rows(d,function(z)c(z$rhat_above_1_01,z$informative,z$max_rhat,z$min_bulk_ess,z$max_mcse_mean)),skip=2)
d<-read("S11_fitting_depth.csv")
check("S11","supplements_r1","tab:DM_implementation_sensitivity",rows(d,function(z)c(z$minutes_per_chain,pair(z,"robust_ess_second"),z$pip_rmse_mean,z$max_rhat_log_posterior,z$datasets_above_1_01,z$datasets)))
d<-ord(read("S12_dm_recovery.csv"),"method",methods)
check("S12","supplements_r1","tab:DM_accuracy",rows(d,function(z)c(pair(z,"pip_rmse"),pair(z,"true_positive_rate_05"),pair(z,"false_positive_rate_05",100),pair(z,"active_coefficient_rmse"),pair(z,"active_coefficient_coverage"))))
d<-ord(read("S13_dm_computation.csv"),"method",methods)
check("S13","supplements_r1","tab:DM_computation_full",rows(d,function(z)c(pair(z,"seconds"),pair(z,"acceptance_flip"),if(z$method=="online_similarity_local")pair(z,"acceptance_local")else numeric(),pair(z,"robust_ess_second"),pair(z,"rhat_model_size"),pair(z,"rhat_log_posterior"))))
d<-ord(read("main3_application_comparison.csv"),"method",app)
check("Main3","main_r1","tab:real_data_comparison",rows(d,function(z)c(z$model_size_rhat,z$other_model_size_rhat,z$all_genera_median_max_spread,z$other_genera_median_max_spread)))
d<-read("main4_application_associations.csv")
check("Main4","main_r1","tab:vs_taxa_covariate",lapply(sort(unique(d$association)),function(i) {
 z<-ord(d[d$association==i,],"method",app);as.vector(t(z[c("min_pip","max_pip")]))
}),skip=2)
d<-read("S14_application_chains.csv");d<-d[order(match(d$method,app),match(d$start,c("null","moderate_5pct","moderate_10pct","dense_15pct"))),]
check("S14","supplements_r1","tab:real_data_chain_diagnostics",rows(d,function(z)c(z$hours,z$mean_model_size,z$cold_exchange_acceptance)),skip=2)
d<-read("S15_application_scalars.csv");d<-d[d$variable!="other_model_size",]
d<-d[order(match(d$method,app),match(d$variable,c("model_size","prev_size","prev_alpha","log_posterior","positive_subject_mean"))),]
check("S15","supplements_r1","tab:real_data_scalar_diagnostics",rows(d,function(z)c(z$rhat,z$bulk_ess,z$tail_ess)),skip=2)
write.csv(do.call(rbind,checks),file.path(out,"table_verification.csv"),row.names=FALSE)
cat("All 16 computational tables agree with the manuscript/supplement at printed precision.\n")
