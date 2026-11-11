#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Rcpp))
options(stringsAsFactors=FALSE,warn=1)
args<-commandArgs(trailingOnly=TRUE)
arg_value<-function(name,default){z<-grep(paste0("^--",name,"="),args,value=TRUE);if(length(z))sub(paste0("^--",name,"="),"",z[[1]])else default}
mode<-arg_value("mode","quick");if(!mode%in%c("quick","full"))stop("--mode must be quick or full")
n<-as.integer(arg_value("n","200"));p<-as.integer(arg_value("p","200"));blocks<-as.integer(arg_value("blocks","10"));block_size<-as.integer(arg_value("block-size","5"))
datasets<-as.integer(arg_value("datasets",if(mode=="quick")"1" else "30"));chains<-4L
runs<-as.integer(arg_value("runs",if(mode=="quick")"2000" else "50000"));burn<-as.integer(arg_value("burn",if(mode=="quick")"500" else "10000"))
rho<-as.numeric(arg_value("rho","0.95"));signal_size<-as.numeric(arg_value("signal-size","0.8"));graph_threshold<-as.numeric(arg_value("graph-threshold","0.7"))
out_dir<-arg_value("output",file.path("New_simulations","outputs",paste0("correlated_surrogate_local_move_",mode,"_matched_v2")))
audit_file<-arg_value("audit",file.path("New_simulations","outputs","local_move_rank_one_equivalence","summary.csv"))
if(p<blocks*block_size||burn>=runs||datasets<1||rho<=0||rho>=1)stop("Invalid experiment configuration")
if(!file.exists(audit_file)||!isTRUE(read.csv(audit_file)$passed[[1]]))stop("Run the local-move rank-one equivalence audit first")
source(file.path("New_simulations","gaussian_revision_support.R"));revision_initialize(out_dir,require_audit=TRUE);revision_require_audit(dirname(audit_file))

Rcpp::sourceCpp(file.path("Simulations","sampler_F_stat.cpp"),rebuild=TRUE,verbose=FALSE);f_sampler<-get("var_sel_distance_cpp",.GlobalEnv)
Rcpp::sourceCpp(file.path("Simulations","sampler_LR.cpp"),rebuild=TRUE,verbose=FALSE);bic_sampler<-get("var_sel_distance_cpp",.GlobalEnv)
methods<-expand.grid(score=c("f","bic"),move=c("flip_only","uniform_swap","informed_swap"),stringsAsFactors=FALSE)
methods$method<-paste(methods$score,methods$move,sep="_")

ess_basic<-function(x){if(length(x)<20||var(x)==0)return(NA_real_);z<-as.numeric(acf(x,lag.max=min(length(x)-1,1000),plot=FALSE)$acf)[-1];m<-length(z)-length(z)%%2;if(m<2)return(length(x));q<-z[seq(1,m,2)]+z[seq(2,m,2)];k<-match(TRUE,q<=0,nomatch=length(q)+1);keep<-if(k==1)numeric()else q[seq_len(k-1)];min(length(x),length(x)/(1+2*sum(keep)))}
rank_split_rhat<-function(chains){L<-min(vapply(chains,length,integer(1)));h<-floor(L/2);if(h<20)return(NA_real_);s<-unlist(lapply(chains,function(x)list(x[seq_len(h)],x[(L-h+1):L])),recursive=FALSE);z<-qnorm((rank(unlist(s),ties.method="average")-.5)/(length(s)*h));M<-matrix(z,nrow=h);W<-mean(apply(M,2,var));B<-h*var(colMeans(M));sqrt(((h-1)/h*W+B/h)/W)}
simulate_data<-function(seed){set.seed(seed);X<-matrix(rnorm(n*p),n,p);for(b in seq_len(blocks)){idx<-((b-1)*block_size+1):(b*block_size);latent<-rnorm(n);X[,idx]<-sqrt(rho)*latent+sqrt(1-rho)*matrix(rnorm(n*block_size),n)};X<-unclass(scale(X));active<-1+(0:(blocks-1))*block_size;beta<-numeric(p);beta[active]<-rep(c(signal_size,-signal_size),length.out=blocks);list(X=X,y=as.numeric(.4+X%*%beta+rnorm(n)),active=active)}
make_graph<-function(X){C<-abs(cor(X));diag(C)<-0;lapply(seq_len(ncol(X)),function(j)which(C[j,]>=graph_threshold))}
make_starts<-function(active){null<-integer(p);truth<-null;truth[active]<-1;surrogate<-null;surrogate[active+1]<-1;set.seed(91);random<-null;random[sample.int(p,blocks)]<-1;list(null,truth,surrogate,random)}

chain_rows<-list();pip_rows<-list();block_rows<-list();adapt_rows<-list();id<-0L;bid<-0L;aid<-0L
for(dataset in seq_len(datasets)){
  dat<-simulate_data(20270000L+dataset);graph_time<-system.time(neigh<-make_graph(dat$X))[["elapsed"]];starts<-make_starts(dat$active);design<-cbind(1,dat$X)
  true_edge<-function(i,j) i<=blocks*block_size&&j<=blocks*block_size&&((i-1)%/%block_size)==((j-1)%/%block_size)
  edge_pairs<-which(upper.tri(matrix(FALSE,p,p)),arr.ind=TRUE);edge_selected<-vapply(seq_len(nrow(edge_pairs)),function(k)edge_pairs[k,2]%in%neigh[[edge_pairs[k,1]]],logical(1));edge_truth<-vapply(seq_len(nrow(edge_pairs)),function(k)true_edge(edge_pairs[k,1],edge_pairs[k,2]),logical(1));graph_recall<-sum(edge_selected&edge_truth)/sum(edge_truth);graph_precision<-if(any(edge_selected))sum(edge_selected&edge_truth)/sum(edge_selected)else NA_real_
  for(mi in seq_len(nrow(methods))){spec<-methods[mi,];traces<-vector("list",chains);block_traces<-vector("list",chains)
    for(chain in seq_len(chains)){set.seed(20271000L+dataset*100+mi*10+chain);sampler<-if(spec$score=="f")f_sampler else bic_sampler
      call<-list(y=dat$y,X=design,neigh=neigh,gamma=starts[[chain]],phi=.05,a=2,d=1,tau=1,a_phi=1,b_phi=3,runs=runs,burn=burn,thin=1,start_adapt=100,end_adapt=burn,lambda_bd=.7,lambda_move=1,trgt_acpt=.234,adapt=TRUE,move_step=spec$move!="flip_only",adaptation_mode="diversity",balanced_directions=TRUE,uniform_mixture=.02,target_effective_fraction=.025,adaptation_gain=.5,adaptation_decay=.6,adaptation_offset=100,proposal_engine="rank_one",move_weight_mode=if(spec$move=="uniform_swap")"uniform" else "informed")
      if(spec$score=="bic")call$score_type<-"bic";tm<-system.time(fit<-do.call(sampler,call))[["elapsed"]];G<-fit$gamma;ms<-rowSums(G);traces[[chain]]<-ms;block_traces[[chain]]<-lapply(seq_len(blocks),function(b)G[,dat$active[b]])
      ham<-if(nrow(G)>1)rowSums(abs(G[-1,,drop=FALSE]-G[-nrow(G),,drop=FALSE])) else NA
      id<-id+1L;chain_rows[[id]]<-data.frame(dataset=dataset,chain=chain,method=spec$method,score=spec$score,move=spec$move,seconds=tm,graph_seconds=graph_time,graph_edge_recall=graph_recall,graph_edge_precision=graph_precision,iterations=runs,retained=nrow(G),acceptance_flip=fit$acc_rate_bd,accepted_swaps_per_iteration=fit$acc_rate_move,swap_attempt_rate=fit$move_attempts/runs,conditional_swap_acceptance=fit$conditional_move_acceptance,final_lambda=tail(fit$lambda_chain,1),ess_model_size=ess_basic(ms),ess_log_target=ess_basic(fit$loglik),ess_model_size_per_second=ess_basic(ms)/tm,mean_hamming=mean(ham,na.rm=TRUE),multi_component_jump_rate=mean(ham>=2,na.rm=TRUE))
      keep<-unique(round(seq(1,runs,length.out=min(250L,runs))));aid<-aid+1L;adapt_rows[[aid]]<-data.frame(dataset=dataset,chain=chain,method=spec$method,iteration=keep,lambda=fit$lambda_chain[keep],acceptance_rate=fit$acc_rate_bd_chain[keep],effective_neighbors=fit$proposal_effective_neighbors[keep])
      pp<-colMeans(G);pip_rows[[id]]<-data.frame(dataset=dataset,chain=chain,method=spec$method,variable=seq_len(p),active=seq_len(p)%in%dat$active,pip=pp)
      for(b in seq_len(blocks)){idx<-((b-1)*block_size+1):(b*block_size);B<-G[,idx,drop=FALSE];representative<-apply(B,1,function(z){w<-which(z>0);if(length(w)==1)w else 0L});bid<-bid+1L;block_rows[[bid]]<-data.frame(dataset=dataset,chain=chain,method=spec$method,block=b,true_pip=mean(B[,1]),surrogate_pip_sum=sum(colMeans(B[,-1,drop=FALSE])),any_inclusion=mean(rowSums(B)>0),multiple_inclusion=mean(rowSums(B)>1),switches_per_1000=1000*mean(diff(representative)!=0),substitutions_per_1000=1000*mean(diff(representative)!=0 & head(representative,-1)>0 & tail(representative,-1)>0))}
      message(sprintf("dataset=%d/%d method=%s chain=%d/%d seconds=%.1f",dataset,datasets,spec$method,chain,chains,tm))
    }
    rhat_model<-rank_split_rhat(traces);rhat_blocks<-vapply(seq_len(blocks),function(b)rank_split_rhat(lapply(block_traces,function(x)x[[b]])),numeric(1));for(k in seq_len(chains)){chain_rows[[id-chains+k]]$rhat_model_size<-rhat_model;chain_rows[[id-chains+k]]$mean_rhat_true_indicator<-mean(rhat_blocks,na.rm=TRUE)}
  }
}
raw<-do.call(rbind,chain_rows);pips<-do.call(rbind,pip_rows);block_metrics<-do.call(rbind,block_rows);adaptation<-do.call(rbind,adapt_rows)
num<-setdiff(names(raw)[vapply(raw,is.numeric,logical(1))],c("dataset","chain"));summary<-aggregate(raw[num],raw[c("method","score","move")],mean,na.rm=TRUE)
write.csv(raw,file.path(out_dir,"chain_metrics.csv"),row.names=FALSE);write.csv(summary,file.path(out_dir,"method_summary.csv"),row.names=FALSE);write.csv(pips,file.path(out_dir,"pip_estimates.csv"),row.names=FALSE);write.csv(block_metrics,file.path(out_dir,"block_metrics.csv"),row.names=FALSE);write.csv(adaptation,file.path(out_dir,"adaptive_trajectories_thinned.csv"),row.names=FALSE)
writeLines(c(paste("mode:",mode),paste("n:",n),paste("p:",p),paste("datasets:",datasets),paste("chains:",chains),paste("runs:",runs),paste("burn:",burn),paste("blocks:",blocks),paste("block_size:",block_size),paste("rho:",rho),paste("graph_threshold:",graph_threshold),"engine: production rank-one flip and swap scoring","adaptation: genuine online diminishing diversity adaptation throughout the run, matching the production experiments","timing: sampler time reported separately from one-off graph construction time"),file.path(out_dir,"manifest.txt"));writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"));saveRDS(list(raw=raw,summary=summary,pips=pips,block_metrics=block_metrics,adaptation=adaptation),file.path(out_dir,"results.rds"),compress="xz");print(summary,row.names=FALSE)
