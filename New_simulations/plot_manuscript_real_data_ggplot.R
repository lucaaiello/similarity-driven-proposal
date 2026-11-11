# Called by make_manuscript_real_data_figures.R after loading saved summaries.
# Uses the same theme and palette as the synthetic-experiment figures.
suppressPackageStartupMessages(library(ggplot2))
source("New_simulations/manuscript_plot_theme.R")
method_colours <- palette[methods]
method_labels <- setNames(method_names,methods)
save_real_data_plot <- function(p,name,width,height) {
  ggsave(file.path(images,paste0(name,".pdf")),p,device=grDevices::cairo_pdf,
    width=width,height=height,units="in",bg="white")
  ggsave(file.path(images,paste0(name,".png")),p,
    width=width,height=height,units="in",dpi=180,bg="white")
}

# Association support: all 30 genera enter the selection criterion.
ar <- pip_ranges
ar$method <- factor(ar$method,levels=methods)
ar$y <- nrow(associations)+1-ar$association+c(.18,0,-.18)[as.integer(ar$method)]
p <- ggplot(ar,aes(colour=method)) +
  geom_segment(aes(x=min_pip,xend=max_pip,y=y,yend=y),linewidth=.5) +
  geom_point(aes(x=mean_pip,y=y),size=2) +
  scale_colour_manual(values=method_colours,labels=method_labels) +
  scale_y_continuous(breaks=seq_len(nrow(associations)),
    labels=rev(paste(associations$genus,associations$nutrient,sep=": ")),
    expand=expansion(add=.5)) +
  scale_x_continuous(breaks=seq(.5,1,.1)) +
  coord_cartesian(xlim=c(.45,1.01)) +
  labs(x="Empirical PIP: range across four starts",y=NULL) + theme_manuscript() +
  guides(colour=guide_legend(nrow=1))
save_real_data_plot(p,"real_data_associations",7.6,4)

# Paired descriptive summaries from the same full-model draws.
cs <- category_stats
cs$method <- factor(cs$method,levels=methods)
stopifnot(all(table(cs$method)==30L))
paired <- rbind(transform(cs,scope="All 30 genera"),
  transform(cs[cs$genus!="Prevotella",],scope="Without Prevotella"))
paired$scope <- factor(paired$scope,levels=c("All 30 genera","Without Prevotella"))
stopifnot(all(table(paired$method,paired$scope)[,1]==30L),
  all(table(paired$method,paired$scope)[,2]==29L))
p <- ggplot(paired,aes(x=method,y=max_spread)) +
  geom_boxplot(aes(fill=method,alpha=scope,group=interaction(method,scope)),
    position=position_dodge(width=.7),width=.6,linewidth=.4,coef=Inf,outlier.shape=NA) +
  scale_fill_manual(values=method_colours,guide="none") +
  scale_alpha_manual(values=c(.25,.65),
    guide=guide_legend(override.aes=list(fill="#666666",colour="#333333"))) +
  scale_x_discrete(labels=method_labels) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2)) +
  labs(x=NULL,y="Maximum PIP spread within genus") + theme_manuscript()
save_real_data_plot(p,"real_data_spreads",7.6,3.3)

# Category-by-category view also includes all 30 genera.
cs$y <- 31-cs$category+c(.18,0,-.18)[as.integer(cs$method)]
p <- ggplot(cs,aes(max_spread,y,colour=method)) +
  geom_point(size=1.8) +
  scale_colour_manual(values=method_colours,labels=method_labels) +
  scale_y_continuous(breaks=1:30,labels=rev(sub("^.*\\.","",x$taxa)),
    expand=expansion(add=.6)) +
  scale_x_continuous(breaks=seq(0,1,.2),limits=c(0,1)) +
  labs(x="Maximum cross-start PIP difference within genus",y=NULL) +
  theme_manuscript() + guides(colour=guide_legend(nrow=1))
save_real_data_plot(p,"real_data_pip_stability",7.6,8)

# Retained parameter traces: same observations and 125-draw block means.
traces <- list(); blocks <- list()
fields <- c("prev_alpha","prev_size")
field_labels <- c(prev_alpha="Prevotella intercept",prev_size="Active Prevotella predictors")
for(m in seq_along(methods))for(k in 1:4)for(field in fields) {
  values <- rr(methods[m])[[k]][[field]]
  id <- length(traces)+1L
  traces[[id]] <- data.frame(iteration=it,value=values,method=methods[m],
    start=start_names[k],quantity=field)
  group <- rep(1:200,each=125)
  blocks[[id]] <- data.frame(iteration=as.numeric(tapply(it,group,mean)),
    value=as.numeric(tapply(values,group,mean)),method=methods[m],start=start_names[k],quantity=field)
}
traces <- bind(traces); blocks <- bind(blocks)
for(name in c("traces","blocks")) {
  d <- get(name)
  d$method <- factor(d$method,levels=methods,labels=method_names)
  d$quantity <- factor(d$quantity,levels=fields,labels=field_labels[fields])
  d$start <- factor(d$start,levels=start_names)
  assign(name,d)
}
# Match vertical ranges across methods for each monitored quantity.
limits <- bind(lapply(levels(traces$quantity),function(q) {
  expand.grid(method=levels(traces$method),quantity=q,
    value=range(traces$value[traces$quantity==q]),stringsAsFactors=FALSE)
}))
limits$method <- factor(limits$method,levels=levels(traces$method))
limits$quantity <- factor(limits$quantity,levels=levels(traces$quantity))
p <- ggplot(traces,aes(iteration,value,colour=start)) +
  geom_line(linewidth=.15,alpha=.12) +
  geom_line(data=blocks,linewidth=.4) +
  geom_blank(data=limits,aes(x=25001,y=value),inherit.aes=FALSE) +
  facet_wrap(vars(method,quantity),ncol=2,scales="free_y") +
  scale_colour_manual(values=setNames(start_colors,start_names)) +
  scale_x_continuous(breaks=c(25000,35000,50000),labels=c("25,000","35,000","50,000")) +
  labs(x="Retained iteration",y=NULL) + theme_manuscript() +
  theme(strip.text=element_text(size=9),panel.spacing=grid::unit(.8,"lines"))
save_real_data_plot(p,"real_data_retained_parameters",9,7)
