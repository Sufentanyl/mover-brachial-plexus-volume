suppressMessages(library(survival))
suppressMessages(library(jsonlite))
OUT <- 'Sys.getenv('OUT_DIR', '.')'
d <- read.csv(file.path(OUT, 'analysis_571.csv'), stringsAsFactors = FALSE)
d$asa3 <- as.numeric(d$asa3)
covs <- c('vol_le15','vol_ge25','age','asa3','surg_fracture','shoulder','anes_general','inpatient')
dd <- d[d$surgery_type != 'av_fistula' & !is.na(d$asa3), ]
fit <- coxph(as.formula(paste('Surv(time_primary,event_primary) ~', paste(covs, collapse='+'))), data=dd)
s <- summary(fit)
tab <- data.frame(term=rownames(s$conf.int),
  hr=round(s$conf.int[,'exp(coef)'],2), lo=round(s$conf.int[,'lower .95'],2),
  hi=round(s$conf.int[,'upper .95'],2), p=round(s$coefficients[,'Pr(>|z|)'],4),
  stringsAsFactors=FALSE)
print(tab)
get2 <- function(tab, term) {
  r <- tab[tab$term==term, c('hr','lo','hi','p')]
  list(hr=r$hr, lo=r$lo, hi=r$hi, p=r$p)
}
grp_e <- function(dd, ev)
  list(le15=sum(ev[dd$vol_le15==1]), ref=sum(ev[dd$vol_le15==0 & dd$vol_ge25==0]), ge25=sum(ev[dd$vol_ge25==1]))
sj <- fromJSON(file.path(OUT,'sensitivity.json'), simplifyVector=FALSE)
idx <- which(vapply(sj$rows, function(r) r$label=='非AV瘘亚组(其他手术)', logical(1)))
sj$rows[[idx]] <- list(
  label='非AV瘘亚组(其他手术)', n=nrow(dd), events=sum(dd$event_primary),
  epv=round(sum(dd$event_primary)/length(covs),1),
  converged=all(is.finite(coef(fit))),
  note='剔除anes_mac(MAC组18人0事件, 分离); 手术类型仅fracture/other_ue两类, 只保留surg_fracture(other_ue为参照), 避免虚拟变量陷阱',
  events_by_group=grp_e(dd, dd$event_primary),
  le15=get2(tab,'vol_le15'), ge25=get2(tab,'vol_ge25'))
write_json(sj, file.path(OUT,'sensitivity.json'), auto_unbox=TRUE, pretty=TRUE, digits=5)
cat('patched. converged =', all(is.finite(coef(fit))), '\n')
