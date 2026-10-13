suppressMessages(library(survival))
suppressMessages(library(jsonlite))
OUT <- 'Sys.getenv('OUT_DIR', '.')'
d <- read.csv(file.path(OUT, 'analysis_571.csv'), stringsAsFactors = FALSE)
d$asa3 <- as.numeric(d$asa3)
COVS10 <- c('vol_le15','vol_ge25','age','asa3','surg_fracture','surg_other',
            'shoulder','anes_general','anes_mac','inpatient')
mkform <- function(tvar, evar, covs)
  as.formula(paste0('Surv(', tvar, ',', evar, ') ~ ', paste(covs, collapse = '+')))
cox_tab <- function(fit) {
  s <- summary(fit)
  data.frame(term = rownames(s$conf.int),
             hr = round(s$conf.int[, 'exp(coef)'], 2),
             lo = round(s$conf.int[, 'lower .95'], 2),
             hi = round(s$conf.int[, 'upper .95'], 2),
             p  = round(s$coefficients[, 'Pr(>|z|)'], 4),
             stringsAsFactors = FALSE)
}
get2 <- function(tab, term) {
  r <- tab[tab$term == term, c('hr','lo','hi','p')]
  list(hr = r$hr, lo = r$lo, hi = r$hi, p = r$p)
}
grp_e <- function(dd, ev, a = 'vol_le15', b = 'vol_ge25')
  list(le15 = sum(ev[dd[[a]] == 1]), ref = sum(ev[dd[[a]] == 0 & dd[[b]] == 0]),
       ge25 = sum(ev[dd[[b]] == 1]))

fits <- list(
  list(label = 'AV瘘亚组',
       dat = d[d$surgery_type == 'av_fistula' & !is.na(d$asa3), ],
       covs = c('vol_le15','vol_ge25','age','asa3','anes_general','anes_mac','inpatient'),
       note = '剔除恒定协变量: surg_fracture/surg_other(手术类型恒定)、shoulder(亚组内全为0)'),
  list(label = '非AV瘘亚组(其他手术)',
       dat = d[d$surgery_type != 'av_fistula' & !is.na(d$asa3), ],
       covs = c('vol_le15','vol_ge25','age','asa3','surg_fracture','surg_other','shoulder','anes_general','inpatient'),
       note = '剔除anes_mac(MAC组18人0事件, 分离致系数无界); 保留fracture vs other_ue')
)
s <- fromJSON(file.path(OUT, 'sensitivity.json'), simplifyVector = FALSE)
for (fc in fits) {
  dd <- fc$dat
  fit <- coxph(mkform('time_primary', 'event_primary', fc$covs), data = dd)
  t <- cox_tab(fit)
  idx <- which(vapply(s$rows, function(r) r$label == fc$label, logical(1)))
  s$rows[[idx]] <- list(
    label = fc$label, n = nrow(dd), events = sum(dd$event_primary),
    epv = round(sum(dd$event_primary) / length(fc$covs), 1),
    converged = all(is.finite(coef(fit))),
    note = fc$note,
    events_by_group = grp_e(dd, dd$event_primary),
    le15 = get2(t, 'vol_le15'), ge25 = get2(t, 'vol_ge25'))
  cat(sprintf('%s n=%d e=%d epv=%.1f | le15 %.2f(%.2f-%.2f)p=%.4f | ge25 %.2f(%.2f-%.2f)p=%.4f\n',
              fc$label, nrow(dd), sum(dd$event_primary),
              sum(dd$event_primary) / length(fc$covs),
              t[t$term=='vol_le15','hr'], t[t$term=='vol_le15','lo'],
              t[t$term=='vol_le15','hi'], t[t$term=='vol_le15','p'],
              t[t$term=='vol_ge25','hr'], t[t$term=='vol_ge25','lo'],
              t[t$term=='vol_ge25','hi'], t[t$term=='vol_ge25','p']))
}
write_json(s, file.path(OUT, 'sensitivity.json'),
           auto_unbox = TRUE, pretty = TRUE, digits = 5)
cat('sensitivity.json patched\n')
