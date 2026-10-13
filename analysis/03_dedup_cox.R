# ============================================================
# rerun5-03b: 同分钟去重后的主模型（稿中去重后 <=15mL HR 1.87）
# 输入：$OUT_DIR/rerun5/dedup_groups.csv（新分组）+ $OUT_DIR/analysis_571.csv
#       （协变量与结局；暴露改为去重后分组）
# 输出：$OUT_DIR/rerun5/dedup_cox.json
# ============================================================
suppressMessages(library(survival))
suppressMessages(library(jsonlite))
options(warn = 1)

OUT_DIR <- Sys.getenv("OUT_DIR", ".")
OUT <- file.path(OUT_DIR, 'rerun5')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(file.path(OUT_DIR, 'analysis_571.csv'),
              stringsAsFactors = FALSE)
d$asa3 <- as.numeric(d$asa3)
g <- read.csv(file.path(OUT, 'dedup_groups.csv'), stringsAsFactors = FALSE)

COVS10 <- c('vol_le15','vol_ge25','age','asa3','surg_fracture','surg_other',
            'shoulder','anes_general','anes_mac','inpatient')

# 用去重后分组替换暴露
d <- merge(d, g[, c('mrn','grp_new')], by = 'mrn', all.x = TRUE)
stopifnot(!any(is.na(d$grp_new)))
d$vol_le15 <- as.integer(d$grp_new == 'le15')
d$vol_ge25 <- as.integer(d$grp_new == 'ge25')

cc <- d[!is.na(d$asa3), ]
n <- nrow(cc); e <- sum(cc$event_primary)
cat(sprintf('deduped: n=%d events=%d\n', n, e))
cat(sprintf('groups: le15=%d ref=%d ge25=%d\n',
            sum(cc$vol_le15 == 1),
            sum(cc$vol_le15 == 0 & cc$vol_ge25 == 0),
            sum(cc$vol_ge25 == 1)))

fit <- coxph(as.formula(paste0('Surv(time_primary,event_primary) ~ ',
                              paste(COVS10, collapse = '+'))),
             data = cc)
s <- summary(fit)
tab <- data.frame(term = rownames(s$conf.int),
                  hr = round(s$conf.int[, 'exp(coef)'], 2),
                  lo = round(s$conf.int[, 'lower .95'], 2),
                  hi = round(s$conf.int[, 'upper .95'], 2),
                  p  = round(s$coefficients[, 'Pr(>|z|)'], 4),
                  stringsAsFactors = FALSE)
print(tab[tab$term %in% c('vol_le15','vol_ge25'), ])
cat('converged:', all(is.finite(coef(fit))), '\n')

get2 <- function(term) {
  r <- tab[tab$term == term, c('hr','lo','hi','p')]
  list(hr = r$hr, lo = r$lo, hi = r$hi, p = r$p)
}
res <- list(
  method = paste0('same 10-covariate cause-specific Cox as main analysis, ',
                  'exposure replaced by post-dedup volume groups; ',
                  'complete cases on asa3'),
  n = n, events = e,
  groups_n = list(le15 = sum(cc$vol_le15 == 1),
                  ref  = sum(cc$vol_le15 == 0 & cc$vol_ge25 == 0),
                  ge25 = sum(cc$vol_ge25 == 1)),
  groups_events = list(
    le15 = sum(cc$event_primary[cc$vol_le15 == 1]),
    ref  = sum(cc$event_primary[cc$vol_le15 == 0 & cc$vol_ge25 == 0]),
    ge25 = sum(cc$event_primary[cc$vol_ge25 == 1])),
  vol_le15 = get2('vol_le15'), vol_ge25 = get2('vol_ge25'),
  converged = all(is.finite(coef(fit))),
  script = file.path(OUT, '03_dedup_cox.R'),
  manuscript_target = list(vol_le15_hr = 1.87)
)
write_json(res, file.path(OUT, 'dedup_cox.json'),
           auto_unbox = TRUE, pretty = TRUE, digits = 5)
cat('wrote dedup_cox.json\n')
