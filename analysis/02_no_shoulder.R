# ============================================================
# rerun5-02: 排除肩部手术（稿中 n=543，<=15mL HR 2.25）
# 方法：在 shoulder==0 的人群中跑主模型，去掉 shoulder 协变量
#       （9 协变量 cause-specific Cox）。
# 输入：$OUT_DIR/analysis_571.csv（571 版主队列）
# 输出：$OUT_DIR/rerun5/no_shoulder.json
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

COVS9 <- c('vol_le15','vol_ge25','age','asa3','surg_fracture','surg_other',
           'anes_general','anes_mac','inpatient')  # 去掉 shoulder

sub <- d[d$shoulder == 0, ]
cc <- sub[!is.na(sub$asa3), ]
n <- nrow(cc); e <- sum(cc$event_primary)
cat(sprintf('no-shoulder: n=%d events=%d (excluded shoulder n=%d)\n',
            n, e, sum(d$shoulder == 1, na.rm = TRUE)))

grp_n <- c(le15 = sum(cc$vol_le15 == 1),
           ref  = sum(cc$vol_le15 == 0 & cc$vol_ge25 == 0),
           ge25 = sum(cc$vol_ge25 == 1))
grp_e <- c(le15 = sum(cc$event_primary[cc$vol_le15 == 1]),
           ref  = sum(cc$event_primary[cc$vol_le15 == 0 & cc$vol_ge25 == 0]),
           ge25 = sum(cc$event_primary[cc$vol_ge25 == 1]))

fit <- coxph(as.formula(paste0('Surv(time_primary,event_primary) ~ ',
                              paste(COVS9, collapse = '+'))),
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
  method = paste0('exclude shoulder==1; cause-specific Cox with 9 covariates ',
                  '(main set minus shoulder); complete cases on asa3'),
  n = n, events = e,
  groups_n = as.list(grp_n), groups_events = as.list(grp_e),
  vol_le15 = get2('vol_le15'), vol_ge25 = get2('vol_ge25'),
  converged = all(is.finite(coef(fit))),
  script = file.path(OUT, '02_no_shoulder.R'),
  manuscript_target = list(n = 541, vol_le15_hr = 2.12)
)
write_json(res, file.path(OUT, 'no_shoulder.json'),
           auto_unbox = TRUE, pretty = TRUE, digits = 5)
cat('wrote no_shoulder.json\n')
