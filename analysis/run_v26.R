# ============================================================
# MOVER 主队列 (n=571; 完整病例 n=561) 全部模型
# Table 2 (models.json) / Table 3 (sensitivity.json) / A9 (firth_ow.json)
# ============================================================
suppressMessages({ library(survival); library(coxphf); library(cmprsk) })
if (!requireNamespace('jsonlite', quietly = TRUE))
  install.packages('jsonlite', repos = 'https://cloud.r-project.org', quiet = TRUE)
suppressMessages(library(jsonlite))
options(warn = 1)

OUT <- Sys.getenv('OUT_DIR', '.')
d <- read.csv(file.path(OUT, 'analysis_571.csv'), stringsAsFactors = FALSE)
d$asa3 <- as.numeric(d$asa3)
d$bmi <- as.numeric(d$bmi)
d$dex_preend <- as.numeric(d$dex_preend)
d$sexF <- as.integer(d$sex == 'Female')

COVS10 <- c('vol_le15','vol_ge25','age','asa3','surg_fracture','surg_other',
            'shoulder','anes_general','anes_mac','inpatient')
COVS8  <- setdiff(COVS10, c('vol_le15','vol_ge25'))   # PS 协变量
PSCOV  <- COVS8

cc <- d[!is.na(d$asa3), ]
N_CC <- nrow(cc); E_CC <- sum(cc$event_primary)
cat(sprintf('cohort: n=%d events=%d | complete: n=%d events=%d\n',
            nrow(d), sum(d$event_primary), N_CC, E_CC))

# T0=术毕: 非正时间修正
npos <- sum(d$t0end_ok == 1 & !is.na(d$time_t0end) & d$time_t0end <= 0, na.rm = TRUE)
if (npos > 0) {
  mp <- min(d$time_t0end[d$t0end_ok == 1 & d$time_t0end > 0], na.rm = TRUE) / 2
  d$time_t0end[d$t0end_ok == 1 & !is.na(d$time_t0end) & d$time_t0end <= 0] <- mp
}
cat('t0end non-positive fixed:', npos, '\n')

mkform <- function(tvar, evar, covs)
  as.formula(paste0('Surv(', tvar, ',', evar, ') ~ ', paste(covs, collapse = '+')))
ok_fit <- function(fit) !is.null(fit) && all(is.finite(coef(fit)))

cox_tab <- function(fit) {
  s <- summary(fit)
  data.frame(term = rownames(s$conf.int),
             hr = round(s$conf.int[, 'exp(coef)'], 2),
             lo = round(s$conf.int[, 'lower .95'], 2),
             hi = round(s$conf.int[, 'upper .95'], 2),
             p  = round(s$coefficients[, 'Pr(>|z|)'], 4),
             stringsAsFactors = FALSE)
}
firth_tab <- function(ff) {
  data.frame(term = names(ff$coefficients),
             hr = round(as.numeric(exp(ff$coefficients)), 2),
             lo = round(as.numeric(ff$ci.lower), 2),
             hi = round(as.numeric(ff$ci.upper), 2),
             p  = round(as.numeric(ff$prob), 4),
             stringsAsFactors = FALSE)
}
get2 <- function(tab, term) {
  r <- tab[tab$term == term, c('hr','lo','hi','p')]
  list(hr = r$hr, lo = r$lo, hi = r$hi, p = r$p)
}
grp_n <- function(dd, a = 'vol_le15', b = 'vol_ge25')
  list(le15 = sum(dd[[a]] == 1), ref = sum(dd[[a]] == 0 & dd[[b]] == 0), ge25 = sum(dd[[b]] == 1))
grp_e <- function(dd, ev, a = 'vol_le15', b = 'vol_ge25')
  list(le15 = sum(ev[dd[[a]] == 1]), ref = sum(ev[dd[[a]] == 0 & dd[[b]] == 0]),
       ge25 = sum(ev[dd[[b]] == 1]))

## ---------------- A3 Table 2 ----------------
cohort <- list(n_total = nrow(d), events_total = sum(d$event_primary),
               n_complete = N_CC, events_complete = E_CC,
               groups_n = grp_n(d), groups_events = grp_e(d, d$event_primary),
               epv_main = round(E_CC / length(COVS10), 1))

# 粗 HR: 单因素 Cox, available-case
crude <- do.call(rbind, lapply(COVS10, function(v) {
  dd <- d[!is.na(d[[v]]), ]
  fit <- coxph(mkform('time_primary', 'event_primary', v), data = dd)
  r <- cox_tab(fit); r$n <- nrow(dd); r$events <- sum(dd$event_primary); r
}))

# 校正 HR: 10 协变量主模型
fit_main <- coxph(mkform('time_primary', 'event_primary', COVS10), data = cc)
adjusted <- cox_tab(fit_main); adjusted$n <- N_CC; adjusted$events <- E_CC

# 逐个剔除
loo <- do.call(rbind, lapply(COVS10, function(v) {
  fit <- coxph(mkform('time_primary', 'event_primary', setdiff(COVS10, v)), data = cc)
  data.frame(dropped = v, le15_hr = round(unname(exp(coef(fit)['vol_le15'])), 3),
             stringsAsFactors = FALSE)
}))
loo_range <- list(min = min(loo$le15_hr, na.rm = TRUE), max = max(loo$le15_hr, na.rm = TRUE),
                  at_min = loo$dropped[which.min(loo$le15_hr)],
                  at_max = loo$dropped[which.max(loo$le15_hr)],
                  note = 'dropping vol_le15 itself gives no le15 estimate (excluded); range over the other 9 covariates')

# Schoenfeld
zph <- cox.zph(fit_main)
schoenfeld <- data.frame(term = rownames(zph$table),
                         p = round(zph$table[, 'p'], 4), stringsAsFactors = FALSE)
rownames(schoenfeld) <- NULL

# 时间分段 Firth (0-6h 保留所有行; >6h landmark)
d6 <- cc
d6$ts <- pmin(d6$time_primary, 6)
d6$es <- ifelse(d6$time_primary <= 6, d6$event_primary, 0)
ff6 <- coxphf(mkform('ts', 'es', COVS10), data = d6, maxit = 500)
t6 <- firth_tab(ff6)
dl <- cc[cc$time_primary > 6, ]; dl$tl <- dl$time_primary - 6
ffl <- coxphf(mkform('tl', 'event_primary', COVS10), data = dl, maxit = 500)
tl <- firth_tab(ffl)
timesplit <- list(
  h0_6    = list(n = nrow(d6), events = sum(d6$es),
                 events_by_group = grp_e(d6, d6$es),
                 le15 = get2(t6, 'vol_le15'), ge25 = get2(t6, 'vol_ge25')),
  h6plus  = list(n = nrow(dl), events = sum(dl$event_primary),
                 events_by_group = grp_e(dl, dl$event_primary),
                 le15 = get2(tl, 'vol_le15'), ge25 = get2(tl, 'vol_ge25'))
)

models <- list(cohort = cohort,
               crude = crude, adjusted = adjusted,
               leave_one_out = list(range = loo_range, details = loo),
               schoenfeld = schoenfeld,
               timesplit_firth = timesplit,
               expanded_mice = NULL,
               mice_available = requireNamespace('mice', quietly = TRUE))
write_json(models, file.path(OUT, 'models.json'),
           auto_unbox = TRUE, pretty = TRUE, digits = 5)
cat('models.json written\n')

## ---------------- A9: Firth 全表 + PS 重叠加权 ----------------
ffm <- coxphf(mkform('time_primary', 'event_primary', COVS10), data = cc, maxit = 500)
firth_full <- firth_tab(ffm); firth_full$n <- N_CC; firth_full$events <- E_CC

psd <- cc[cc$vol_ge25 == 0, ]; psd$trt <- psd$vol_le15
cat(sprintf('PS subset: n=%d events=%d (le15 n=%d)\n',
            nrow(psd), sum(psd$event_primary), sum(psd$trt == 1)))
psf <- glm(as.formula(paste('trt ~', paste(PSCOV, collapse = '+'))),
           data = psd, family = binomial)
ps <- as.numeric(predict(psf, type = 'response'))
w <- ifelse(psd$trt == 1, 1 - ps, ps)
wfit <- coxph(Surv(time_primary, event_primary) ~ vol_le15, data = psd,
              weights = w, robust = TRUE)
sw <- summary(wfit); b <- sw$coefficients['vol_le15', ]
ps_le15 <- list(hr = round(exp(b['coef']), 2),
                lo = round(exp(b['coef'] - 1.96 * b['robust se']), 2),
                hi = round(exp(b['coef'] + 1.96 * b['robust se']), 2),
                p  = round(b['Pr(>|z|)'], 4))
smd_w <- function(x, trt, w) {
  m1 <- weighted.mean(x[trt == 1], w[trt == 1]); m0 <- weighted.mean(x[trt == 0], w[trt == 0])
  v1 <- sum(w[trt == 1] * (x[trt == 1] - m1)^2) / sum(w[trt == 1])
  v0 <- sum(w[trt == 0] * (x[trt == 0] - m0)^2) / sum(w[trt == 0])
  if ((v1 + v0) == 0) return(0)
  abs(m1 - m0) / sqrt((v1 + v0) / 2)
}
smds <- vapply(PSCOV, function(v) smd_w(psd[[v]], psd$trt, w), numeric(1))
ps_ow <- list(report_n = N_CC, report_events = E_CC,
              subset_n = nrow(psd), subset_events = sum(psd$event_primary),
              le15 = ps_le15,
              max_abs_smd_weighted = round(max(smds), 3),
              smd_by_covariate = round(smds, 3),
              method = 'logistic PS (8 covariates) -> overlap weights w=1-PS(trt=1), w=PS(trt=0); weighted Cox (vol_le15 only), robust SE')
write_json(list(firth_full = firth_full, ps_overlap_weighting = ps_ow),
           file.path(OUT, 'firth_ow.json'), auto_unbox = TRUE, pretty = TRUE, digits = 5)
cat('firth_ow.json written\n')
cat(sprintf('PS-OW le15: HR=%.2f (%.2f-%.2f) p=%.4f max|SMD|=%.3f\n',
            ps_le15$hr, ps_le15$lo, ps_le15$hi, ps_le15$p, max(smds)))

## ---------------- A4 Table 3 ----------------
sens_row <- function(label, dat, tvar = 'time_primary', evar = 'event_primary',
                     covs = COVS10, le = 'vol_le15', ge = 'vol_ge25',
                     ga = 'vol_le15', gb = 'vol_ge25', note = '') {
  cols <- unique(c(tvar, evar, covs, le, ge, ga, gb))
  dd <- dat[complete.cases(dat[, cols]), ]
  n <- nrow(dd); e <- sum(dd[[evar]])
  r <- list(label = label, n = n, events = e, epv = round(e / length(covs), 1),
            converged = TRUE, note = note,
            events_by_group = grp_e(dd, dd[[evar]], ga, gb))
  if (e < 5) {
    r$converged <- FALSE
    r$note <- paste0(note, '; 事件数<5, 未建模')
    return(r)
  }
  fit <- tryCatch(suppressWarnings(coxph(mkform(tvar, evar, covs), data = dd)),
                  error = function(err) NULL)
  if (!ok_fit(fit)) {
    r$converged <- FALSE
    r$note <- paste0(note, '; 模型未收敛/系数非有限')
    return(r)
  }
  t <- cox_tab(fit)
  r$le15 <- get2(t, le); r$ge25 <- get2(t, ge)
  r
}

rows <- list()
rows[[length(rows) + 1]] <- sens_row('主分析(10协变量校正)', cc)
# Firth 行
rows[[length(rows) + 1]] <- list(label = 'Firth惩罚', n = N_CC, events = E_CC,
  epv = round(E_CC / length(COVS10), 1), converged = TRUE,
  note = 'coxphf, profile-likelihood CI',
  events_by_group = grp_e(cc, cc$event_primary),
  le15 = get2(firth_full, 'vol_le15'), ge25 = get2(firth_full, 'vol_ge25'))
# PS-OW 行 (统一 73/563 口径)
rows[[length(rows) + 1]] <- list(label = 'PS重叠加权', n = N_CC, events = E_CC,
  converged = TRUE,
  note = sprintf('overlap-weighted Cox; 实际估计用 ≤15mL vs >15-<25mL 子集 (n=%d, 事件=%d); 加权后最大|SMD|=%.3f',
                 ps_ow$subset_n, ps_ow$subset_events, ps_ow$max_abs_smd_weighted),
  le15 = ps_ow$le15, ge25 = list(hr = NA, lo = NA, hi = NA, p = NA))
# Fine-Gray
mm <- model.matrix(as.formula(paste('~', paste(COVS10, collapse = '+'))), data = cc)[, -1]
fg <- crr(ftime = cc$fg_time, fstatus = cc$fg_status, cov1 = mm)
sfg <- summary(fg)
fgt <- data.frame(term = rownames(sfg$coef),
                  hr = round(sfg$conf.int[, 'exp(coef)'], 2),
                  lo = round(sfg$conf.int[, '2.5%'], 2),
                  hi = round(sfg$conf.int[, '97.5%'], 2),
                  p  = round(sfg$coef[, 'p-value'], 4), stringsAsFactors = FALSE)
rows[[length(rows) + 1]] <- list(label = 'Fine-Gray(出院为竞争事件)', n = N_CC, events = E_CC,
  converged = isTRUE(fg$converged),
  note = '竞争事件=出院(491); 96h删失(8); 主要事件(74)',
  events_by_group = grp_e(cc, cc$event_primary),
  le15 = get2(fgt, 'vol_le15'), ge25 = get2(fgt, 'vol_ge25'))
# 最大单次剂量
COVS10B <- c('max_le15','max_ge25', COVS8)
rows[[length(rows) + 1]] <- sens_row('最大单次剂量定义', cc, covs = COVS10B,
  le = 'max_le15', ge = 'max_ge25', ga = 'max_le15', gb = 'max_ge25',
  note = '暴露=最大单次经神经周围给药量分组')
# 参照限精确20mL
cc20 <- cc[!(cc$cum60 > 15 & cc$cum60 < 25 & abs(cc$cum60 - 20) > 1e-9), ]
cat(sprintf('ref20: excluded %d (15<cum60<25 & !=20mL); n=%d events=%d; ref(=20mL) n=%d\n',
            N_CC - nrow(cc20), nrow(cc20), sum(cc20$event_primary),
            sum(abs(cc20$cum60 - 20) < 1e-9)))
rows[[length(rows) + 1]] <- sens_row('参照限精确20mL', cc20,
  note = '参照组=恰为20mL; 剔除15-25mL中非20mL者')
# 24h截断
rows[[length(rows) + 1]] <- sens_row('24h截断', cc, tvar = 'time_24', evar = 'event_24')
# T0=术毕
d_te <- d[d$t0end_ok == 1, ]
rows[[length(rows) + 1]] <- sens_row('T0=术毕', d_te, tvar = 'time_t0end', evar = 'event_t0end',
  note = sprintf('12例缺AN_STOP剔除; %d例非正时间取半最小正值', npos))
# 疼痛子集
rows[[length(rows) + 1]] <- sens_row('有疼痛评分子集', d[d$has_pain == 1, ])
# 住院亚组
rows[[length(rows) + 1]] <- sens_row('住院亚组', d[d$inpatient == 1, ],
  covs = setdiff(COVS10, 'inpatient'), note = 'inpatient恒定, 剔除该协变量')
# AV瘘亚组
rows[[length(rows) + 1]] <- sens_row('AV瘘亚组', d[d$surgery_type == 'av_fistula', ],
  covs = setdiff(COVS10, c('surg_fracture','surg_other')), note = '手术类型恒定, 剔除手术协变量')
# 非AV瘘亚组
rows[[length(rows) + 1]] <- sens_row('非AV瘘亚组(其他手术)', d[d$surgery_type != 'av_fistula', ],
  note = 'fracture vs other_ue')
# PRN限定: 不可行
rows[[length(rows) + 1]] <- list(label = 'PRN限定', feasible = FALSE,
  note = '不可行: ORDER_CLASS_NM取值为Inpatient/ePrescribe/Normal等, 无法区分PRN与常规医嘱(既往核查结论); 主要结局采用疼痛锚定法构建, 无PRN标识可用')
# 次要结局: 首次术后阿片
rows[[length(rows) + 1]] <- sens_row('次要结局(首次术后阿片)', cc,
  tvar = 'time_opioid', evar = 'event_opioid',
  note = '术后=AN_STOP之后首次Given阿片(缺AN_STOP者用T0); 删失=出院/96h')

sensitivity <- list(rows = rows, expanded_mice_pending = !models$mice_available)
write_json(sensitivity, file.path(OUT, 'sensitivity.json'),
           auto_unbox = TRUE, pretty = TRUE, digits = 5)
cat('sensitivity.json written\n')

## 控制台摘要
for (r in rows) {
  if (isTRUE(r$feasible == FALSE)) { cat(sprintf('%-28s 不可行: %s\n', r$label, r$note)); next }
  le <- r$le15; ge <- r$ge25
  cat(sprintf('%-28s n=%d e=%d | le15 %.2f(%.2f-%.2f)p=%.4f | ge25 %.2f(%.2f-%.2f)p=%.4f %s\n',
              r$label, r$n, r$events,
              le$hr, le$lo, le$hi, le$p, ge$hr, ge$lo, ge$hi, ge$p,
              ifelse(isTRUE(r$converged), '', '[未收敛]')))
}
cat('ALL DONE\n')
