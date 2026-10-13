# ============================================================
# rerun5-01: MICE 扩展模型（Table 3 "Extended model" 行，稿中 HR 1.97）
#
# 稿中定义（经手稿正文核对）：
#   扩展模型 = 主模型 10 协变量 + sex + bmi + 静脉地塞米松(dex_preend)，
#   共 13 协变量；缺失 BMI 用 chained equations 多重插补（m=20）。
#   dex_preend 缺失（无给药记录）按"未给药"=0 处理（与稿中
#   "249 patients (43.6%)" 的描述口径一致：249/571）。
#   asa3 缺失 10 例同样参与插补。全队列 n=571，73 事件。
#
# 方法：mice 包在运行环境中不可用（编译失败），改用 base R 手工实现 FCS：
#   asa3（二分类）→ 贝叶斯 logistic 插补；
#   bmi（连续）   → 贝叶斯线性回归插补；
#   插补模型均纳入其余全部协变量 + 结局（事件指示 + Nelson-Aalen
#   累积风险）；每轮插补 10 次 chained 迭代后保留；m=20，种子固定；
#   Rubin 法则合并 log(HR)。
#
# 输入：$OUT_DIR/analysis_571.csv（571 版主队列）
# 输出：$OUT_DIR/rerun5/mice_extended.json
# ============================================================
suppressMessages(library(survival))
suppressMessages(library(jsonlite))
options(warn = 1)

SEED <- 20261011
M <- 20
NITER <- 10
OUT_DIR <- Sys.getenv("OUT_DIR", ".")
OUT <- file.path(OUT_DIR, 'rerun5')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(file.path(OUT_DIR, 'analysis_571.csv'),
              stringsAsFactors = FALSE)
d$asa3 <- as.numeric(d$asa3)
d$bmi <- as.numeric(d$bmi)
d$sexF <- as.integer(d$sex == 'Female')
# dex 缺失按未给药=0（稿中 249/571 = 43.6% 口径）
d$dex <- ifelse(is.na(d$dex_preend) | d$dex_preend == '', 0,
                as.numeric(d$dex_preend))
stopifnot(sum(is.na(d$dex)) == 0)

COVS13 <- c('vol_le15','vol_ge25','age','asa3','surg_fracture','surg_other',
            'shoulder','anes_general','anes_mac','inpatient',
            'sexF','bmi','dex')

miss_tab <- data.frame(
  variable = COVS13,
  n_missing = vapply(COVS13, function(v) sum(is.na(d[[v]])), integer(1)),
  stringsAsFactors = FALSE)
miss_tab$pct_missing <- round(100 * miss_tab$n_missing / nrow(d), 2)
print(miss_tab)
cat(sprintf('dex: 1=%d, 0(含缺失按0)=%d\n',
            sum(d$dex == 1), sum(d$dex == 0)))
stopifnot(sum(is.na(d$time_primary)) == 0, sum(is.na(d$event_primary)) == 0)

# ---- 结局信息（Nelson-Aalen 累积风险） ----
sf <- survfit(Surv(time_primary, event_primary) ~ 1, data = d)
d$H <- approx(sf$time, sf$cumhaz, xout = d$time_primary, rule = 2)$y

pred_all <- c('vol_le15','vol_ge25','age','surg_fracture','surg_other',
              'shoulder','anes_general','anes_mac','inpatient',
              'sexF','dex','event_primary','H')

# ---- 贝叶斯插补单步（miss：原始缺失位置的逻辑向量，每轮只重插补这些位置） ----
imp_binary <- function(y, miss, X) {
  cc <- !miss  # 用观测行拟合（首轮初始化后 miss 行已有初值，不参与拟合）
  yfit <- y
  yfit[miss] <- NA
  fit <- glm(yfit ~ ., data = X, family = binomial)
  bhat <- coef(fit); V <- vcov(fit)
  bstar <- bhat + as.numeric(t(chol(V)) %*% rnorm(length(bhat)))
  newX <- model.matrix(delete.response(terms(fit)), data = X[miss, ])
  pr <- as.numeric(plogis(newX %*% bstar))
  y[miss] <- rbinom(sum(miss), 1, pr)
  y
}
imp_cont <- function(y, miss, X) {
  yfit <- y
  yfit[miss] <- NA
  fit <- lm(yfit ~ ., data = X)
  bhat <- coef(fit); V <- vcov(fit); s2 <- summary(fit)$sigma^2
  cc <- !miss
  ncc <- sum(cc); pp <- length(bhat)
  # sigma2* ~ scaled-inverse-chisq
  s2star <- s2 * (ncc - pp) / rchisq(1, ncc - pp)
  bstar <- bhat + as.numeric(t(chol(V * s2star / s2)) %*% rnorm(pp))
  newX <- model.matrix(delete.response(terms(fit)), data = X[miss, ])
  y[miss] <- as.numeric(newX %*% bstar) + rnorm(sum(miss), 0, sqrt(s2star))
  y
}

mk_imp_X <- function(dd, exclude) {
  # 插补设计矩阵：其余协变量 + 结局
  cols <- setdiff(pred_all, exclude)
  as.data.frame(dd[, cols])
}

set.seed(SEED)
miss_asa3 <- is.na(d$asa3)
miss_bmi <- is.na(d$bmi)
imp_list <- vector('list', M)
for (m in seq_len(M)) {
  dd <- d
  # 简单初始化
  dd$asa3[miss_asa3] <- rbinom(sum(miss_asa3), 1, mean(d$asa3, na.rm = TRUE))
  dd$bmi[miss_bmi] <- mean(d$bmi, na.rm = TRUE)
  # chained 迭代：每轮只重插补原始缺失位置
  for (it in seq_len(NITER)) {
    dd$asa3 <- imp_binary(dd$asa3, miss_asa3, mk_imp_X(dd, 'asa3'))
    dd$bmi  <- imp_cont(dd$bmi, miss_bmi, mk_imp_X(dd, 'bmi'))
  }
  imp_list[[m]] <- dd
}
cat(sprintf('FCS done: m=%d, %d iterations each\n', M, NITER))

# ---- 每个插补集跑扩展模型 ----
mkform <- function()
  as.formula(paste0('Surv(time_primary,event_primary) ~ ',
                    paste(COVS13, collapse = '+')))
ests <- t(vapply(imp_list, function(dd) {
  s <- summary(coxph(mkform(), data = dd))
  setNames(s$coefficients[, 'coef'], rownames(s$coefficients))
}, numeric(length(COVS13))))
colnames(ests) <- COVS13
ses <- t(vapply(imp_list, function(dd) {
  s <- summary(coxph(mkform(), data = dd))
  setNames(s$coefficients[, 'se(coef)'], rownames(s$coefficients))
}, numeric(length(COVS13))))
colnames(ses) <- COVS13

rubin <- function(b, se) {
  qbar <- mean(b); W <- mean(se^2); B <- var(b)
  Tvar <- W + (1 + 1/M) * B
  df <- (M - 1) * (1 + W / ((1 + 1/M) * B))^2
  ci <- qbar + c(-1, 1) * qt(0.975, df) * sqrt(Tvar)
  p <- 2 * pt(abs(qbar) / sqrt(Tvar), df, lower.tail = FALSE)
  list(hr = unname(exp(qbar)), lo = unname(exp(ci[1])), hi = unname(exp(ci[2])),
       p = unname(p), df = unname(df),
       fmi = unname((B + B/M) / Tvar))
}
pooled <- lapply(COVS13, function(v) rubin(ests[, v], ses[, v]))
names(pooled) <- COVS13

for (v in c('vol_le15','vol_ge25','bmi','dex','sexF')) {
  r <- pooled[[v]]
  cat(sprintf('%-10s HR=%.2f (%.2f-%.2f) p=%.4f fmi=%.3f\n',
              v, r$hr, r$lo, r$hi, r$p, r$fmi))
}

res <- list(
  method = paste0('manual FCS multiple imputation (mice package unavailable: ',
                  'compilation failed on ARM); asa3 (binary, n=10) via Bayesian ',
                  'logistic, bmi (continuous, n=209) via Bayesian linear regression; ',
                  'imputation models include all other covariates + event indicator + ',
                  'Nelson-Aalen cumulative hazard; dex_preend missing (n=264) coded 0 ',
                  '(no administration record); 10 chained iterations per imputation; ',
                  'm=20, seed=20261011; cause-specific Cox with 13 covariates per ',
                  'imputed dataset; Rubin pooling'),
  seed = SEED, m = M, n_iter = NITER,
  n_total = nrow(d), events_total = sum(d$event_primary),
  covariates = COVS13,
  missingness = miss_tab,
  pooled = lapply(pooled, function(r) list(
    hr = round(r$hr, 2), lo = round(r$lo, 2), hi = round(r$hi, 2),
    p = round(r$p, 4), fmi = round(r$fmi, 3))),
  script = file.path(OUT, '01_mice_extended.R'),
  manuscript_target = list(label = 'Table 3 Extended model (+sex, BMI, IV dex), 73/571',
                           le15_hr = 2.06, le15_ci = c(0.99, 4.25))
)
write_json(res, file.path(OUT, 'mice_extended.json'),
           auto_unbox = TRUE, pretty = TRUE, digits = 5)
cat('wrote mice_extended.json\n')
