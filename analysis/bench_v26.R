suppressMessages({library(survival); library(coxphf)})
d <- read.csv('Sys.getenv('OUT_DIR', '.')/analysis_571.csv', stringsAsFactors=FALSE)
d$asa3 <- as.numeric(d$asa3)
cc <- d[!is.na(d$asa3), ]
cat('n cc:', nrow(cc), ' events:', sum(cc$event_primary), '\n')
covs <- c('vol_le15','vol_ge25','age','asa3','surg_fracture','surg_other','shoulder','anes_general','anes_mac','inpatient')
f <- as.formula(paste('Surv(time_primary, event_primary) ~', paste(covs, collapse='+')))
fit <- coxph(f, data=cc)
s <- summary(fit)
i <- which(rownames(s$conf.int)=='vol_le15')
cat(sprintf('COX le15: HR=%.3f (%.2f-%.2f) p=%.4f\n', s$conf.int[i,'exp(coef)'], s$conf.int[i,'lower .95'], s$conf.int[i,'upper .95'], s$coefficients[i,'Pr(>|z|)']))
ff <- coxphf(f, data=cc, maxit=500)
ci <- exp(confint(ff))
cat(sprintf('FIRTH le15: HR=%.3f (%.2f-%.2f) p=%.4f\n', exp(ff$coefficients['vol_le15']), ci['vol_le15',1], ci['vol_le15',2], ff$prob['vol_le15']))
d6 <- cc; d6$time_s <- pmin(d6$time_primary, 6); d6$event_s <- ifelse(d6$time_primary<=6, d6$event_primary, 0)
cat('0-6h events:', sum(d6$event_s), '\n')
cat('0-6h events by group le15/ref/ge25:',
    sum(d6$event_s[d6$vol_le15==1]), '/',
    sum(d6$event_s[d6$vol_le15==0 & d6$vol_ge25==0]), '/',
    sum(d6$event_s[d6$vol_ge25==1]), '\n')
f6 <- as.formula(paste('Surv(time_s, event_s) ~', paste(covs, collapse='+')))
ff6 <- coxphf(f6, data=d6, maxit=500)
ci6 <- exp(confint(ff6))
cat(sprintf('FIRTH 0-6h le15: HR=%.3f (%.2f-%.2f) p=%.4f\n', exp(ff6$coefficients['vol_le15']), ci6['vol_le15',1], ci6['vol_le15',2], ff6$prob['vol_le15']))
dl <- cc[cc$time_primary>6, ]; dl$time_l <- dl$time_primary-6
cat('>6h landmark n:', nrow(dl), ' events:', sum(dl$event_primary), '\n')
cat('>6h events by group le15/ref/ge25:',
    sum(dl$event_primary[dl$vol_le15==1]), '/',
    sum(dl$event_primary[dl$vol_le15==0 & dl$vol_ge25==0]), '/',
    sum(dl$event_primary[dl$vol_ge25==1]), '\n')
fl <- as.formula(paste('Surv(time_l, event_primary) ~', paste(covs, collapse='+')))
ffl <- coxphf(fl, data=dl, maxit=500)
cil <- exp(confint(ffl))
cat(sprintf('FIRTH >6h le15: HR=%.3f (%.2f-%.2f) p=%.4f\n', exp(ffl$coefficients['vol_le15']), cil['vol_le15',1], cil['vol_le15',2], ffl$prob['vol_le15']))
