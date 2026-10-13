
suppressMessages({library(survival); library(jsonlite)})
d <- read.csv("Sys.getenv('OUT_DIR', '.')/analysis_571.csv", stringsAsFactors=FALSE)
d[["asa3"]] <- as.numeric(d[["asa3"]])
cc <- d[!is.na(d[["asa3"]]), ]
cc[["av"]] <- as.integer(cc[["surgery_type"]] == "av_fistula")
covs <- c("vol_le15","vol_ge25","age","asa3","surg_fracture","surg_other","shoulder","anes_general","anes_mac","inpatient")
f0 <- as.formula(paste("Surv(time_primary,event_primary) ~", paste(covs, collapse="+")))
f1 <- as.formula(paste("Surv(time_primary,event_primary) ~", paste(c(covs,"vol_le15:av","vol_ge25:av"), collapse="+")))
m0 <- coxph(f0, data=cc)
m1 <- coxph(f1, data=cc)
lrt <- anova(m0, m1)
chi2v <- as.numeric(lrt[[2]][2])
p2dfv <- as.numeric(lrt[[4]][2])
sm <- summary(m1)
ci <- sm[["conf.int"]]
co <- sm[["coefficients"]]
gr <- function(term) {
  i <- which(rownames(ci)==term)
  list(HR=as.numeric(round(ci[i,"exp(coef)"],3)), lo=as.numeric(round(ci[i,"lower .95"],3)),
       hi=as.numeric(round(ci[i,"upper .95"],3)), p=as.numeric(round(co[i,"Pr(>|z|)"],4)))
}
out <- list(n=nrow(d), events=sum(d[["event_primary"]]), model_n=nrow(cc), model_events=sum(cc[["event_primary"]]),
  interaction_terms=list(le15_av=gr("vol_le15:av"), ge25_av=gr("vol_ge25:av")),
  LRT_2df_chi2=round(chi2v,3), LRT_2df_p=round(p2dfv,4))
write_json(out, "Sys.getenv('OUT_DIR', '.')/interaction.json", auto_unbox=TRUE, pretty=TRUE, digits=5)
cat("LRT 2df p =", round(p2dfv,4), "
")
cat("le15:av p =", out, "
")
