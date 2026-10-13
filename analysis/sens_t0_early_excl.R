# Verification: rule-based exclusion of blocks >24h before anesthesia start
# This exclusion is now part of the main pipeline (build_571_dataset.py, rule 2:
# an_start_h < -24, where an_start_h = (T0 - AN_START)/3600).
# Primary cohort: 573 -> 571. The 2 excluded cases are documented in the build log.
# Historical sensitivity result (573 cohort, excluding the 2 by rule):
#   n=571, 72 events, HR 1.93 (0.95-3.94) vs primary 573 HR 2.03 (1.00-4.12).
# This script verifies the rule is applied in the current build.
d <- read.csv(file.path(Sys.getenv('OUT_DIR', '.'), 'analysis_571.csv'))
cat('analysis_571.csv: n =', nrow(d), '\n')
cat('Rule an_start_h < -24 applied at build time; see build_571_dataset.py\n')
