# Power analysis to estimate number of samples needed for Covid in non-detects
set.seed(123)
rm(list = ls())

# Data from Ben 30 October 2024
obs_data <- data.frame(
  age_group = rep(c("0-14", "15-24", "25-64", "65+"), 2),
  pathogen  = rep(c("no_path_det", "path_det"), each = 4),
  counts    = c(21, 31, 511, 321, 51, 51, 392, 210)
)


# Data from https://www.gov.uk/government/statistics/national-flu-and-covid-19-surveillance-reports-2024-to-2025-season 
# Covid Pillar 1 data for 7-day moving average typically roughly 5 to 15% for
# all adults, ignoring age, sex, ethnicity etc.
# Population level covid
min_pop <- 5
max_pop <- 15

# Define % increase in Covid in non-detects, on the assumption that this Covid
# is actually causing their IID
non_det_xs_pct <- 10 # % increase in Covid
non_det_xs <- non_det_xs_pct/100 + 1
# Sample range to check
min_samp <- 250
max_samp <- 2000
stp_samp <- 50
samp_range <- seq(min_samp, max_samp, by = stp_samp)
# Draws per sample (not implemented yet)
drw_samp <- 100

for(trial in samp_range){
  det_thresh <- runif(1, min = min_pop/100, max = max_pop/100)
  det_with_covid <- ifelse(runif(trial, 0, 1) > det_thresh, 0, 1)
  nondet_with_covid <- ifelse(runif(trial, 0, 1) > (det_thresh * non_det_xs ), 0, 1)
  # Assemble into data.frame
  dat_for_check <- data.frame(
    freq = c(det_with_covid, nondet_with_covid),
    source = c(rep("det", trial), rep("nondet", trial))
  )
  # Structure of data not right or maybe need a chi-sq
  mod_res <- MASS::loglm(freq ~ source, data = dat_for_check)
  # Count whether mod_res is significant and tally for different trial sizes
}
