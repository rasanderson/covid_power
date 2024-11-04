# Power analysis to estimate number of samples needed for Covid in non-detects
set.seed(123)
library(ggplot2)
library(vcd)

rm(list = ls())

# Data from Ben 30 October 2024
obs_data <- data.frame(
  age_group = rep(c("0-14", "15-24", "25-64", "65+"), 2),
  pathogen  = rep(c("no_path_det", "path_det"), each = 4),
  counts    = c(21, 31, 511, 321, 51, 51, 392, 210)
)

# Check for differences in age groups
# Create a contingency table
contingency_table <- xtabs(counts ~ age_group + pathogen, data = obs_data)
# Perform the chi-squared test
chi_squared_test <- chisq.test(contingency_table)
print(chi_squared_test)
# Perform pairwise comparisons
pairwise_test <- pairwise.prop.test(contingency_table, p.adjust.method = "bonferroni")
print(pairwise_test)

# Bar plot of counts by age group and pathogen status
ggplot(obs_data, aes(x = age_group, y = counts, fill = pathogen)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Counts by Age Group and Pathogen Status", x = "Age Group", y = "Counts") +
  theme_classic()


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
max_samp <- 3000
stp_samp <- 50
samp_range <- seq(min_samp, max_samp, by = stp_samp)
# Draws per sample 
drw_samp <- 10

# Ignore age structure initially
# Repeat the process for each sample size based on drw_samp
p_val_tally <- data.frame(drw_samp = numeric(0), trial = numeric(0), p_val = numeric(0))
for(trial in samp_range){
  for(i in 1 : drw_samp){
    det_thresh <- runif(1, min = min_pop/100, max = max_pop/100)
    det_with_covid <- ifelse(runif(trial, 0, 1) > det_thresh, 0, 1)
    nondet_with_covid <- ifelse(runif(trial, 0, 1) > (det_thresh * non_det_xs ), 0, 1)
    # Assemble into data.frame
    dat_for_check <- data.frame(
      freq = c(det_with_covid, nondet_with_covid),
      source = c(rep("det", trial), rep("nondet", trial))
    )
    mod_res <- glm(freq ~ source, data = dat_for_check, family = "binomial")
    p_val <- summary(mod_res)$coefficients[2, 4]
    p_val_tally <- rbind(p_val_tally, data.frame(drw_samp = i, trial = trial, p_val = p_val))
  }
}

ggplot(p_val_tally, aes(x = trial, y = p_val)) +
  geom_smooth() +
  ylim(0, NA) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  labs(title = "Power Analysis for Covid in Non-Detects", subtitle = "Asssume 10% excess; no age structure", x = "Number of Samples", y = "P-value") +
  theme_classic()


# Again ignore age structure, but give range of % increase in Covid in non-detects
# is actually causing their IID
non_det_xs_pct <- c(5, 10, 15, 20, 25) # % increase in Covid
non_det_xs <- non_det_xs_pct/100 + 1
# Sample range to check
min_samp <- 250
max_samp <- 3000
stp_samp <- 50
samp_range <- seq(min_samp, max_samp, by = stp_samp)
# Draws per sample
drw_samp <- 10

# Ignore age structure initially
# Repeat the process for each sample size based on drw_samp
p_val_tally <- data.frame(non_det_xs = numeric(0), drw_samp = numeric(0), trial = numeric(0), p_val = numeric(0))
for(xs in non_det_xs){
  for(trial in samp_range){
    for(i in 1 : drw_samp){
      det_thresh <- runif(1, min = min_pop/100, max = max_pop/100)
      det_with_covid <- ifelse(runif(trial, 0, 1) > det_thresh, 0, 1)
      nondet_with_covid <- ifelse(runif(trial, 0, 1) > (det_thresh * xs ), 0, 1)
      # Assemble into data.frame
      dat_for_check <- data.frame(
        freq = c(det_with_covid, nondet_with_covid),
        source = c(rep("det", trial), rep("nondet", trial))
      )
      mod_res <- glm(freq ~ source, data = dat_for_check, family = "binomial")
      p_val <- summary(mod_res)$coefficients[2, 4]
      p_val_tally <- rbind(p_val_tally, data.frame(non_det_xs = xs, drw_samp = i, trial = trial, p_val = p_val))
    }
  }
}

# Smoothed ggplot with different lines for each non_det_xs
p_val_tally$non_det_xs <- as.factor(p_val_tally$non_det_xs)
ggplot(p_val_tally, aes(x = trial, y = p_val, color = non_det_xs)) +
  geom_smooth() +
  ylim(0, NA) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  labs(title = "Power Analysis for Covid in Non-Detects", subtitle = "Asssume 5-25% excess; no age structure", x = "Number of Samples", y = "P-value") +
  theme_classic()


# Now include age structure
# Data from Ben 30 October 2024
# Age groupings
# 0-14, 15-24, 25-64, 65+
age_props <- rowSums(contingency_table) / sum(rowSums(contingency_table))
non_det_xs_pct <- c(5, 10, 15, 20, 25) # % increase in Covid
non_det_xs <- non_det_xs_pct/100 + 1
# Sample range to check
min_samp <- 250
max_samp <- 3000
stp_samp <- 50
samp_range <- seq(min_samp, max_samp, by = stp_samp)
# Draws per sample
drw_samp <- 10

# Repeat the process for each sample size based on drw_samp
p_val_tally <- data.frame(non_det_xs = numeric(0), drw_samp = numeric(0), trial = numeric(0), p_val = numeric(0))
for(xs in non_det_xs){
  for(trial in samp_range){
    for(i in 1 : drw_samp){
      det_thresh <- runif(1, min = min_pop/100, max = max_pop/100)
      det_with_covid <- ifelse(runif(trial, 0, 1) > det_thresh, 0, 1)
      nondet_with_covid <- ifelse(runif(trial, 0, 1) > (det_thresh * xs ), 0, 1)
      age_class <- ifelse(runif(trial, 0, 1) < age_props[1], "0-14",
                          ifelse(runif(trial, 0, 1) < sum(age_props[1:2]), "15-24",
                                 ifelse(runif(trial, 0, 1) < sum(age_props[1:3]), "25-64", "65+")))
      # Assemble into data.frame
      dat_for_check <- data.frame(
        freq = c(det_with_covid, nondet_with_covid),
        source = c(rep("det", trial), rep("nondet", trial)),
        age_class = age_class
      )
      mod_res <- glm(freq ~ source + age_class, data = dat_for_check, family = "binomial")
      p_val <- summary(mod_res)$coefficients[2, 4]
      p_val_tally <- rbind(p_val_tally, data.frame(non_det_xs = xs, drw_samp = i, trial = trial, p_val = p_val))
    }
  }
}

# Smoothed ggplot with different lines for each non_det_xs
p_val_tally$non_det_xs <- as.factor(p_val_tally$non_det_xs)
ggplot(p_val_tally, aes(x = trial, y = p_val, color = non_det_xs)) +
  geom_smooth() +
  ylim(0, NA) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  labs(title = "Power Analysis for Covid in Non-Detects", subtitle = "Asssume 5-25% excess; Age structure included", x = "Number of Samples", y = "P-value") +
  theme_classic()


