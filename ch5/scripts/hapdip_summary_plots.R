#!/usr/bin/env Rscript
library(tidyverse)

# --- Functions to get file paths ---

#return the data directory appended to path
get_data_dir <- function(path){
  return(paste0("../data/", path))
}

#return file name of the roc tsv
#given the prefix. For example, 0_0
#returns "../data/0_0.snp_roc.tsv.gz"
get_summary_path <- function(prefix){
  return(get_data_dir(paste0(prefix, ".summary.txt")))
}

# --- Import Dataframe ---
import_fnrfpr_tsvs <- function(fnr, fpr){
  prefixes <- unlist(map(fnr, paste0, '_', fpr))
  summary_files <- get_summary_path(prefixes)
  dfs <- map(summary_files, read_tsv,
             col_names = c("mode","type","operator","value"))
  names(dfs) <- prefixes
  # notempty <- map_lgl(dfs, ~dim(.)[1] != 0)
  # dfs <- dfs[notempty]
  fnr_var <- sym('%FNR')
  fdr_var <- sym('%FDR')
  bind_rows(dfs, .id = 'FNR_FPR') %>%
    pivot_wider(names_from = operator, values_from = value) %>%
    separate(FNR_FPR, into = c("FalseNegativeRate","FalsePositiveRate"),
             sep = "_", convert = TRUE) %>%
    rename(FNR = all_of(fnr_var)) %>%
    rename(FDR = all_of(fdr_var)) %>%
    mutate(FDR = FDR/100) %>%
    mutate(FNR = FNR/100)
}

# --- Manipulate Dataframe ---
only_snps <- function(df){
  df %>%
    filter(mode == "rtgeval-gt") %>%
    filter(type == "SNP") %>%
    select(-type, -mode)
}

get_f_statistic <- function(df){
  df %>%
    group_by(FalseNegativeRate, FalsePositiveRate) %>%
    # mutate(precision = 1 - FDR) %>%
    # mutate(recall = 1 - FNR) %>%
    mutate(recall = TPc / (TPc + FN)) %>%
    mutate(precision = TPc / (TPc + FP)) %>%
    mutate(f = 2 * precision * recall / (precision + recall))
}

calc_f_and_filter <- function(df){
  df %>% only_snps() %>% filter(FalsePositiveRate != 100) %>% get_f_statistic()
}

# --- Main ---
fnr <- seq(0,100,20)
df <- import_fnrfpr_tsvs(fnr,fnr) %>%
  calc_f_and_filter()

# --- Plotting ---

#Better calibrated data has more TPs
pdf("../figures/tp_heatmap.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
  geom_raster(aes(fill = TP)) +
  scale_fill_viridis_c("TP Calls") +
  scale_x_continuous("Variable Sites FNR") +
  scale_y_continuous("Variable Sites FPR") +
  ggtitle("Better Calibration Increases\nTrue Positive Calls") +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

#and more FPs
pdf("../figures/fp_heatmap.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
  geom_raster(aes(fill = FP)) +
  scale_fill_viridis_c("FP Calls") +
  scale_x_continuous("Variable Sites FNR") +
  scale_y_continuous("Variable Sites FPR") +
  ggtitle("Better Calibration Increases\nFalse Positive Calls") +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

#Shown together
pdf("../figures/tp_fp_plot.pdf", width = 9, height = 7)
df %>% ggplot(aes(TP, FP)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  coord_fixed(ratio = 1) +
  geom_segment(
    aes(x = 267750, y = 7200, xend = 268050, yend = 7200 + (268050-267750))) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis') +
  ggtitle("Better Calibration Produces\nMore Positive Calls") +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0),
        axis.text = element_text(size = rel(.5)))
dev.off()

#This leads to an increase in sensitivity
pdf("../figures/sensitivity.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = recall)) +
  scale_fill_viridis_c('Sensitivity', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Better Calibration Increases Sensitivity') +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

#But a decrease in precision
pdf("../figures/precision.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = precision)) +
  scale_fill_viridis_c('Precision', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Better Calibration Reduces Precision') +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

#Shown together:
pdf("../figures/sens_precision.pdf", width = 9, height = 7)
df %>% ggplot(aes(precision, recall)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  coord_fixed(ratio = 1) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis') +
  ggtitle("Better Calibration Improves\nSensitivity At A Cost") +
  scale_x_continuous("Precision") +
  scale_y_continuous("Sensitivity") +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0),
        axis.text = element_text(size = rel(.5)))
dev.off()
# Taken together, this leads to a reduced F-statistic for the well-calibrated
# data.
# F line plot 
pdf("../figures/f_plot.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalsePositiveRate, f)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  geom_line(aes(color = factor(FalseNegativeRate))) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = "viridis") +
  ggtitle("Better Calibration Yields\nLower Call F-statistic") + 
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()
# F Heatmap
pdf("../figures/f_heatmap.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = f)) +
  scale_fill_viridis_c('F', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Combined Effect on F-statistic') +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()
#Import filtered DF
# --- Import Dataframe ---
import_flt_tsvs <- function(fnr, fpr){
  prefixes <- unlist(map(fnr, paste0, '_', fpr))
  summary_files <- get_summary_path(paste0('flt/',prefixes,'_flt'))
  dfs <- map(summary_files, read_tsv,
             col_names = c("mode","type","operator","value"))
  names(dfs) <- prefixes
  # notempty <- map_lgl(dfs, ~dim(.)[1] != 0)
  # dfs <- dfs[notempty]
  fnr_var <- sym('%FNR')
  fdr_var <- sym('%FDR')
  bind_rows(dfs, .id = 'FNR_FPR') %>%
    pivot_wider(names_from = operator, values_from = value) %>%
    separate(FNR_FPR, into = c("FalseNegativeRate","FalsePositiveRate"),
             sep = "_", convert = TRUE) %>%
    rename(FNR = all_of(fnr_var)) %>%
    rename(FDR = all_of(fdr_var)) %>%
    mutate(FDR = FDR/100) %>%
    mutate(FNR = FNR/100)
}

# -- Main ---
fltdf <- import_flt_tsvs(fnr,fnr) %>%
  calc_f_and_filter()

#The patterns in the data do not significantly change after filtration;
# the better-calibrated data maintains more TPs, more FPs, higher sensitivity
# and lower precision.

#Better calibrated data has more TPs
pdf("../figures/flt_tp_heatmap.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
  geom_raster(aes(fill = TP)) +
  scale_fill_viridis_c("TP Calls") +
  scale_x_continuous("Variable Sites FNR") +
  scale_y_continuous("Variable Sites FPR") +
  ggtitle("Better Calibration Increases True Positive Calls") +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

#and more FPs
pdf("../figures/flt_fp_heatmap.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
  geom_raster(aes(fill = FP)) +
  scale_fill_viridis_c("FP Calls") +
  scale_x_continuous("Variable Sites FNR") +
  scale_y_continuous("Variable Sites FPR") +
  ggtitle("Better Calibration Increases False Positive Calls") +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

#Shown together
pdf("../figures/tp_fp_plot.pdf", width = 10, height = 6)
fltdf %>% ggplot(aes(TP, FP)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  coord_fixed(ratio = 1) +
  # geom_segment(
    # aes(x = 217500, y = 2800, xend = 219750, yend = 2800 + (219750-217500))) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis') +
  ggtitle("Better Calibration Produces More Positive Calls") +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0),
        axis.text = element_text(size = rel(.5)))
dev.off()

#This leads to an increase in sensitivity
pdf("../figures/flt_sensitivity.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = recall)) +
  scale_fill_viridis_c('Sensitivity', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Better Calibration Increases Sensitivity') +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

#But a decrease in precision
pdf("../figures/flt_precision.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = precision)) +
  scale_fill_viridis_c('Precision', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Better Calibration Reduces Precision') +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

#Shown together:
pdf("../figures/flt_sens_precision.pdf", width = 9, height = 7)
fltdf %>% ggplot(aes(precision, recall)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  coord_fixed(ratio = 1) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis') +
  ggtitle("Better Calibration Improves Sensitivity At A Cost") +
  scale_x_continuous("Precision") +
  scale_y_continuous("Sensitivity") +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0),
        axis.text = element_text(size = rel(.5)))
dev.off()

# In contrast to the pre-filtration data, the best-calibrated data now
# has the best F-statstic
pdf("../figures/flt_f_plot.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalsePositiveRate, f)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  geom_line(aes(color = factor(FalseNegativeRate))) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = "viridis") +
  ggtitle("After Filtering,\nBetter Calibration Yields\nHigher Call F-statistic") + 
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

# F Heatmap
pdf("../figures/flt_f_heatmap.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = f)) +
  scale_fill_viridis_c('F', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle("After Filtering, Better Calibration\nYields Higher Call F-statistic") +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))
dev.off()

