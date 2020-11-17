#!/usr/bin/env Rscript
library(tidyverse)
library(colorblindr)

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
  prefixes <- c(prefixes, "raw","cbbq")
  summary_files <- get_summary_path(prefixes)
  dfs <- map(summary_files, read_tsv,
             col_names = c("mode","type","operator","value"))
  prefixes[(length(prefixes)-1):length(prefixes)] = c("raw_raw","kbbq_kbbq")
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
    mutate(FNR = FNR/100) %>%
    mutate(across(c(FalseNegativeRate,FalsePositiveRate),
          ~factor(., levels = c("raw","kbbq",0,20,40,60,80,100))))
}

import_bcftools_tsvs <- function(fnr, fpr){
  prefixes <- unlist(map(fnr, paste0, '_', fpr))
  prefixes <- c(prefixes, "raw","cbbq")
  summary_files <- get_summary_path(paste0("bcftools/",prefixes))
  dfs <- map(summary_files, read_tsv,
             col_names = c("mode","type","operator","value"))
  prefixes[(length(prefixes)-1):length(prefixes)] = c("raw_raw","kbbq_kbbq")
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
    mutate(FNR = FNR/100) %>%
    mutate(across(c(FalseNegativeRate,FalsePositiveRate),
                  ~factor(., levels = c("raw","kbbq",0,20,40,60,80,100))))
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

#set theme
theme_set(theme_minimal(base_size = 22, base_family = 'Times') +
            theme(plot.margin = margin(0,0,0,0))
)

# --- Plotting ---

#Better calibrated data has more TPs
pdf("../figures/tp_heatmap.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
  geom_raster(aes(fill = TP)) +
  scale_fill_viridis_c("TP Calls") +
  scale_x_discrete("Variable Sites FNR") +
  scale_y_discrete("Variable Sites FPR") +
  # ggtitle("Better Calibration Increases\nTrue Positive Calls") +
  ggtitle("Unfiltered True Positive Calls")
dev.off()

#and more FPs
pdf("../figures/fp_heatmap.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
  geom_raster(aes(fill = FP)) +
  scale_fill_viridis_c("FP Calls") +
  scale_x_discrete("Variable Sites FNR") +
  scale_y_discrete("Variable Sites FPR") +
  ggtitle("Unfiltered False Positive Calls")
  # ggtitle("Better Calibration Increases\nFalse Positive Calls") +
dev.off()

#Shown together
pdf("../figures/tp_fp_plot.pdf", width = 9, height = 7)
df %>% ggplot(aes(TP, FP)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  # coord_fixed(ratio = 1) +
  # geom_segment(
  #   aes(x = 267750, y = 7200, xend = 268050, yend = 7200 + (268050-267750))) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis') +
  ggtitle("Unfiltered Positive Calls") +
  # ggtitle("Better Calibration Produces\nMore Positive Calls") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#This leads to an increase in sensitivity
pdf("../figures/sensitivity.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = recall)) +
  scale_fill_viridis_c('Sensitivity', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  # ggtitle('Better Calibration Increases Sensitivity') +
  ggtitle('Unfiltered Sensitivity') +
  coord_fixed(ratio = 1)
dev.off()

#But a decrease in precision
pdf("../figures/precision.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = precision)) +
  scale_fill_viridis_c('Precision', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Unfiltered Precision') +
  coord_fixed(ratio = 1)
dev.off()

#Shown together:
pdf("../figures/sens_precision.pdf", width = 9, height = 7)
df %>% ggplot(aes(precision, recall)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  # coord_fixed(ratio = 1) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis') +
  ggtitle("Unfiltered Sensitivity and Precision") +
  scale_x_continuous("Precision") +
  scale_y_continuous("Sensitivity")
dev.off()
# Taken together, this leads to a reduced F-statistic for the well-calibrated
# data.
# F line plot 
pdf("../figures/f_plot.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalsePositiveRate, f)) +
  geom_point(aes(color = FalseNegativeRate)) +
  geom_line(aes(color = FalseNegativeRate, group = FalseNegativeRate)) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = "viridis") +
  ggtitle("Unfiltered F-statistic")
dev.off()
# F Heatmap
pdf("../figures/f_heatmap.pdf", width = 9, height = 7)
df %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = f)) +
  scale_fill_viridis_c('F', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Unfiltered F-statistic') +
  coord_fixed(ratio = 1)
dev.off()
#Import filtered DF
# --- Import Dataframe ---
import_flt_tsvs <- function(fnr, fpr){
  prefixes <- unlist(map(fnr, paste0, '_', fpr))
  prefixes <- c(prefixes, "raw","cbbq")
  summary_files <- get_summary_path(paste0('flt/',prefixes,'_flt'))
  dfs <- map(summary_files, read_tsv,
             col_names = c("mode","type","operator","value"))
  prefixes[(length(prefixes)-1):length(prefixes)] = c("raw_raw","kbbq_kbbq")
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
    mutate(FNR = FNR/100) %>% 
    mutate(across(c(FalseNegativeRate,FalsePositiveRate),
                  ~factor(., levels = c("raw","kbbq",0,20,40,60,80,100))))
}

import_bcftools_flt_tsvs <- function(fnr, fpr){
  prefixes <- unlist(map(fnr, paste0, '_', fpr))
  prefixes <- c(prefixes, "raw","cbbq")
  summary_files <- get_summary_path(paste0('bcftools_flt/',prefixes,'_flt'))
  dfs <- map(summary_files, read_tsv,
             col_names = c("mode","type","operator","value"))
  prefixes[(length(prefixes)-1):length(prefixes)] = c("raw_raw","kbbq_kbbq")
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
    mutate(FNR = FNR/100) %>% 
    mutate(across(c(FalseNegativeRate,FalsePositiveRate),
                  ~factor(., levels = c("raw","kbbq",0,20,40,60,80,100))))
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
  scale_x_discrete("Variable Sites FNR") +
  scale_y_discrete("Variable Sites FPR") +
  ggtitle("Filtered True Positive Calls") +
  coord_fixed(ratio = 1)
dev.off()

#and more FPs
pdf("../figures/flt_fp_heatmap.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
  geom_raster(aes(fill = FP)) +
  scale_fill_viridis_c("FP Calls") +
  scale_x_discrete("Variable Sites FNR") +
  scale_y_discrete("Variable Sites FPR") +
  ggtitle("Filtered False Positive Calls") +
  coord_fixed(ratio = 1)
dev.off()

#Shown together
pdf("../figures/flt_tp_fp_plot.pdf", width = 10, height = 6)
fltdf %>% ggplot(aes(TP, FP)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis') +
  ggtitle("Filtered Positive Calls")
dev.off()

#This leads to an increase in sensitivity
pdf("../figures/flt_sensitivity.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = recall)) +
  scale_fill_viridis_c('Sensitivity', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Filtered Sensitivity') +
  coord_fixed(ratio = 1)
dev.off()

#But a decrease in precision
pdf("../figures/flt_precision.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = precision)) +
  scale_fill_viridis_c('Precision', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle('Filtered Precision') +
  coord_fixed(ratio = 1)
dev.off()

#Shown together:
pdf("../figures/flt_sens_precision.pdf", width = 9, height = 7)
fltdf %>% ggplot(aes(precision, recall)) +
  geom_point(aes(color = factor(FalseNegativeRate))) +
  # coord_fixed(ratio = 1) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis') +
  ggtitle("Filtered Sensitivity And Precision") +
  scale_x_continuous("Precision") +
  scale_y_continuous("Sensitivity")
dev.off()

#Show before/after filtering of sensitivity and precision:
bothdf <- bind_rows(before = df, after = fltdf, .id = "filtering") %>%
  pivot_wider(id_cols = c(FalseNegativeRate, FalsePositiveRate),
              names_from = filtering, values_from = c(recall, precision))
bothdf %>% ggplot(aes(precision_before, recall_before)) +
  geom_segment(aes(xend = precision_after, yend = recall_after,
                   color = factor(FalseNegativeRate)),
               arrow = arrow()) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = 'viridis')

# In contrast to the pre-filtration data, the best-calibrated data now
# has the best F-statstic
pdf("../figures/flt_f_plot.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalsePositiveRate, f)) +
  geom_point(aes(color = FalseNegativeRate)) +
  geom_line(aes(color = FalseNegativeRate, group = FalseNegativeRate)) +
  scale_color_viridis_d("Variable\nSites\nFNR", option = "viridis") +
  ggtitle("Filtered F-statistic")
dev.off()

# F Heatmap
pdf("../figures/flt_f_heatmap.pdf", width = 9, height = 7)
fltdf %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = f)) +
  scale_fill_viridis_c('F', option = "viridis") + 
  xlab("Variable Sites FNR") +
  ylab("Variable Sites FPR") +
  ggtitle("Filtered F-statistic") +
  coord_fixed(ratio = 1)
dev.off()

# ---- Simulated Data Summaries -----
get_sim_f <- function(df){
  df %>%
    group_by(CalibrationMethod) %>%
    mutate(recall = TPc / (TPc + FN)) %>%
    mutate(precision = TPc / (TPc + FP)) %>%
    mutate(f = 2 * precision * recall / (precision + recall))
}

import_sim_tsvs <- function(){
  datanames <- c("ngm","ngm.recal","initial-calls.recal","kbbq-ngm.recal")
  nicenames <- c("Raw","GATK","Initial-calls","KBBQ")
  summary_files <- paste0("../data/sims/",datanames,".vcf.gz.summary.txt")
  
  dfs <- map(summary_files, read_tsv,
             col_names = c("mode","type","operator","value"))
  names(dfs) <- nicenames
  
  fnr_var <- sym('%FNR')
  fdr_var <- sym('%FDR')
  bind_rows(dfs, .id = 'CalibrationMethod') %>%
    pivot_wider(names_from = operator, values_from = value) %>%
    rename(FNR = all_of(fnr_var)) %>%
    rename(FDR = all_of(fdr_var)) %>%
    mutate(FDR = FDR/100) %>%
    mutate(FNR = FNR/100) %>%
    mutate(CalibrationMethod = factor(CalibrationMethod, levels = nicenames)) %>%
    only_snps() %>%
    get_sim_f()
}

sim_tsvs <- import_sim_tsvs()
pdf("../figures/sims_sens_precision.pdf", width = 9, height = 7)
sim_tsvs %>% ggplot(aes(precision, recall)) +
  geom_point(aes(color = CalibrationMethod), size = 3) +
  ggtitle("Sensitivity and Precision of Calls from Simulated Data") +
  scale_x_continuous("Precision") +
  scale_y_continuous("Sensitivity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_OkabeIto(name = 'Calibration Method', use_black = T, drop = FALSE)
dev.off()

pdf("../figures/sims_fstat.pdf", width = 9, height = 7)
sim_tsvs %>%
  ggplot(aes(CalibrationMethod, f)) +
  geom_point(aes(color = CalibrationMethod), size = 3) +
  ylab("F-statistic") +
  ggtitle("F-statistic of Calls from Simulated Data") +
  scale_color_OkabeIto(name = 'Calibration Method', use_black = T, drop = FALSE)
dev.off()

sim_table <- sim_tsvs %>%
  select(CalibrationMethod, TPc, FP, recall, precision,f) %>%
  rename(TP = TPc, Sensitivity = recall, Precision = precision, "F-statistic" = f)

write_tsv(sim_table, "../tables/sims_summary.txt")
