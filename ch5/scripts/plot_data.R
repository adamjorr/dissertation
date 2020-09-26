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
  df %>% only_snps() %>% get_f_statistic() %>% filter(f > .05)
}

# --- Plotting ---
  
df %>%
  unite("FNR_FPR", c(FalseNegativeRate,FalsePositiveRate)) %>%
  ggplot(aes(false_positives,true_positives_call)) +
  geom_line()

df %>% group_by(FalseNegativeRate,FalsePositiveRate) %>%
  add_max_fscore() %>%
  ggplot(aes(false_positives, true_positives_call, color = factor(FalseNegativeRate))) +
  geom_line() +
  coord_fixed(ratio = 1) +
  facet_grid(cols = vars(FalsePositiveRate), rows = NULL,
             switch = "both",
             as.table = FALSE) +
  scale_color_viridis_d('FNR', option = "viridis")

df %>% group_by(FalseNegativeRate,FalsePositiveRate) %>%
  add_max_fscore() %>%
  ggplot(aes(sensitivity, precision, color = factor(FalsePositiveRate))) +
  geom_line() +
  facet_grid(cols = vars(FalseNegativeRate), rows = NULL,
             switch = "both",
             as.table = FALSE) +
  scale_color_viridis_d(option = "viridis")

# F line plot
df %>% calc_f_and_filter() %>%
  ggplot(aes(FalseNegativeRate, f)) +
  geom_point(aes(color = factor(FalsePositiveRate))) +
  geom_line(aes(color = factor(FalsePositiveRate))) +
  scale_color_viridis_d("FPR", option = "viridis")

# F Heatmap
df %>% calc_f_and_filter() %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
    geom_raster(aes(fill = f)) +
    scale_fill_viridis_c('F', option = "viridis") + 
    xlab("False Negative Rate") +
    ylab("False Positive Rate") +
    ggtitle('Combined Effect on F-statistic') +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 18) + 
    theme(plot.margin = margin(0,0,0,0))

#Recall Heatmap (= sensitivity)
#Better calibration *slightly* increases sensitivity!
df %>% calc_f_and_filter() %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = recall)) +
  scale_fill_viridis_c('Recall', option = "viridis") + 
  xlab("False Negative Rate") +
  ylab("False Positive Rate") +
  ggtitle('Combined Effect on Recall') +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))

#Precision
#Precision INCREASES with worse calibration! That is,
#a higher proportion of variants are correct with worse calibration!
df %>% calc_f_and_filter() %>%
  ggplot(aes(FalseNegativeRate, FalsePositiveRate)) +
  geom_raster(aes(fill = precision)) +
  scale_fill_viridis_c('Precision', option = "viridis") + 
  xlab("False Negative Rate") +
  ylab("False Positive Rate") +
  ggtitle('Combined Effect on Precision') +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) + 
  theme(plot.margin = margin(0,0,0,0))

#Need to investigate DP, Qual, and QD. ROCs for Qual and QD seem promising.
# I predict proper scores will have an outsized effect at low DP, so try that.
# I want to look into BaseQRankSum more. It's not straightforward to plot the ROC
# since it's a Z score, but we should be able to filter on some arbitrary +/- threshold
# and see if we see an effect.
# At the end, let's combine all the filters we find and see what the combined
# effect of them all is.
# I may need to plot some distributions based on the annotations, like distribution
# of Depth given TP, FP and FN and compare these.

