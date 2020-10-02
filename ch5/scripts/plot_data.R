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

get_calibration_path <- function(prefix){
  return(paste0(get_data_dir("calibration"),paste0("/", prefix, ".calibration.txt.gz")))
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

import_calibration_df <- function(fnr,fpr){
  prefixes <- unlist(map(fnr, paste0, '_', fpr))
  cal_files <- get_calibration_path(prefixes)
  dfs <- map(cal_files, read_tsv, na = c("","NA","."), col_types = cols(
    QUAL = "d",
    GQ = "i",
    DP = "i",
    QD = "d",
    BaseQRankSum = "d",
    PL = "c",
    CALL = "c"
  ))
  names(dfs) <- prefixes
  # notempty <- map_lgl(dfs, ~dim(.)[1] != 0)
  # dfs <- dfs[notempty]
  bind_rows(dfs, .id = "FNR_FPR") %>%
    separate(FNR_FPR, into = c("FalseNegativeRate","FalsePositiveRate"),
             sep = "_", convert = TRUE) %>%
    separate(PL, into = c("PL1","PL2","PL3"), sep = ",", convert = TRUE) %>%
    mutate(PLSum = PL1 + PL2 + PL3) %>%
    select(-PL1, -PL2, -PL3)
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

### --- ABOVE THIS LINE IS OLD - FROM RTG DFs

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
# Qual is −10log10prob(no variant), so we can actually make a calibration plot
# ~= -10log10 P(FP)
# FP_CA = false positive GT, but the call contains the correct ALT allele
caldf <- import_calibration_df(fnr,fnr)

#turn a probability to a phred-scaled quality score
p_to_q <- function(p, maxscore = 2000){
  return(if_else(p != 0, floor(-10*log10(p)), maxscore))
}

q_to_p <- function(q){
  return(10 ** (q / -10))
}

#given phred-scaled QUAL values, average them in probability space and return
# the values in PHRED space
calc_binned_qual <- function(q){
  p <- q_to_p(q)
  p_avg <- mean(p)
  q_avg <- p_to_q(p_avg)
}

#a histogram for QUAL given TP; shows a slight shift but not too large
caldf %>% group_by(FalseNegativeRate, FalsePositiveRate) %>%
  filter(CALL == "TP") %>%
  ggplot(aes(QUAL)) + geom_freqpoly(bins = 100, aes(color = factor(FalseNegativeRate))) + 
  facet_grid(cols = vars(FalsePositiveRate), as.table = FALSE) + 
  scale_color_viridis_d()

#GQ is basically useless wrt a real "Q" value. It's better than QUAL though.
caldf %>% group_by(FalsePositiveRate,FalseNegativeRate,GQ) %>%
  summarize(realq = p_to_q(sum(CALL != "TP")/length(CALL))) %>%
  filter(realq != 2000) %>%
  ggplot(aes(GQ,realq)) +
  geom_line(aes(color = factor(.group))) +
  geom_abline(slope = 1)
# 

cs <- caldf %>% mutate(BQRP = 2*pnorm(abs(BaseQRankSum), lower.tail = FALSE)) %>% group_by(FalseNegativeRate,FalsePositiveRate) %>% arrange(desc(BQRP), .by_group = TRUE) %>% mutate(cumFP = cumsum(CALL == "FP"), cumTP = cumsum(CALL == "TP"))
cs %>% slice(seq(1,n(),1000)) %>% ggplot(aes(cumFP,cumTP)) + geom_line(aes(color = factor(.group))) + scale_color_viridis_d()

#So it seems most of these distributions are pretty useless, and QUAL is not
# really related to the definition of QUAL. 
calc_roc <- function(df, variable){
  df %>%
    arrange(desc({{variable}}), .by_group = TRUE) %>%
    mutate(cumFP = cumsum(CALL != "TP"), cumTP = cumsum(CALL == "TP"))
}

plot_roc <- function(df){
  df %>%
    slice(seq(1,n(),1000)) %>%
    ggplot(aes(cumFP,cumTP)) +
    facet_grid(rows = vars(FalsePositiveRate), as.table = FALSE) +
    geom_line(aes(color = factor(FalseNegativeRate))) +
    scale_color_viridis_d(option = 'viridis')
}

calc_roc_auc <- function(df){
  df %>% summarize(AUC = sum(cumTP))
}

plot_roc_auc <- function(df){
  df %>%
    ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
    geom_raster(aes(fill = AUC)) +
    scale_fill_viridis_c()
}

#ROCs on ROCs on ROCs
gcaldf <- caldf %>% 
  group_by(FalseNegativeRate,FalsePositiveRate) %>%
  filter(FalsePositiveRate != 100)

#qualroc
qualroc <- gcaldf %>% calc_roc(QUAL)
qualroc_auc <- qualroc %>% calc_roc_auc()
qualroc %>% plot_roc()
qualroc_auc %>% plot_roc_auc()

#gqroc
gqroc <- gcaldf %>% calc_roc(GQ)
gqroc_auc <- gqroc %>% calc_roc_auc()
gqroc %>% plot_roc()
gqroc_auc %>% plot_roc_auc()

#dproc
dproc <- gcaldf %>% calc_roc(DP)
dproc_auc <- dproc %>% calc_roc_auc()
dproc %>% plot_roc()
dproc_auc %>% plot_roc_auc()

#qd
qdroc <- gcaldf %>% calc_roc(QD)
qdroc_auc <- qdroc %>% calc_roc_auc()
qdroc %>% plot_roc()
qdroc_auc %>% plot_roc_auc()

#pl
plroc <- gcaldf %>% calc_roc(PLSum)
plroc_auc <- plroc %>% calc_roc_auc()
plroc %>% plot_roc()
plroc_auc %>% plot_roc_auc()

#The less calibrated data simply has fewer TPs
gcaldf %>% summarise(TP = sum(CALL=="TP")) %>% ggplot(aes(FalseNegativeRate,FalsePositiveRate)) + geom_raster(aes(fill = TP)) + scale_fill_viridis_c()

#AND fewer FPs
gcaldf %>% summarise(FP = sum(CALL!="TP")) %>% ggplot(aes(FalseNegativeRate,FalsePositiveRate)) + geom_raster(aes(fill = FP)) + scale_fill_viridis_c()


