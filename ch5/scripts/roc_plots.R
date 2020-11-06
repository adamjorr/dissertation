#!/usr/bin/env Rscript

#ROCs on ROCs on ROCs!!

library(tidyverse)
library(colorblindr)

# --- Import Data ---
#import rtg-eval roc files.

get_data_dir <- function(path){
  return(paste0("../data/roc/", path))
}

#return file name of the roc tsv
#given the prefix and type. For example, 0_0 QUAL
#returns "../data/roc/0_0_QUAL.snp_roc.tsv.gz"
get_roc_path <- function(prefix, type){
  return(get_data_dir(paste0(prefix, "_", type, ".snp_roc.tsv.gz")))
}

# --- Get Data Frame ---
import_rtg_rocs <- function(fnr, fpr, type){
  prefixes <- unlist(map(fnr, paste0, '_', fpr))
  prefixes <- c(prefixes, "raw","cbbq")
  roc_files <- get_roc_path(prefixes, type)
  dfs <- map(roc_files, read_tsv, comment = "#",
             col_names = c(
               type,"TPb","FP","TPc", "FN", "precision", "sensitivity", "f"))
  prefixes[(length(prefixes)-1):length(prefixes)] = c("raw_raw","kbbq_kbbq")
  names(dfs) <- prefixes
  notempty <- map_lgl(dfs, ~dim(.)[1] != 0)
  dfs <- dfs[notempty]
  bind_rows(dfs, .id = 'FNR_FPR') %>%
  separate(FNR_FPR, into = c("FalseNegativeRate","FalsePositiveRate"),
           sep = "_", convert = TRUE) %>%
  mutate(across(c(FalseNegativeRate,FalsePositiveRate),
                ~factor(., levels = c("raw","kbbq",0,20,40,60,80,100)))) %>%
  group_by(FalseNegativeRate,FalsePositiveRate) %>%
  filter(FalsePositiveRate != 100)
}

rate_labeller <- labeller(
  # .cols = function(str){str_replace(str,"^","FNR: 0")},
  # .rows = function(str){str_replace(str,'^(\\d+)$','FPR:\\1')}
  .rows = function(str){str_replace(str,'^0$','FPR:0')}
)

# --- Make plots ---

r_factor <- function(v){factor({{v}}, levels = c("kbbq","raw",0,20,40,60,80,100))}

plot_roc <- function(df){
  just_empirical <- df %>%
    ungroup() %>%
    filter(FalseNegativeRate == "kbbq" | FalseNegativeRate == "raw") %>%
    select(FalseNegativeRate, FPR, sensitivity)
  
  df %>%
    unite("FNR_FPR", c(FalseNegativeRate,FalsePositiveRate), remove = FALSE) %>%
    mutate(FalsePositiveRate = factor(FalsePositiveRate)) %>%
    mutate(FalsePositiveRate = fct_collapse(FalsePositiveRate,
      empirical = c('kbbq','raw')
    )) %>%
    ggplot(aes(FPR,sensitivity)) +
    # facet_grid(rows = vars(FalsePositiveRate), as.table = FALSE, 
               # margins = TRUE, labeller = rate_labeller) +
    facet_wrap(facets = vars(FalsePositiveRate), as.table = TRUE, 
              labeller = rate_labeller) +
    geom_line(aes(color = factor(FalseNegativeRate, levels = c("raw","kbbq",0,20,40,60,80,100)), group = FNR_FPR)) +
    # geom_line(aes(color = FalseNegativeRate, group = FNR_FPR)) +
    scale_color_viridis_d('Variable\nSites\nFNR', option = 'viridis') +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_fixed(ratio = max(df$FPR)) +
    geom_line(data = just_empirical, aes(color = FalseNegativeRate, group = FalseNegativeRate)) #add empirical on top of every plot
}

plot_single_roc <- function(df){
  df %>%
    unite("FNR_FPR", c(FalseNegativeRate,FalsePositiveRate), remove = FALSE) %>%
    mutate(FalsePositiveRate = factor(FalsePositiveRate)) %>%
    mutate(FalsePositiveRate = fct_collapse(FalsePositiveRate,
                                            empirical = c('kbbq','raw')
    )) %>%
    ggplot(aes(FPR,sensitivity)) +
    geom_line(aes(color = factor(FalseNegativeRate, levels = c("raw","kbbq",0,20,40,60,80,100)), group = FNR_FPR)) +
    scale_color_viridis_d('Variable\nSites\nFNR', option = 'viridis') +
    # scale_color_OkabeIto(name = 'Variable\n Sites\nFNR') +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    theme(strip.text.y = element_text(angle = 0)) +
    coord_fixed(ratio = max(df$FPR))
}

plot_sen_prec <- function(df){
  df %>% ggplot(aes(sensitivity,precision)) +
    facet_grid(rows = vars(FalsePositiveRate), as.table = FALSE) +
    geom_line(aes(color = factor(FalseNegativeRate))) +
    scale_color_viridis_d(option = 'viridis')
}

plot_roc_maxf <- function(df){
  df %>%
    summarize(maxf = max(f)) %>%
    ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
    geom_raster(aes(fill = maxf)) +
    scale_fill_viridis_c() + 
    coord_fixed(ratio = 1)
}

#the value of the statistic where F is maximized
plot_roc_argmaxf <- function(df, var){
  df %>%
    slice_max(order_by = f, with_ties = FALSE) %>%
    ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
    geom_raster(aes(fill = {{var}})) +
    scale_fill_viridis_c() +
    xlab("Variable Sites FNR") +
    ylab("Variable Sites FPR") +
    coord_fixed(ratio = 1)
}

plot_roc_avgf <- function(df){
  df %>%
    summarize(avgf = median(f)) %>%
    ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
    geom_raster(aes(fill = avgf)) +
    scale_fill_viridis_c() + 
    xlab("Variable Sites FNR") +
    ylab("Variable Sites FPR") +
    coord_fixed(ratio = 1)
}

#emulate filtering on an (exclusive) maximum and minimum value.
# The ROC lists the cumulative TP, FP, FN in DESCENDING ORDER.
#
# Thus, above the maximum:
# FPs and TPs are removed.
# FNs need to be added.
#
# Below the minimum:
# FPs and TPs are removed.
# FNs are calculated correctly.
#
# Precision, sensitivity, and F are then recalculated.
emulate_filter <- function(df, var, min, max){
  ceil <- df %>%
    filter({{var}} > max, .preserve = TRUE) %>%
    slice_min(order_by = {{var}}, with_ties = FALSE) %>%
    select(FalseNegativeRate,FalsePositiveRate,TPb:FN) %>%
    rename(TPbMAX = TPb, FPMAX = FP, TPcMAX = TPc) %>%
    mutate(FNMAX = max(FN)) %>%
    select(-FN)
  
  df %>% select(FalseNegativeRate:FN) %>%
    filter({{var}} <= max, .preserve = TRUE) %>%
    left_join(ceil) %>%
    mutate(FP = FP - FPMAX,
           TPb = TPb - TPbMAX,
           TPc = TPc - TPcMAX) %>%
    filter({{var}} >= min, .preserve = TRUE) %>%
    mutate(FN = if_else(FN == max(FN),FNMAX,FN)) %>%
    select(!TPbMAX:FNMAX) %>%
    mutate(sensitivity = TPb/(TPb + FN), precision = TPc/(TPc + FP)) %>%
    mutate(f = 2 * precision * sensitivity / (precision + sensitivity))
}


# --- Main ---
fpr <- seq(0,100,20)

#zcat ansites.bed.gz | awk '{x += ($3-$2)}END{print x}'
neg_sites <- 213684509

#set theme
theme_set(theme_minimal(base_size = 22, base_family = 'Times') +
            theme(plot.margin = margin(0,0,0,0))
          )

# I plot the ROC, argmaxf, and avgf

#qualroc
qualroc <- import_rtg_rocs(fpr,fpr,"QUAL") %>% mutate(FPR = FP/neg_sites)

pdf("../figures/qualroc.pdf", width = 9, height = 7)
qualroc %>% plot_roc() + #xlim(0 , 2e-5) +
  ggtitle("QUAL ROC")
dev.off()

pdf("../figures/0fpr_roc.pdf", width = 9, height = 7)
qualroc %>% filter(FalsePositiveRate == 0 | FalsePositiveRate == "kbbq" | FalsePositiveRate == "raw") %>% plot_single_roc() +
  ggtitle("QUAL ROC")
dev.off()

    pdf("../figures/qualmaxf.pdf", width = 9, height = 7)
qualroc %>% plot_roc_argmaxf(QUAL) +
  ggtitle("QUAL with Best F-statistic")
dev.off()
# qualroc %>% plot_roc_avgf() +
#   ggtitle("Better Calibration Increases Average QUAL F-statistic")

#dproc
# I decided not to use this since it's not well-defined as a score:
#   larger numbers != better. Depth can probably be transformed to
#   a score assuming a Poisson distribution, but that is outside of
#   the scope of this work.

# dproc <- import_rtg_rocs(fpr,fpr,"DP")
# dproc %>% plot_roc()
# dproc %>% plot_roc_maxf()

#gqroc
gqroc <- import_rtg_rocs(fpr,fpr,"GQ") %>% mutate(FPR = FP/neg_sites)
pdf("../figures/gqroc.pdf", width = 9, height = 7)
gqroc %>% plot_roc() +
  ggtitle("GQ ROC")
dev.off()
# gqroc %>% plot_roc_argmaxf(GQ) +
#   ggtitle("")
# gqroc %>% plot_roc_avgf()

# Better calibration reduces the power of GQ, but increases the power of QUAL
bcftools_data_dir <- function(path){
  return(paste0("../data/bcftools_roc/", path))
}

bcftools_roc_path <- function(prefix, type){
  return(bcftools_data_dir(paste0(prefix, "_", type, ".snp_roc.tsv.gz")))
}

import_bcftools_rocs <- function(fnr, fpr, type){
  prefixes <- unlist(map(fnr, paste0, '_', fpr))
  prefixes <- c(prefixes, "raw","cbbq")
  roc_files <- bcftools_roc_path(prefixes, type)
  dfs <- map(roc_files, read_tsv, comment = "#",
             col_names = c(
               type,"TPb","FP","TPc", "FN", "precision", "sensitivity", "f"))
  prefixes[(length(prefixes)-1):length(prefixes)] = c("raw_raw","kbbq_kbbq")
  names(dfs) <- prefixes
  notempty <- map_lgl(dfs, ~dim(.)[1] != 0)
  dfs <- dfs[notempty]
  bind_rows(dfs, .id = 'FNR_FPR') %>%
    separate(FNR_FPR, into = c("FalseNegativeRate","FalsePositiveRate"),
             sep = "_", convert = TRUE) %>%
    mutate(across(c(FalseNegativeRate,FalsePositiveRate),
                  ~factor(., levels = c("raw","kbbq",0,20,40,60,80,100)))) %>%
    group_by(FalseNegativeRate,FalsePositiveRate) %>%
    filter(FalsePositiveRate != 100)
}



