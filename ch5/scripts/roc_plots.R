    #!/usr/bin/env Rscript

#ROCs on ROCs on ROCs!!

library(tidyverse)

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
  .rows = function(str){str_replace(str,"^80$","FPR:\n80")}
)

# --- Make plots ---

r_factor <- function(v){factor({{v}}, levels = c("kbbq","raw",0,20,40,60,80,100))}

plot_roc <- function(df){
  just_empirical <- df %>%
    ungroup() %>%
    filter(FalseNegativeRate == "kbbq" | FalseNegativeRate == "raw") %>%
    select(FalseNegativeRate, FP, TPc)
  
  df %>%
    unite("FNR_FPR", c(FalseNegativeRate,FalsePositiveRate), remove = FALSE) %>%
    mutate(FalsePositiveRate = factor(FalsePositiveRate)) %>%
    mutate(FalsePositiveRate = fct_collapse(FalsePositiveRate,
      empirical = c('kbbq','raw')
    )) %>%
    ggplot(aes(FP,TPc)) +
    facet_grid(rows = vars(FalsePositiveRate), as.table = FALSE, 
               margins = TRUE, labeller = rate_labeller) +
    # geom_line(aes(color = factor(FalseNegativeRate, levels = c("raw",0,20,40,60,80,100)), group = FNR_FPR)) +
    geom_line(aes(color = FalseNegativeRate, group = FNR_FPR)) +
    scale_color_viridis_d('Variable\nSites\nFNR', option = 'viridis') +
    xlab("False Positive Calls") +
    ylab("True Positive Calls") +
    theme_minimal(base_size = 18) + 
    theme(plot.margin = margin(0,0,0,0),
          strip.text.y = element_text(angle = 0)) +
    geom_line(data = just_empirical, aes(color = FalseNegativeRate, group = FalseNegativeRate)) #add empirical on top of every plot
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
    scale_fill_viridis_c()
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
    theme_minimal(base_size = 18) + 
    theme(plot.margin = margin(0,0,0,0))
}

plot_roc_avgf <- function(df){
  df %>%
    summarize(avgf = median(f)) %>%
    ggplot(aes(FalseNegativeRate,FalsePositiveRate)) +
    geom_raster(aes(fill = avgf)) +
    scale_fill_viridis_c() + 
    xlab("Variable Sites FNR") +
    ylab("Variable Sites FPR") +
    theme_minimal(base_size = 18) + 
    theme(plot.margin = margin(0,0,0,0))
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

# I plot the ROC, argmaxf, and avgf

#qualroc
qualroc <- import_rtg_rocs(fpr,fpr,"QUAL")

pdf("../figures/qualroc.pdf", width = 9, height = 7)
qualroc %>% plot_roc() + xlim(0,3000) +
  ggtitle("Better Calibration Makes QUAL A Better Classifier")
dev.off()

pdf("../figures/qualmaxf.pdf", width = 9, height = 7)
qualroc %>% plot_roc_argmaxf(QUAL) +
  ggtitle("Better Calibration Makes Optimal QUAL Filter Larger")
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
gqroc <- import_rtg_rocs(fpr,fpr,"GQ")
pdf("../figures/gqroc.pdf", width = 9, height = 7)
gqroc %>% plot_roc() +
  ggtitle("Better Calibration Reduces Classifying Power of GQ")
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



