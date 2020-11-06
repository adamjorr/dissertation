library(tidyverse)
library(colorblindr)

rg_val <- 24
assigned_q_val <- 30
covs <- c("Context", "Cycle")
color_levels <- c("Read Group", "Assigned Q", "Context", "Cycle", "Recalibrated")

data_df <- data.frame(
  group = factor(c(1, 1, 2, 2)),
  x = c(.9,1.1,1.9,2.1),
  covariate = factor(c(covs,covs), levels = color_levels),
  quality = c(31, 28, 25, 32),
  rg_qual = rg_val,
  assigned_qual = assigned_q_val
)

label_df <- data.frame(label = c("Read Group Estimate","Assigned Q Estimate for this RG"),
                      # group = factor(c(3)),
                      color = c("Read Group","Assigned Q"),
                      quality = c(rg_val, assigned_q_val)
                      )

label2_df <- data.frame(label = c("Context Estimate for this RG and Q","Cycle Estimate for this RG and Q"),
                        color = c("Context","Cycle"),
                        quality = c(25,31.5))

adding_df <- data.frame(label = c("30 + 1 - 2","30 + 2 - 5"),
                        color = c("black","black"),
                        x = c(1,2),
                        y = c(29, 27))

pdf("./figures/recalibration_explainer.pdf", width = 11, height = 7)
#show how to go from assigned Q to a new quality score
data_df %>% ggplot(aes(group, quality, family = 'Times')) +
  # geom_pointrange(aes(ymin = assigned_qual, ymax = quality), size = 1) +
  geom_segment(aes(x = x, xend = x, y = assigned_q_val, yend = quality, color = covariate), size = 1, arrow = arrow(length = unit(.15,"inches"))) +
  geom_text(mapping = aes(x = 3.5, quality+.25, color = color, label = label, hjust = 1),
            data = label_df, show.legend = F, size = 6) +
  geom_text(mapping = aes(x = 3.5, quality+.25, color = color, label = label, hjust = 1), data = label2_df, show.legend = FALSE, size = 6) +
  # geom_text(mapping = aes(x = 3.5, y = 25.25, color = "Context", label = "Context Estimate", hjust = 1), data = data.frame()) +
  # geom_text(mapping = aes(x = 3.5, y = 31.75, color = "Cycle", label = "Cycle Estimate", hjust = 1)) +
  geom_hline(aes(yintercept = quality, color = color), data = label_df, size = 1) +
  geom_segment(aes(x = .5, xend = .5, y = rg_qual, yend = assigned_qual,
                   color = "Read Group"),
               arrow = arrow(length = unit(.15,"inches")), show.legend = FALSE, size = 1) +
  # geom_curve(aes(x = .9, xend = .95, y = assigned_qual, yend = 29.05),
  #              arrow = arrow(length = unit(.15,"inches")), show.legend = FALSE, curvature = .5, size = 1, color = "black") +
  # geom_curve(aes(x = 1.9, xend = 1.95, y = assigned_qual, yend = 27.05),
  #              arrow = arrow(length = unit(.15,"inches")), show.legend = FALSE, curvature = .5, size = 1, color = "black") +
  geom_segment(aes(x = 1, xend = 1, y = assigned_qual, yend = 29.075),
             arrow = arrow(length = unit(.15,"inches")), show.legend = FALSE, size = 1, color = "black") +
  geom_segment(aes(x = 2, xend = 2, y = assigned_qual, yend = 27.075),
             arrow = arrow(length = unit(.15,"inches")), size = 1, color = "black") +
  geom_point(x = 1, y = 29, color = "black", size = 3) +
  geom_point(x = 2, y = 27, color = "black", size = 3) +
  geom_text(mapping = aes(x = x + .15, y = y, label = label), color = "black", data = adding_df, hjust = 0, show.legend = FALSE, size = 6) +
  xlab("Base") + ggtitle("Two Bases Assigned The Same Quality") + ylab("Quality") +
  scale_color_OkabeIto(name = '', use_black = T, order = c(1,2,3,7,8), drop = FALSE)

dev.off()
