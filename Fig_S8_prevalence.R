library(dplyr)
library(tidyr)
library(readr)

# Path to your metrics files
files <- list.files("E:/Tao_Zuo_Lab/linux_analyse_output/bwa-mem2 for phage abundance/output", 
                    pattern = "_metrics.tsv$", full.names = TRUE)

# Read all files and add sample name
metrics_list <- lapply(files, function(f) {
  sample <- gsub("_metrics.tsv$", "", basename(f))   # strip suffix
  df <- read_tsv(f, col_types = cols())
  df <- df %>% select(genome, breadth) %>% mutate(Sample = sample)
  return(df)
})

# Combine into one long table
metrics_long <- bind_rows(metrics_list)

# Pivot to wide format: rows=genome, cols=Sample, values=breadth
metrics_wide <- metrics_long %>%
  pivot_wider(
    names_from = Sample,
    values_from = breadth,
    values_fill = list(breadth = 0)   # replace NA with 0
  )

# Write to file
write_tsv(metrics_wide, "phage_breadth_matrix.tsv")

library(tidyverse)
df <- metrics_wide   # 你的数据
threshold <- 0.30   # breadth > 30% 认为存在

# 计算每个噬菌体的检出率（prevalence）
prev <- df %>%
  pivot_longer(-genome, names_to = "Sample", values_to = "breadth") %>%
  mutate(present = breadth > threshold) %>%
  group_by(genome) %>%
  summarise(
    prevalence = mean(present),      # 检出率
    mean_breadth = mean(breadth),    # 平均覆盖度（备用排序）
    .groups = "drop"
  ) %>%
  arrange(desc(prevalence), desc(mean_breadth))

# 取出 Top 10（严格按照检出率排序）
top10 <- prev %>% slice_max(prevalence, n = 10) %>% pull(genome)

# 只保留 Top 10
plot_data <- df %>%
  filter(genome %in% top10) %>%
  pivot_longer(-genome, names_to = "Sample", values_to = "breadth") %>%
  mutate(
    present = breadth > threshold,
    genome = factor(genome, levels = rev(top10))  # 垂直顺序：检出率从下到上递增
  )

# 样本排序：可按组排序（如果有分组信息就用，没有就按样本名）
# 如果你有分组（比如 GZ-CD vs GZ-HC），取消下面注释并改成你的分组列名
# plot_data <- plot_data %>%
#   left_join(metadata %>% select(Sample, Group), by = "Sample") %>%
#   arrange(Group, Sample) %>%
#   mutate(Sample = factor(Sample, levels = unique(Sample)))

# 没分组就直接按名字排序
plot_data <- plot_data %>%
  arrange(Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# ====================================
p <- ggplot(plot_data, aes(x = Sample, y = genome)) +
  geom_tile(aes(fill = present), color = "white", size = 0.2, height = 0.85, width = 0.95) +
  scale_fill_manual(
    values = c("FALSE" = "gray88", "TRUE" = "#D7191C"),  # 经典 Nature 红
    labels = c("Absent", "Present (breadth > 30%)"),
    name = NULL
  ) +
  labs(
    title = "Distribution of Top 10 Most Prevalent GPIC Phages in Human Gut",
    subtitle = paste0("n = ", n_distinct(plot_data$Sample), 
                      " samples | Red = genome breadth coverage > 30%"),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5),
    axis.text.y = element_text(face = "italic", size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    panel.grid = element_blank(),
    legend.position = "top",
    plot.margin = margin(15, 15, 15, 15)
  )

# ===================================
#ggsave("Top10_GPIC_Phage_Distribution.png", p, width = 20, height = 7, dpi = 400, bg = "white")

#ggsave("Top10_GPIC_Phage_Distribution.pdf", p, width = 20, height = 7, device = cairo_pdf)

print(p)

library(patchwork)   

 