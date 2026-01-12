# Load libraries
library(phyloseq)
library(vegan)
library(BiodiversityR)
library(tidyverse)

# Set working directory
setwd("E:/Tao_Zuo_Lab/enterotype script")

# -------------------- Load data --------------------
species_RPKM <- read.table("otu_virus.txt", header = TRUE, row.names = 1, sep = "\t") %>%
  t() %>% as.data.frame()

# Transpose for analysis (samples as rows)
species_RPKM_t <- t(species_RPKM)

# -------------------- Diversity indices --------------------
diversity <- tibble(
  shannon = diversity(species_RPKM_t, index = "shannon"),
  simpson = diversity(species_RPKM_t, index = "simpson"),
  sample  = rownames(species_RPKM_t)
)

# -------------------- Richness (Chao1) --------------------
species_counts <- round(species_RPKM * 10000) %>% as.data.frame()
otu_tab <- otu_table(species_counts, taxa_are_rows = TRUE)
chao1 <- estimate_richness(otu_tab, measures = "Chao1") %>%
  rownames_to_column("sample")

# -------------------- Metadata --------------------
subject_enterotype <- read.csv("metadata.csv") %>%
  mutate(
    enterotype = recode(Cluster, `1` = "Enterotype 2", `2` = "Enterotype 1"),
    enterotype = factor(enterotype, levels = c("Enterotype 1", "Enterotype 2"))
  )

# -------------------- Merge diversity + richness + metadata --------------------
merged_diversity <- diversity %>%
  left_join(subject_enterotype, by = c("sample" = "Sample")) %>%
  left_join(chao1, by = "sample")

# -------------------- Evenness --------------------
evenness <- diversityresult(species_RPKM_t, index = "Jevenness", method = "each site") %>%
  as.data.frame() %>%
  rownames_to_column("sample")

# Merge evenness with metadata
evenness2 <- evenness %>%
  left_join(subject_enterotype, by = c("sample" = "Sample"))

# -------------------- Statistical tests --------------------
# Richness
t.test(Chao1 ~ enterotype, data = merged_diversity)#p=0.83711
wilcox.test(Chao1 ~ enterotype, data = merged_diversity)#p=0.89778

# Shannon
t.test(shannon ~ enterotype, data = merged_diversity)#p-value < 2.22e-16
wilcox.test(shannon ~ enterotype, data = merged_diversity)#p-value < 2.22e-16

# Simpson
t.test(simpson ~ enterotype, data = merged_diversity)#p-value < 2.22e-16
wilcox.test(simpson ~ enterotype, data = merged_diversity)#p-value < 2.22e-16

# Evenness
t.test(Jevenness ~ enterotype, data = evenness2)#p-value < 2.22e-16
wilcox.test(Jevenness ~ enterotype, data = evenness2)#p-value < 2.22e-16

library(ggplot2)
library(ggpubr)

# Define colors
cols <- c("#1773B6","#972735")

# -------------------- Shannon --------------------
p_shannon <- ggplot(merged_diversity, aes(x = enterotype, y = shannon, color = enterotype)) +
  geom_jitter(shape = 19, size = 5, alpha = 0.6, width = 0.3) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.colour = NA, color = "black",
               fill = alpha(cols, 0.4)) +
  scale_color_manual(values = cols) +
  theme_bw() +
  labs(y = "Shannon diversity", x = "") +
  ylim(0, 3.0) +
  stat_compare_means(method = "wilcox.test", label = "p.format")
p_shannon
# -------------------- Simpson --------------------
p_simpson <- ggplot(merged_diversity, aes(x = enterotype, y = simpson, color = enterotype)) +
  geom_jitter(shape = 19, size = 5, alpha = 0.6, width = 0.3) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.colour = NA, color = "black",
               fill = alpha(cols, 0.4)) +
  scale_color_manual(values = cols) +
  theme_bw() +
  labs(y = "Simpson diversity", x = "") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
p_simpson
# -------------------- Chao1 --------------------
p_Richness <- ggplot(merged_diversity, aes(x = enterotype, y = Chao1, color = enterotype)) +
  geom_jitter(shape = 19, size = 5, alpha = 0.6, width = 0.3) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.colour = NA, color = "black",
               fill = alpha(cols, 0.4)) +
  scale_color_manual(values = cols) +
  theme_bw() +
  labs(y = "Richness", x = "") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
p_Richness
# -------------------- Evenness --------------------
evenness <- read.table("evenness.txt", header = TRUE, sep = "\t")
evenness2 <- merge(evenness, merged_diversity[,c("sample","enterotype")],
                   by.x = "Row.names", by.y = "sample")
colnames(evenness2)[2] <- "Jevenness"

p_evenness <- ggplot(evenness2, aes(x = enterotype, y = Jevenness, color = enterotype)) +
  geom_jitter(shape = 19, size = 5, alpha = 0.6, width = 0.3) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.colour = NA, color = "black",
               fill = alpha(cols, 0.4)) +
  scale_color_manual(values = cols) +
  theme_bw() +
  labs(y = "Evenness (J)", x = "") +
  stat_compare_means(method = "wilcox.test", label = "p.format")

# -------------------- Print plots --------------------
print(p_shannon)
print(p_simpson)
print(p_chao1)
print(p_evenness)

# -------------------- Statistics --------------------
wilcox.test(shannon ~ enterotype, data = merged_diversity)
wilcox.test(simpson ~ enterotype, data = merged_diversity)
wilcox.test(Chao1 ~ enterotype, data = merged_diversity)
wilcox.test(Jevenness ~ enterotype, data = evenness2)

library(ggplot2)
library(ggpubr)

cols <- c("#1773B6","#972735")

# Function to make consistent plots
make_plot <- function(df, yvar, ylab_text, ylim_max){
  ggplot(df, aes(x = enterotype, y = !!sym(yvar), color = enterotype)) +
    geom_jitter(shape = 19, size = 6, alpha = 0.6, width = 0.3) +
    geom_boxplot(alpha = 0.5, width = 0.5, outlier.colour = NA, color = "black",
                 fill = alpha(cols, 0.4)) +
    scale_color_manual(values = cols) +
    labs(y = ylab_text, x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    # Add extra space above for p-value bars
    expand_limits(y = ylim_max) +
    # Clean theme: no grey grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 16))
}

# Shannon
p_shannon <- make_plot(merged_diversity, "shannon", "Shannon diversity", max(merged_diversity$shannon)*1.2)

# Simpson
p_simpson <- make_plot(merged_diversity, "simpson", "Simpson diversity", max(merged_diversity$simpson)*1.2)

# Chao1
p_chao1 <- make_plot(merged_diversity, "Chao1", "Chao1 richness", max(merged_diversity$Chao1)*1.2)

# Evenness
evenness <- read.table("evenness.txt", header = TRUE, sep = "\t")
evenness2 <- merge(evenness, merged_diversity[,c("sample","enterotype")],
                   by.x = "sample", by.y = "sample")
colnames(evenness2)[2] <- "Jevenness"

p_evenness <- make_plot(evenness2, "Jevenness", "Evenness (J)", max(evenness2$Jevenness)*1.2)

# Print plots
print(p_shannon)
print(p_simpson)
print(p_chao1)
print(p_evenness)

