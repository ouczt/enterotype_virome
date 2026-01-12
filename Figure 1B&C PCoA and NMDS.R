#Fig1 B&C genus level

# Create base PCoA plot
create_pcoa_plot <- function(data, eig_percent) {
  ggplot(data, aes(x = PCoA1, y = PCoA2, color = enterotype)) +
    geom_point(aes(shape = enterotype), size = 4, alpha = 0.5) +
    labs(x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
         y = paste("PCoA 2 (", eig_percent[2], "%)", sep = "")) +
    scale_colour_manual(values = c("#1773B6", "#972735")) +
    theme_minimal() +
    theme(
      legend.position = c(0.9, 0.69),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(color = 'black', fill = 'transparent'),
      axis.text = element_text(color = "black", size = 10)
    ) +
    geom_hline(yintercept = 0, colour = "#BEBEBE", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "#BEBEBE", linetype = "dashed")
}

# Create basic plot
p1 <- create_pcoa_plot(pcoa_result, eig_percent)

# Add ellipses
p2 <- p1 + 
  stat_ellipse(geom = "polygon", 
               level = 0.9, 
               linetype = 2, 
               linewidth = 0.5, 
               aes(fill = enterotype), 
               alpha = 0.3, 
               show.legend = TRUE) +
  scale_fill_manual(values = c("#a1d5b9", "#e1abbc"))
library(ggExtra)
# Create marginal plots
marginal_plot <- ggMarginal(
  p2,
  type = c('density'),
  margins = 'both',
  size = 3.5,
  groupColour = FALSE,
  groupFill = TRUE
)
marginal_plot
# Perform PERMANOVA
permanova_result <- adonis2(t_data ~ Cluster, 
                            data = pcoa_result, 
                            permutations = 999, 
                            method = "bray")

# Print results
print(permanova_result)
print(marginal_plot)

# Save plots if needed
ggsave("pcoa_plot.pdf", p2, width = 10, height = 8)
ggsave("E:/Tao_Zuo_Lab/enterotype manuscript/revision/FigS2/pcoa_plot_with_margins.pdf", marginal_plot, width = 12, height = 10)



#NMDS
p1 <- ggplot(data = nmds_result, aes(x=MDS1, y=MDS2, color=enterotype)) +
  geom_point(aes(color=enterotype,shape=enterotype),size=5,alpha=0.5)+
  labs(x=paste("NMDS 1"),
       y=paste("NMDS 2"),
       caption= paste('Stress =', stress))  +# 也可用title、caption
  scale_colour_manual(values = c("#1773B6","#972735"))+
  theme(legend.position = c(0.9,0.8),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.title=element_text(hjust=0),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.text = element_text(color = "black",size=10))+
  geom_hline(aes(yintercept=0), colour="#BEBEBE", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="#BEBEBE", linetype="dashed")


p2=p1+stat_ellipse(data=nmds_result,
                       geom="polygon",
                       level=0.99,
                       linetype=2,
                       linewidth=0.5,
                       aes(fill=enterotype),
                       alpha=0.3,
                       show.legend=T)+
  scale_fill_manual(values=c("#a1d5b9","#e1abbc"))#这里添加或修改颜色

p3 <-ggMarginal(
  p2,
  type=c('boxplot'),
  margins='both',
  size=3.5,
  groupColour=F,
  groupFill=T
)

rownames(t_data) <- nmds_result$Sample_ID

# Compute Euclidean distance (as in your code)
eucl_dist <- vegdist(t_data, method = "euclidean")

# PERMANOVA test (enterotype ~ community structure)
permanova_result <- adonis2(eucl_dist ~ enterotype, 
                            data = nmds_result, 
                            permutations = 999)

# Print results (R² = explained variance, Pr(>F) = p-value)
print(permanova_result)#R²=0.0863 p=0.001

