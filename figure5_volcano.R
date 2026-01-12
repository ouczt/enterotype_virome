p <- ggplot(data=df,
            aes(x=log2FC,
                y=-log10(pvalue)))+  
  geom_point(alpha=0.6, aes(size = -log10(pvalue), color=log2FC))+
  scale_color_gradientn(colours = c("#2166ac", "#4393c3","#ffffbf", "#de77ae", "#c51b7d"),
                        values = seq(0, 1, 0.2)) + 
  geom_point(data = df %>% 
               dplyr::filter(Name %in% c("rutC", "rpsM", "fusA")),  # Filter for specific genes
             aes(x = log2FC, y = -log10(pvalue), size = -log10(pvalue)), 
             shape = 21, color = "#000000", fill = "#ff7f00") + 
  geom_text_repel(data = df %>% 
                    dplyr::filter(Name %in% c("rutC", "rpsM", "fusA")),  # Filter for specific genes
                  aes(x = log2FC, y = -log10(pvalue), label = Name),
                  nudge_x = 0.5,
                  box.padding = 0.5,
                  nudge_y = 1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20
  ) + 
  scale_size(range = c(2,10),
             guide = guide_legend(override.aes = list(fill = NA)))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) + 
  xlim(c(-3, 3)) + 
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) + 
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) + 
  xlab('log2 fold change')+
  ylab('-log10 pvalue')+
  theme_bw() + 
  theme(plot.title = element_text(size = 15,hjust = 0.5),
        legend.background = element_roundrect(color = '#808080',linetype = 1),
        axis.text = element_text(size = 12.5, color = "#000000"),
        axis.title = element_text(size = 15, color = "#000000")
  ) + 
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#74add1")
    ), 
    xmin = -1.5, 
    xmax = -3,
    ymin = 135,
    ymax = 142
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Down\n(7)",
      gp = grid::gpar(col = "#74add1")
    ),
    xmin = -1.5, 
    xmax = -3,
    ymin = 135,
    ymax = 142
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
      gp = grid::gpar(lwd = 3, col = "#d73027")
    ), 
    xmin = 1.5, 
    xmax = 3,
    ymin = 135,
    ymax = 142
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Up\n(18)",
      gp = grid::gpar(col = "#d73027")
    ),
    xmin = 1.5, 
    xmax = 3,
    ymin = 135,
    ymax = 142
  ) 
print (p)
ggsave(filename="Output/violent_plot.pdf",
       plot=p,
       height=8,
       width=8.5)

library(ggplot2)  # Ensure ggplot2 package is loaded

# Assuming you want to create a basic scatter plot
p <- ggplot(data = df, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(alpha = 0.6, aes(size = -log10(pvalue), color = log2FC)) +
  scale_color_gradientn(colours = c("#2166ac", "#4393c3", "#ffffbf", "#de77ae", "#c51b7d"),
                        values = seq(0, 1, 0.2)) +
  xlab('log2 fold change') +
  ylab('-log10 p-value') +
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        legend.background = element_rect(color = '#808080', linetype = 1),
        axis.text = element_text(size = 12.5, color = "#000000"),
        axis.title = element_text(size = 15, color = "#000000"))

# Print or save the plot
print(p)

# If you want to save it to a PDF file
ggsave(filename = "Output/volcano_plot.pdf",
       plot = p,
       height = 8,
       width = 8.5)