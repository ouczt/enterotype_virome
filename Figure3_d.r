#figure 3d
# Ensure mantel has the categorized columns rd and pd
mantel <- mantel %>%
  mutate(
    rd = cut(r,
             breaks = c(-Inf, 0.1, 0.2, Inf),
             labels = c("<= 0.1", "0.1 to 0.2", ">= 0.2")),
    pd = cut(p,
             breaks = c(-Inf, 0.05, Inf),
             labels = c("<0.05", ">= 0.05"))
  )

# Create the plot
qcorrplot(rcorr(cor_E2_subset), type = "upper", 
          diag = FALSE, grid_size = 0.4, grid_col = "lightgray") +
  geom_square(linetype = 0) +
  geom_couple(data = mantel, aes(colour = pd, size = rd),
              node.shape = 21, node.size = 2, node.colour = "#0000ff", node.fill = "#9ba9f1",
              curvature = 0.1) + 
  scale_fill_gradientn(limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.5),
                       colours = c("#5BC2CD", "white", "#FF8040")) +
  scale_size_manual(values = c(0.2, 0.5, 1.0),  # 3 values for 3 levels of rd
                    breaks = c("<= 0.1", "0.1 to 0.2", ">= 0.2")) +
  scale_colour_manual(values = c("#87CEEB", "gray"),  # 2 values for 2 levels of pd
                      breaks = c("<0.05", ">= 0.05")) +
  guides(size = guide_legend(title = expression("Mantel's" ~ italic("r")),
                             override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = expression("Mantel's" ~ italic("p")),
                               override.aes = list(size = 2), order = 1),
         fill = guide_colorbar(title = expression("Pearson's" ~ italic("r")), order = 3)) +
  geom_mark(only_mark = T, size = 2, sig_level = c(0.05, 0.01, 0.001), sig_thres = 0.05, colour = "black") +
  theme(axis.text.x = element_text(size = 6),  # Adjust x-axis text size
        axis.text.y = element_text(size = 6))  # Adjust y-axis text size