#figure 3g-i
library(ggplot2)
library(ggpubr) # 用于stat_compare_means函数

# 假设dat已经按照enterotype和CD_HC进行了分组和排序
dat <- meta_bacteria %>% 
  group_by(enterotype, CD_HC) %>%
  arrange(desc(CD_HC))

# 将 Blautia_obeum 列转换为数值型
dat$Blautia_obeum <- as.numeric(as.character(dat$Blautia_obeum))
dat$Blautia_obeum <-  log10(dat$Blautia_obeum+0.0001)
dat$Bacteroides_thetaiotaomicron <- as.numeric(as.character(dat$Bacteroides_thetaiotaomicron))
dat$Bacteroides_thetaiotaomicron <-  log10(dat$Bacteroides_thetaiotaomicron+0.0001)

# 生成 Group 列，假设 enterotype 和 CD_HC 的组合形成了 Group
dat$Group <- with(dat, paste0(enterotype, "_", CD_HC))
dat$Group <- factor(dat$Group, levels = c("1_HC", "1_CD", "2_HC", "2_CD"))
# 创建ggplot图表
p <- ggplot(dat, aes(x=Group, y=Blautia_obeum, color=Group)) +
  geom_jitter(alpha=0.6, size=5, 
              position=position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 0.01)) +
  geom_boxplot(alpha=0.8, width=0.5,
               position=position_dodge(width=2),
               linewidth=0.8, outlier.colour = NA, color="black", fill=alpha(c("#1773B6","#1773B6","#972735","#972735"))) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) + # 外边框粗细
  theme(legend.position="none") +
  theme(text = element_text(size=1)) +
  scale_fill_manual(values=alpha(c('#3CAF20','#F9874E','#0A5A60','#F73D0E'), 0.8)) +
  scale_color_manual(values=c("#1773B6","#1773B6","#972735","#972735")) +
  ylab("Abundance of Wulfhauvirus") +
  xlab("") +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA), 
        panel.grid.minor = element_blank(),) +
  theme(plot.title = element_text(hjust = 0.5, size = 30, face="bold")) +
  theme(axis.text = element_text(size = 20, color="black"), 
        axis.ticks.x=element_line(color="black", linewidth=1), 
        axis.line.x=element_line(linetype=1, color="black", linewidth=1),
        axis.ticks.y=element_line(color="black", linewidth=1), 
        axis.line.y=element_line(linetype=1, color="black", linewidth=1),
        axis.title= element_text(size=20),
        legend.position = 'none')


print(p)

p <- ggplot(dat, aes(x=Group, y=Bacteroides_thetaiotaomicron, color=Group, fill=Group)) +
  geom_violin(alpha=0.4, width=0.8, color="black", linewidth=0.8, scale = "width", trim = FALSE) +
  geom_jitter(alpha=0.6, size=5, 
              position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  theme(legend.position="none") +
  scale_fill_manual(values=alpha(c("#1773B6","#1773B6","#972735","#972735"), 0.3)) +
  scale_color_manual(values=c("#1773B6","#1773B6","#972735","#972735")) +
  ylab("Abundance of Wulfhauvirus") +
  xlab("") +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA), 
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 30, face="bold")) +
  theme(axis.text = element_text(size = 20, color="black"), 
        axis.ticks.x=element_line(color="black", linewidth=1), 
        axis.line.x=element_line(linetype=1, color="black", linewidth=1),
        axis.ticks.y=element_line(color="black", linewidth=1), 
        axis.line.y=element_line(linetype=1, color="black", linewidth=1),
        axis.title= element_text(size=20))

p <- p + 
  stat_compare_means(comparisons = list(c("1_HC", "1_CD"), c("2_HC", "2_CD")),
                     method = "wilcox.test", 
                     label = "p.signif", # or "p.format"
                     tip.length = 0.01, 
                     size = 6)
print(p)

