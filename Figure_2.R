library(ggplot2)
library(reshape)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

load("Figure_2.Rdata")
dats$Pops <- factor(dats$Pops, levels = c("SAS", "EAS", "AFR"))
qs <- list()
##Plot for Scenario 1
dats1 <- dats[dats$Scenario == "Scenario1",]
dat1 <- melt(dats1, id = c("Pops", "Scenario"))
qs[[1]] <- ggplot(data = dat1, aes(x = Pops, y = value), group = variable) +
  geom_boxplot(position=position_dodge(0.8), aes(color = variable), outlier.alpha = 0.5, outlier.size = 0.5, width = 0.7) +
  geom_point(position=position_dodge(0.8), aes(color = variable), size = 0.5, alpha = 0.5) +
  theme_classic() +
  theme(
    text=element_text(size=16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5,vjust = 1.2),
    legend.position = c(0.22,0.35),
    legend.text = element_text(size = 14),
    legend.key.height=unit(1.2,"line"),
    legend.key.width=unit(1.2,"line"),
    legend.background = element_blank(),
    legend.title = element_blank(),
    axis.line=element_blank(),
    legend.key = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  scale_y_continuous(expand = c(0,0)) +
  geom_segment( aes(x=0.2, xend=0.2, y=-0.0005, yend= 100.0005), size = 0.5) +
  geom_segment( aes(x=0.995, xend=3.005, y=-0.0003, yend=-0.0003), size = 0.5) +  
  scale_color_manual(values = cbPalette ) +
  ylab(paste0("Relative Accuracy (%)")) +
  scale_x_discrete(labels = c( "SAS", "EAS", "AFR")) 


##Plot for Scenario 2 & 3
for(ss in 2:3){
  dats1 <- dats[dats$Scenario == paste0("Scenario", ss),]
  dat1 <- melt(dats1, id = c("Pops", "Scenario"))
  qs[[ss]] <- ggplot(data = dat1, aes(x = Pops, y = value), group = variable) +
    geom_boxplot(position=position_dodge(0.8), aes(color = variable), outlier.alpha = 0.5, outlier.size = 0.5, width = 0.7) +
    geom_point(position=position_dodge(0.8), aes(color = variable), size = 0.5, alpha = 0.5) +
    theme_classic() +
    theme(
      text=element_text(size=16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5,vjust = 1.2),
      legend.position = "none",
      axis.line=element_blank(),
      axis.title.x = element_blank(),
      axis.ticks = element_line(size = 0.5),
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    scale_y_continuous(expand = c(0,0)) +
    geom_segment( aes(x=0.2, xend=0.2, y=-0.0005, yend= 100.0005), size = 0.5) +
    geom_segment( aes(x=0.995, xend=3.005, y=-0.0003, yend=-0.0003), size = 0.5) +  
    scale_color_manual(values = cbPalette) +
    ylab(paste0("Relative Accuracy (%)")) +
    scale_x_discrete(labels = c( "SAS", "EAS", "AFR")) 
}

les <- c("A", "B", "C")
for(ll in 1:3){
  qs[[ll]] <-  arrangeGrob(qs[[ll]], top = textGrob(les[ll], x = unit(0, "npc")  , y   = unit(0.8, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=16))) 
}
ggsave("Figure_2.pdf", grid.arrange(qs[[1]], qs[[2]], qs[[3]], nrow = 1), height = 4, width = 14, dpi = 300)
