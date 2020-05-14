library(ggplot2)
library(reshape)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

load("Figure_1.Rdata")

dat2 <- melt(dats, id = c("Pops", "h2", "M"))

p <- ggplot(dat = dat2, aes(x = Pops, y = value), group = variable) +
  geom_boxplot(position=position_dodge(0.8), aes(color = variable), outlier.alpha = 0.5, outlier.size = 0.5, width = 0.7) +
  geom_point(position=position_dodge(0.8), aes(color = variable), size = 0.5, alpha = 0.5) +
  facet_grid(h2 ~ M, labeller = label_bquote(cols = italic(M[C]) *"="* .(M), rows =  italic(h)^2 *"="* .(h2))) +
  theme_bw() +
  theme(
    text=element_text(size=16),
    plot.title = element_text(hjust = 0.5,vjust = 1.2),
    legend.position = c(0.06,0.65),
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
  scale_color_manual(values = cbPalette )+
  ylab("Relative Accuracy (%)") +
  scale_x_discrete(labels = c( "SAS \n(0.0234)", "EAS \n(0.1046)", "AFR \n(0.1360)"))

ggsave("Figure_1.pdf", p, height = 8, width = 14, dpi = 300)

