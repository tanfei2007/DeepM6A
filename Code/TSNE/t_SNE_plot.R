rm(list=ls())
library(ggplot2)
setwd("~/Documents/Research/m6A/Result/t_SNE/Celegan/")
dat = read.csv("layer_dropout_6_seq_tsne2.csv", stringsAsFactors=F, check.names=F)
dat$label[dat$label==0] = "Negative"
dat$label[dat$label==1] = "Positive"
dat$label = factor(dat$label, levels=c("Positive", "Negative"))
dat$fraction = "None"
dat$fraction[0<dat$frac & dat$frac<=0.2] = "Lowly"
dat$fraction[0.2<dat$frac & dat$frac<=0.8] = "Intermediate"
dat$fraction[0.8<dat$frac] = "Highly"
pdf("Celegan_t_SNE.pdf", width=6, height=5, onefile=F)
ggplot(dat, aes(x=`x-tsne`, y=`y-tsne`, color=fraction)) + 
  geom_point(alpha=0.8, size=0.5) + 
  ggtitle("tSNE dimensions colored by 6mA") +
  theme_bw() +
  scale_color_manual(values=c("blue", "pink", "red", "darkred"), limits=c("None", "Lowly", "Intermediate", "Highly"))
dev.off()