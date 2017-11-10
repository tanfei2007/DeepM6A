rm(list=ls())
library(stringr)
library(openxlsx)
library(ggplot2)
library(gridExtra)


plot_auc = function(file1, file2, file3, file4, output){

	# Ara
	ara.LR.auc.dat = read.xlsx(file1)
	ara.LR.auc.nb.dat = read.xlsx(file2)
	ara.DP.auc.dat = read.csv(file3)
	ara.DP.auc.nb.dat = read.csv(file4)

	ara.LR.auc.table = data.frame(LR.AUC=c(ara.LR.auc.dat$len7, ara.LR.auc.nb.dat$len7), 
	                          distance=rep(c('200', 'Closest'), each=10), 
	                          index=rep(1:10, 2))
	ara.DP.auc.table = data.frame(DP.AUC=c(ara.DP.auc.dat$auc, ara.DP.auc.nb.dat$auc), 
	                          distance=rep(c('200', 'Closest'), each=10), 
	                          index=rep(1:10, 2))

	ara.auc.table = merge(ara.LR.auc.table, ara.DP.auc.table, by=c("distance","index"))


	ara.auc.table1 = data.frame(AUC=c(ara.auc.table$LR.AUC, ara.auc.table$DP.AUC),
	                        distance=rep(ara.auc.table$distance, 2),
	                        method=rep(c("LR", "DeepM6A"), each=nrow(ara.auc.table)))
	ara.auc.table1$method = factor(ara.auc.table1$method, c('DeepM6A', 'LR'))



	ara.p = ggplot(ara.auc.table1, aes(x=method, y=AUC, fill=distance)) + geom_boxplot(position=position_dodge(1)) +
	  theme(legend.title=element_blank(), panel.background = element_rect(fill = 'white', colour = 'black'), 
	        panel.grid.major = element_line(colour = "gray90"), panel.grid.minor = element_line(colour = "gray90"), 
	        plot.title = element_text(hjust = 0.5, size=20),
	        legend.key = element_rect(fill = 'white', colour = 'white'), legend.position="bottom") +
	  ylim(0.5,1) + scale_fill_manual(values=c("deepskyblue", "green")) +
	  stat_summary(fun.y="mean", geom="point", size=1, shape=5, position=position_dodge(1), aes(colour=distance)) +
	  scale_colour_manual(values = c("plum1", "red")) +
	  ggtitle(output) + xlab("")

	return(ara.p)  
}

ara.p = plot_auc("LR/Ara_length_rep10_auc.xlsx", "LR/Ara_neighbour_rep10_auc.xlsx", "DL/Ara.csv", "DL/Ara_nb.csv", output = "A. thaliana")
dr.p = plot_auc("LR/Dr3_length_rep10_auc.xlsx", "LR/Dr3_neighbour_rep10_auc.xlsx", "DL/Dr3.csv", "DL/Dr3_nb.csv", output = "D. melanogaster")
ecoli.p = plot_auc("LR/Ecoli_length_rep10_auc.xlsx", "LR/Ecoli_neighbour_rep10_auc.xlsx", "DL/Ecoli.csv", "DL/Ecoli_nb.csv", output = "E. coli")



#require(cowplot)
pdf("nb_AUC.pdf", width=12, height=4, onefile=TRUE)
#multiplot(ara.p, dr.p, ecoli.p, cols=3)
#plot_grid(ara.p, dr.p, ecoli.p, labels = c('A', 'B', 'C'))
grid.arrange(ara.p, dr.p, ecoli.p, ncol=3, nrow=1)
dev.off()








