rm(list=ls())
library(stringr)
library(openxlsx)
library(ggplot2)
library(gridExtra)


plot_mcc = function(file1, file2, file3, file4, output){

	# Ara
	ara.LR.mcc.dat = read.xlsx(file1)
	ara.LR.mcc.nb.dat = read.xlsx(file2)
	ara.DP.mcc.dat = read.csv(file3)
	ara.DP.mcc.nb.dat = read.csv(file4)

	ara.LR.mcc.table = data.frame(LR.mcc=c(ara.LR.mcc.dat$len7, ara.LR.mcc.nb.dat$len7), 
	                          distance=rep(c('200', 'Closest'), each=10), 
	                          index=rep(1:10, 2))
	ara.DP.mcc.table = data.frame(DP.mcc=c(ara.DP.mcc.dat$mcc, ara.DP.mcc.nb.dat$mcc), 
	                          distance=rep(c('200', 'Closest'), each=10), 
	                          index=rep(1:10, 2))

	ara.mcc.table = merge(ara.LR.mcc.table, ara.DP.mcc.table, by=c("distance","index"))


	ara.mcc.table1 = data.frame(MCC=c(ara.mcc.table$LR.mcc, ara.mcc.table$DP.mcc),
	                        distance=rep(ara.mcc.table$distance, 2),
	                        method=rep(c("LR", "DeepM6A"), each=nrow(ara.mcc.table)))
	ara.mcc.table1$method = factor(ara.mcc.table1$method, c('DeepM6A', 'LR'))



	ara.p = ggplot(ara.mcc.table1, aes(x=method, y=MCC, fill=distance)) + geom_boxplot(position=position_dodge(1)) +
	  theme(legend.title=element_blank(), panel.background = element_rect(fill = 'white', colour = 'black'), 
	        panel.grid.major = element_line(colour = "gray90"), panel.grid.minor = element_line(colour = "gray90"), 
	        plot.title = element_text(hjust = 0.5, size=20),
	        legend.key = element_rect(fill = 'white', colour = 'white'), legend.position="bottom") +
	  ylim(0,1) + scale_fill_manual(values=c("deepskyblue", "green")) +
	  stat_summary(fun.y="mean", geom="point", size=1, shape=5, position=position_dodge(1), aes(colour=distance)) +
	  scale_colour_manual(values = c("plum1", "red")) +
	  ggtitle(output) + xlab("")

	return(ara.p)  
}

ara.p = plot_mcc("LR/Ara_length_rep10_mcc.xlsx", "LR/Ara_neighbour_rep10_mcc.xlsx", "DL/Ara.csv", "DL/Ara_nb.csv", output = "A. thaliana")
dr.p = plot_mcc("LR/Dr3_length_rep10_mcc.xlsx", "LR/Dr3_neighbour_rep10_mcc.xlsx", "DL/Dr3.csv", "DL/Dr3_nb.csv", output = "D. melanogaster")
ecoli.p = plot_mcc("LR/Ecoli_length_rep10_mcc.xlsx", "LR/Ecoli_neighbour_rep10_mcc.xlsx", "DL/Ecoli.csv", "DL/Ecoli_nb.csv", output = "E. coli")



#require(cowplot)
pdf("nb_MCC.pdf", width=12, height=4, onefile=TRUE)
#multiplot(ara.p, dr.p, ecoli.p, cols=3)
#plot_grid(ara.p, dr.p, ecoli.p, labels = c('A', 'B', 'C'))
grid.arrange(ara.p, dr.p, ecoli.p, ncol=3, nrow=1)
dev.off()








