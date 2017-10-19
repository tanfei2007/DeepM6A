rm(list=ls())
library(ggplot2)
library(plotROC)

#******************
#plot roc curves
#******************
roc_plot <- function(filename_rslt, filename_auc, len = 30,  species){
  
  #read prediction results
  files = list.files(filename_rslt, pattern="*.csv")
  files = paste(filename_rslt, files, sep = "")

  
  for(f in files) {
    if(!exists("df_rslt")) {
      df_rslt = read.csv(f, header = T, stringsAsFactors = F)
      df_rslt$index = substr(f, 1, nchar(f)-4)
    } else {
      temp = read.csv(f, header = T, stringsAsFactors = F)
      temp$index = substr(f, 1, nchar(f)-4)
      df_rslt = rbind(df_rslt, temp)
    }
  }
  colnames(df_rslt) <- c('pred', 'groundtruth', 'index')
  
  #read auc results
  auc = read.csv(filename_auc, header = T, stringsAsFactors = F)
  auc = subset(auc, length == len, select = auc)
  
  basicplot = ggplot(df_rslt, aes(d = groundtruth, m = pred, color = index)) + geom_roc(n.cuts = 0, size=0.2) + 
              style_roc(xlab = "False Positive Rate", ylab = "True Positive Rate") +
              theme(legend.position="none") + scale_colour_manual(values=colorRampPalette(c("darkgrey", "black"))(10))
  ROCplot = basicplot + ggtitle(species) + 
            annotate("text", x = .7, y = .2, 
                     label = paste("Average AUC = ", round(mean(auc$auc), 3), '(', round(sd(auc$auc), 3) ,')', sep = ""))
  
  rm(list = 'df_rslt')
  return(ROCplot)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#************************
# plot different species
#************************


#main paper
Cele = roc_plot(filename_rslt = 'Celegan_rslt30/', filename_auc = './C.elegan.csv',
                len = 30, species = 'C. elegans')

Ara = roc_plot(filename_rslt = 'Ara_rslt30/', filename_auc = './Ara.csv',
               len = 30, species = 'A. thaliana')

Dr = roc_plot(filename_rslt = 'Dr3_rslt30/', filename_auc = './Dr3.csv',
              len = 30, species = 'D. melanogaster')

pdf("ROC_Species_Main.pdf", width = 12, height = 3)
multiplot(Cele, Ara, Dr, cols=3)
dev.off()

#supplementary materials

# Yeast = roc_plot(filename_rslt = 'Yeast_rslt30/', filename_auc = './Yeast.csv', 
#                  len = 30, species = 'Yeast')
# 
# Toly = roc_plot(filename_rslt = 'Toly_rslt30/', filename_auc = './Toly.csv', 
#                 len = 30, species = 'Tolypocladium ')

Ecoli = roc_plot(filename_rslt = 'Ecoli_rslt30/', filename_auc = './Ecoli.csv', 
                 len = 30, species = 'E. coli')

Zebrafish = roc_plot(filename_rslt = 'Zebrafish_rslt30/', filename_auc = './Zebrafish.csv', 
                 len = 30, species = 'Zebrafish')

pdf("ROC_Species_Supp.pdf", width = 8, height = 3)
multiplot(Ecoli, Zebrafish, cols=2)
#multiplot(Ecoli, cols=1)
dev.off()



