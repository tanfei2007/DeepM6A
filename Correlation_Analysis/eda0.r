rm(list = ls())

#load data
meth_level <- read.csv("meth_level.csv", header = T, stringsAsFactors = F)
meth_level <- as.numeric(meth_level$X0)
idx <- which(!is.na(meth_level))
meth_level <- meth_level[idx]
  
meth_prob <- read.csv("meth_prob.csv", header = T, stringsAsFactors = F)
meth_prob <- meth_prob$pred_prob
meth_prob <- meth_prob[idx]

#correlation analysis
cor(meth_level, meth_prob, method = "spearman")
#0.3461575

#boxplot
grade <- ceiling(meth_level*10)
df <- data.frame(score = meth_prob, level = grade)
df <- df[df$level != 0,]

pdf('Correlation.pdf')
boxplot(score~level, data=df, notch=TRUE, 
        col=(c("gold","darkgreen")),
        main="Caenorhabditis elegans", 
        xlab="Methylation Level", ylab = "Predicted Score")
dev.off()

