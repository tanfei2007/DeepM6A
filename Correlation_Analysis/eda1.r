rm(list = ls())

#load data
meth_level <- read.csv("meth_level.csv", header = T, stringsAsFactors = F)
meth_level <- as.numeric(meth_level$X0)
idx <- which(!is.na(meth_level))
meth_level <- meth_level[idx]
  
meth_prob <- read.csv("meth_prob.csv", header = T, stringsAsFactors = F)
meth_prob <- meth_prob$pred_prob
meth_prob <- meth_prob[idx]

#histgram
grade <- ceiling(meth_level*10)
df <- data.frame(score = meth_prob, level = grade)


plot_level_score <- function(df, cutoff = 0.5){
  
  print(cutoff)
  idx <- which(df$score > cutoff)
  level_score <- df[idx,'level']
  
  x <- table(level_score)
  y <- table(df$level)

  len = length(y) - length(x)
  x <- c(rep(0, len), x)
  
  ratio <- x / y 
  
  pdf(paste(cutoff,'.pdf', sep = ''))
  barplot(ratio, main = paste('cutoff=', cutoff, sep=''), ylim = c(0,1))
  dev.off()

  print(ratio)
  return(ratio)
}



cutoffs <- c(seq(0.5,0.9,0.1),
             seq(0.91, 0.98, 0.01))

rslts = matrix(0, length(cutoffs), 11)
i = 1
for (value in cutoffs){
  df <- df[df$level!=0,]
  rslt= plot_level_score(df, cutoff = value)
  common_idx = intersect(1:10, names(rslt))
  common_idx = as.numeric(common_idx) + 1
  rslts[i, common_idx] = rslt
  rslts[i, 1] = value
  i = i + 1
}

rslts <- as.data.frame(rslts)
colnames(rslts) <- c('cutoff_probability', paste('level', 1:10, sep = ''))
write.csv(file='celegans.csv', rslts, row.names = F, quote = T)

