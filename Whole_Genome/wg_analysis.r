
#**********
# libs
#**********
rm(list = ls())
library(readr)
source('Evaluation.r')

#****************
# read data
#****************
file_names <- list.files('./data/', pattern = '*.csv')
j = 1
for(file_name in file_names){
  x = read_csv(paste('./data/', file_name, sep = ''))
  
  if(j == 1){
    rslt = x
  }else{
    rslt = rbind(rslt, x)
  }
  
  j = j +1 
}

save(file = 'rslt.RData', rslt, compression_level = 9)

#*******************
# analyze performance
#******************

cutoffs = c(seq(0.5, 0.9, 0.1), seq(0.91, 0.99, 0.01))
j = 1
for(i in cutoffs){
  perf <- eval(post = rslt$pred_prob, ground.truth = rslt$ground_truth, thres = i)
  
  if(j == 1)
      perfs = perf
  else{
      perfs = rbind(perfs, perf)
  }
  
  print(j)
  j = j + 1
}

write.csv(file = 'perfs.csv', perfs, row.names = F, quote = T)


#conduct analysis on cutoff 
load('rslt.RData')
proportion = c(seq(0.1, 0.02, -0.01),
               seq(0.01, 0.001, -0.001))
N = nrow(rslt)
rslt = sort(rslt[,2], decreasing = T, method = 'quick')
ns = c()
for (p in proportion){
  ns = c(ns, ceiling(N*p))
}
cutoff_prob = rslt[ns]
cutoff_prop_prob = data.frame(prop = proportion, prob = cutoff_prob)
write.csv(file='Celegans.csv', cutoff_prop_prob, row.names = F)


# prop      prob
# 1  0.100 0.5359123
# 2  0.090 0.5876898
# 3  0.080 0.6408077
# 4  0.070 0.6935574
# 5  0.060 0.7449306
# 6  0.050 0.7938480
# 7  0.040 0.8399874
# 8  0.030 0.8832679
# 9  0.020 0.9236561
# 10 0.010 0.9614733
# 11 0.009 0.9651591
# 12 0.008 0.9688153
# 13 0.007 0.9724771
# 14 0.006 0.9761189
# 15 0.005 0.9797641
# 16 0.004 0.9834235
# 17 0.003 0.9871500
# 18 0.002 0.9909271
# 19 0.001 0.9949517


