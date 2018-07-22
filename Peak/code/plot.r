
require(ggplot2)


x = read.csv('peak_prob.csv', header=TRUE)

pl <- ggplot(x, aes(x = "", y = prob)) + xlab('') + ylab('Prob') +
  geom_boxplot()
ggsave('Peak.pdf', pl)