library(ggplot2)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
prefix = args[1]

df <- read.table(paste(prefix, '.txt', sep=''), header = F)
n <- nrow(df)-1
df$position <- (-n/2):(n/2)
colnames(df) <- c('score', 'position')

ggplot(data=df, aes(x=position, y=score, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point(color="blue", size=2)+
  xlab("position")+
  scale_x_discrete(limit = seq(-n/2, n/2, 5))+
  ggtitle("C.elegans")	

file_name = paste(prefix, '.pdf', sep='')
ggsave(file_name)
