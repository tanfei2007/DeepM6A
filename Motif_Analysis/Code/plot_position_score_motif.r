
source('~/PositionScore/PositionScore.R')
source('~/PositionScore/pwm.R')
library(grid)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
prefix = args[1]

m <- read.table(paste(prefix, '.txt', sep=''), header = F)
postscript(file = paste(prefix, '.eps', sep=''), width=20, height = 3)
positionScore(m, ic.scale=FALSE, xaxis=T, yaxis=T, xfontsize=15, yfontsize=15, title = 'Celegans')
dev.off()
