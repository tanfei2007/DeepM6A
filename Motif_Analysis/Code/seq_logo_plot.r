require(seqLogo)


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
prefix = args[1]

m <- read.table(paste(prefix, '.txt', sep=''), header = F)
pwm <- makePWM(m)
postscript(file = paste(prefix, '.eps', sep=''), width=20, height = 3)
seqLogo(pwm, ic.scale=FALSE, xaxis=T, yaxis=T, xfontsize=15, yfontsize=15)
dev.off()
