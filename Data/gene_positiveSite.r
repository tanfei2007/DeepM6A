
rm(list = ls)

gene_posSite<- function(input_file, output_file){
  bed_file = read.table(input_file, sep = '\t')
  posSite = bed_file[,-c(3,5)]
  posSite[,2] = posSite[,2] + 10
  colnames(posSite) = c('chr', 'coordinate', 'strand', 'IPD_ratio', 'methy_level', 'low_95', 'up_95')
  write.csv(file = output_file, posSite, row.names = F, quote = F)


# for D. melanogaster
gene_posSite()  
  
  
  
