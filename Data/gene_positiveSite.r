
rm(list = ls)

gene_posSite<- function(input_file, output_file){
  bed_file = read.table(input_file, sep = '\t')
  posSite = bed_file[,-c(3,5)]
