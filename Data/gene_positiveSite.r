
rm(list = ls)

gene_posSite<- function(input_file, output_file){
  bed_file = read.table(input_file, sep = '\t')
  posSite = bed_file[,-c(3,5)]
  posSite[,2] = posSite[,2] + 10
  colnames(posSite) = c('chr', 'coordinate', 'strand', 'IPD_ratio', 'methy_level', 'low_95', 'up_95')
  write.csv(file = output_file, posSite, row.names = F, quote = F)


# for D. melanogaster
gene_posSite(input_file = '825_m6A_25x_novar_dev03.bed', output_file = 'm6A_dm_pos.csv')  

# for A. thaliana
gene_posSite(input_file = 'm6A_TAIR10.bed', output_file = 'm6A_ara_pos.csv')
  
# for A. thaliana
gene_posSite(input_file = 'm6A_ecoli.bed', output_file = 'm6A_ecoli_pos.csv')  
  
