#***************************************
# convert .csv file to .html file
# Xiurui Hou
# xh256@njit.edu
#***************************************
require(DT)
Args <- commandArgs()
file <- Args[6]
html_file <- Args[7]
setwd('/cstor/xsede/users/xs-xh256/6ma2/Ara/WG/merge/csv_to_html/csv_parts/')
csv <- read.csv(file)
y<-datatable(csv)
setwd('/cstor/xsede/users/xs-xh256/6ma2/Ara/WG/merge/csv_to_html/html_parts/')
DT::saveWidget(y, html_file)
