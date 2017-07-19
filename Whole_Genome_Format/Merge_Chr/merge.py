#*****************************************************
# combine fraction information and predict probability
# Xiurui Hou
# xh256@njit.edu
#*****************************************************

import os

gff_file = open('ara_m6A.gff')
gff_line = gff_file.readline()
count1 = 0
count2 = 0
count3 = 0
gff_dict = {}
seq_chr = {}
seq_chr['chr1'] = 'NC_003070.9'
seq_chr['chr2'] = 'NC_003071.7'
seq_chr['chr3'] = 'NC_003074.8'
seq_chr['chr4'] = 'NC_003075.7'
seq_chr['chr5'] = 'NC_003076.8'
while gff_line:
	count3 += 1
	parts = gff_line.split()
	chr_name = seq_chr[parts[0]]
	chr_site = parts[3]
	chr_strand = parts[6]
	if ';frac=' in gff_line:
		count1 += 1
		chr_frac = parts[-1].split(';')[3].split('=')[1]
	else:
		count2 += 1
		chr_frac = 'na'
	if chr_name not in gff_dict.keys():
		gff_dict[chr_name] = {}
	if chr_strand not in gff_dict[chr_name].keys():
		gff_dict[chr_name][chr_strand] = {}
	if chr_site not in gff_dict[chr_name][chr_strand].keys():
		gff_dict[chr_name][chr_strand][chr_site] = chr_frac
	gff_line = gff_file.readline()

# print 'gff file has ' + str(count1)
# print 'gff file has no ' + str(count2)
# print 'gff total ' + str(count3)

file_dir = '/cstor/xsede/users/xs-xh256/6ma2/Ara/WG/9col/comb_rslt_fa_output/'
output_file = open('ara_merge.csv', 'w+')
output_file.write('Chromsome,Strand,Site,Label,Fraction,Probability\n')
count1 = 0
count2 = 0
count3 = 0
file_list = os.listdir(file_dir)
for i in range(0, len(file_list)):
	if True:
		csv_file = open(file_dir + file_list[i])
		csv_line = csv_file.readline()
		csv_line = csv_file.readline()
		while csv_line:
			parts = csv_line.split(',')
			chr_name = parts[0]
			chr_strand = parts[1]
			chr_site = parts[2]
			truth = parts[3]
			prob = parts[4]
			count3 += 1
			if chr_site in gff_dict[chr_name][chr_strand].keys():
				count1 += 1
				frac = gff_dict[chr_name][chr_strand][chr_site]
			else:
				count2 += 1
				frac = 'na'
			data = chr_name+','+chr_strand+','+chr_site+','+truth+','+frac+','+prob
			output_file.write(data)
			csv_line = csv_file.readline()

# print 'rslt file has ' + str(count1)
# print 'rslt file has no ' + str(count2)
# print 'rslt total ' + str(count3)



