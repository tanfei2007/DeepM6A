#**************************************************************
# combine predict result and sequence information(Example: Ara)
# Xiurui Hou
# xh256@njit.edu
#*************************************************************
import csv

chr_name = ['NC_003070.9','NC_003071.7','NC_003074.8','NC_003075.7','NC_003076.8']
chr_nature = ['neg', 'pos']

csv_file_dir = '/cstor/xsede/users/xs-tanfei/m6a/ara/whole_genome/rslt/'
fa_file_dir = '/cstor/xsede/users/xs-ttgump/Ara/genome/'
output_dir = './comb_rslt_fa_output/'
for name in chr_name:
    for nature in chr_nature:
        output = open(output_dir + name+'.'+nature+'_output.csv', 'w+')
        csv_file = open(csv_file_dir + name + '.' + nature + '.csv', 'r')
        fa_file = open(fa_file_dir + name + '.' + nature + '.fa', 'r')
        csv_line = csv_file.readline()
        csv_line = csv_file.readline()
        fa_line = fa_file.readline()
        output.write('chr_name, neg/pos, position, ground truth, pred_prob\n')
        while csv_line and fa_line:
            parts = csv_line.split(',')
            if True:
                seq_range = fa_line.split(':')[1]
                number = seq_range.split('-')
                if len(number) == 3:
                    middle_pos = 30 - int(number[1])
                else:
                    middle_pos = 30 + int(number[0])
                if nature == 'pos':
                    data = name+',+,'+str(middle_pos)+','+parts[0]+','+parts[1]
                elif nature =='neg':
                    data = name+',-,'+str(middle_pos)+','+parts[0]+','+parts[1]
                else: 
                    print 'Error!'
                output.write(data)
                print data
            csv_line = csv_file.readline()
            fa_line = fa_file.readline()
            fa_line = fa_file.readline()
