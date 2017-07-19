#**************************************************************
# extract result whose probability is larger than 0.5 and whose 
# truth is 1 but probability is less than 0.5(Example: Ara)
# Xiurui Hou
# xh256@njit.edu
#***************************************************************
chr_name = ['NC_003070.9','NC_003071.7','NC_003074.8','NC_003075.7','NC_003076.8']
chr_nature = ['neg', 'pos']

count1 = 0
count2 = 0
count3 = 0
input_dir = './comb_rslt_fa_output/'
for name in chr_name:
    for nature in chr_nature:
        input_file = open(input_dir + name+'.'+nature+'_output.csv', 'r')
        output_file = open('./cut_off/'+name+'.'+nature+'_output.csv', 'w+')
        output_file.write('chr_name, neg/pos, position, ground truth, pred_prob\n')
        output_file2 = open('./truth1/'+name+'.'+nature+'_output.csv', 'w+')
        output_file2.write('chr_name, neg/pos, position, ground truth, pred_prob\n')
        line = input_file.readline()
        line = input_file.readline()
        while line:
            count1 += 1
            parts = line.split(',')
            if float(parts[4]) > 0.5:
                count2 += 1
               # print line
                output_file.write(line)
            if float(parts[4]) < 0.5 and int(parts[3]) == 1:
                count3 += 1
                output_file2.write(line)
            line = input_file.readline()
print count1
print count2
print count3
