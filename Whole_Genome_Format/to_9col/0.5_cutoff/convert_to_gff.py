#**************************************************************
# comvert cutoff .csv file to .gff file
# Xiurui Hou
# xh256@njit.edu
#***************************************************************
chr_name = ['NC_003070.9','NC_003071.7','NC_003074.8','NC_003075.7','NC_003076.8']
chr_nature = ['neg', 'pos']

cut_off_dir = './cut_off/'
cut_off_output = './9col_cut_off/'
truth1_dir = './truth1/'
truth1_output = './9col_truth1/'
for name in chr_name:
    for nature in chr_nature:
        input_file = open(cut_off_dir+name+'.'+nature+'_output.csv', 'r')
        output_file = open(cut_off_output+name+'.'+nature+'_output.gff', 'w+')
        #input_file = open(truth1_dir+name+'.'+nature+'_output.csv', 'r')
        #output_file = open(truth1_output+name+'.'+nature+'_output.gff', 'w+')
        line = input_file.readline()
        line = input_file.readline()
        while line:
            parts = line.split(',')
            data = name + '\tprediction\t.\t' + parts[2] + '\t' + parts[2] + '\t' +parts[4].replace('\n','') + '\t'+parts[1]+'\t'+parts[3]+'\t.\n'
            output_file.write(data)
            print line
            line = input_file.readline()

