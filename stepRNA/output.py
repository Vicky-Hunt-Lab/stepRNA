#Python modules
import csv

def write_to_bam(line, left_type, right_type, prefix):
    '''Take a pysam.AlignmentFile record and write it to a file
    
    line [PYSAM_AlignFile_Record] - a pysam.AlignmentFile individual record
    left_type [STR] - type of left overhang [LE, LQ, LR]
    right_type [STR] - type of right overhang [RE, RQ, RR]
    prefix [STR] - filepath prefix

    Returns: Nothing but creates appends to a BAM file
    '''
    with open(prefix + '_' + left_type + '_' + right_type + '.bam', 'a+') as bam_out:
        bam_out.write(line.to_string() + '\n')

def make_csv(dics, csv_name, headers, show=True):
    '''Make a csv continaing the overhang information

    dics [LST] - a list containingng the [right_dic, left_dic] order
    csv_name [STR] - filename for the csv
    headers [LST] - a list of the header information for the csv
    show [True/False] - whether to show on command line
    
    Returns: Nothing, outputs a csv file to csv_name filepath'''
    keys = set()
    for dic in dics:
        for key in dic:
            keys.add(key)
    keys = list(keys)
    keys.sort()
    with open(csv_name, 'w') as csv_out:
        writer = csv.writer(csv_out, delimiter=',')
        writer.writerow(headers)
        if show:
            print('\t'.join(headers))
        for key in keys:
            if dics[0][key] == None:
                dics[0][key] = 0
            if dics[1][key] == None:
                dics[1][key] = 0
            writer.writerow([key, dics[1].get(key), dics[0].get(key)])
            if show:
                print('{}\t{}\t{}'.format(key, dics[1].get(key), dics[0].get(key)))

def make_type_csv(dic, csv_name, headers, show=True, sort=False):
    '''Make a csv from a dictionary with the keys forming the row infomration anthe values forming the count information 

    dic [DICT] - A dictionary containing key: count pairs
    csv_name [STR] - filename for the csv
    headers [LST] - a list of the header information for the csv
    show [True/False] - whether to show on command line
    
    Returns: Nothing, outputs a csv file to csv_name filepath'''
    if sort:
        keys = set()
        for key in dic:
            keys.add(key)
        keys = list(keys)
        keys.sort()
    else:
        keys=dic.keys()
    with open(csv_name, 'w') as csv_out:
        writer = csv.writer(csv_out, delimiter = ',')
        writer.writerow(headers)
        if show:
            print('\t'.join(headers))
        for key in keys:
            writer.writerow([key, dic[key]])
            if show:
                print('{}\t{}'.format(key, dic[key]))

def print_hist(density_dic, keys):
    '''Print a basic histogram to the terminal from a dictionary of density values
    
    density_dic [DIC] - a dictionary containing key: density pairs
    keys [LST] - a list of the keys for the dictionary'''
    for item in range(len(density_dic)):
        try:
            length = int(density_dic[item]) * '.'
            print('{}\t| {}'.format(keys[item], length))
        except ValueError:
            print('{}\t{}'.format(keys[item], ''))

