#Python modules
import csv
import pysam
from collections import defaultdict
import numpy as np
import sys

def refs_counts(bamfile, unique=False):
    '''Take a BAM file and count the number of unique references within it'''
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    if unique:
        return len(set(samfile.header.references))
    else:
        return len(samfile.header.references)

def read_len_counts(bamfile):
    '''Take a BAM file and count the number of reads with each length.
    
    Returns a dictionary of {length : count}'''
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    dic = defaultdict(lambda:0) 
    for record in samfile:
        dic[record.query_length] += 1
    return dic

def ref_read_counts(bamfile, ref_list):
    '''For a list of reference names count their occurences in a BAM file'''
    for ref in ref_list:
        dic[ref] = 0
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    for record in samfile:
        dic[record.reference_name] += 1
    return dic



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

#Needs changing for logodds
def make_csv(dics, csv_name, headers, logger, show=True):
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
            logger.write('\t'.join(headers))
        for key in keys:
            if dics[0].get(key) == None:
                threeprime = [0,'NA', 'NA']
            else:
                threeprime = dics[0].get(key)
            if dics[1].get(key) == None:
                fiveprime = [0, 'NA', 'NA']
            else:
                fiveprime = dics[1].get(key)
#            try:
#                fiveprime = dics[1].get(key)
#            except KeyError:
#                print('KeyError')
#                fiveprime = [0, 'NA', 'NA']
#            try:
#                threeprime = dics[0].get(key)
#            except KeyError:
#                print('Key error')
#                threeprime = [0, 'NA', 'NA']
            writer.writerow([key, fiveprime[0], threeprime[0], fiveprime[1], threeprime[1], fiveprime[2], threeprime[2]])
            if show:
                logger.write('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(key, fiveprime[0], threeprime[0], fiveprime[1], threeprime[1], fiveprime[2], threeprime[2]))

def make_type_csv(dic, csv_name, headers, logger, show=True, sort=False):
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
            logger.write('\t'.join(headers))
        for key in keys:
            writer.writerow([key, dic[key]])
            if show:
                logger.write('{}\t{}'.format(key, dic[key]))

def print_hist(density_dic, keys, logger):
    '''Print a basic histogram to the terminal from a dictionary of density values
    
    density_dic [DIC] - a dictionary containing key: density pairs
    keys [LST] - a list of the keys for the dictionary'''
    for item in range(len(density_dic)):
        try:
            length = int(density_dic[item]) * '.'
            logger.write('{}\t| {}'.format(keys[item], length))
        except ValueError:
            logger.write('{}\t{}'.format(keys[item], ''))

class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, 'w')

    def write(self, message):
        self.terminal.write(message + '\n')
        self.log.write(message + '\n')

    def close(self):
        self.log.close()

def oddsratio(counts):
    '''Input is a oh_length : count dictionary; output is a oh_length : [count, logodds, zscore] dictionary'''
    total = 0
    est_std = 0
    for key in counts:
        total += counts[key]
        est_std += 1 / counts[key]
    est_std = np.sqrt(est_std)
    xbar = total / len(counts)

    outdictionary = {}
    for key in counts:
        phit = counts[key] / total
        qhit = 1 - phit
        poff = xbar / total
        qoff = 1 - poff

        count = counts[key]
        odds = (phit / qhit) / (poff / qoff)
        logodds = np.log(odds)
        zscore = logodds / est_std

        outdictionary[key] = [count, logodds, zscore]

    return outdictionary

