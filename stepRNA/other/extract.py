


import pysam
from collections import defaultdict

class Error(Exception):
    """Base class for exceptions"""
    pass

class LeftOverhangError(Error):
    """Exception raised for 3' Overhang not being at the end of the reference"""
    def __init__(self, expression, message):
        self.epression = expression
        self.message = message

class RightOverhangError(Error):
    """Exception raised for 3' Overhang not being at the end of the reference"""
    def __init__(self, expression, message):
        self.epression = expression
        self.message = message

    

def left_overhang(sam_entry, ref_positions):
    '''Get the length of the left overhang; 0 = no overhang, -ve = reference overhang, +ive = query overhang'''
    ref_positions = sam_entry.get_reference_positions(full_length=True) 
    if ref_positions[0] == 0:
        return 0
    elif ref_positions[0] == None:
        if line.reference_start != 0:
            raise Error(line.reference_start, 'Overhang does not start at the end of the read')
        else:
            return sam_entry.query_alignment_start
    else:
        #If reference_position[0] > 0
        return -ref_positions[0]



def right_overhang(samfile, sam_entry, ref_positions):
    '''Get the length of the right overhang; 0 = no overhang, -ve = reference overhang, +ive = query overhang'''
    ref_length = samfile.lengths[samfile.get_tid(line.reference_name)]
    if ref_positions[-1] == ref_length - 1:
        return 0
    elif ref_positions[-1] == None:
        if sam_entry.reference_end != ref_length:
            raise RightOverhangError(line.reference_end, 'Overhang does not start at the end of the read')
        else:
            return ref_length - sam_entry.query_alignment_end
    else:
        #If reference_position[-1] < ref_length
        return ref_positions[-1] - (ref_length - 1)

# Extracting the overhangs...
samfile = pysam.AlignmentFile('Cel_M_passenger.bam', 'rb')
right_oh_dic = defaultdict(lambda:0)
for line in samfile:
    if line.cigarstring != None:
        if ('D' or 'I') not in line.cigarstring:
            ref_pos = line.get_reference_positions(full_length = True)
            try:
                right = str(right_overhang(samfile, line, ref_pos))
                right_oh_dic[right] += 1
            except RightOverhangError:
                continue

# Having a look at the overhang frequencies...
right_oh_lst = list(right_oh_dic)
right_oh_lst.sort()
temp = right_oh_lst[:8]
temp.reverse()
right_oh_lst = temp + right_oh_lst[8:]
right_oh_lst_values = []
for item in right_oh_lst:
    right_oh_lst_values.append(right_oh_dic.get(item))

for item in range(len(right_oh_lst)):
    print(right_oh_lst[item],'has', right_oh_lst_values[item])


total = 0
for key in right_oh_lst:
    total += right_oh_dic[key]


density = []
for item in right_oh_lst:
    density.append(100 * right_oh_dic.get(item)/total)

# With 0 in the range!
for item in range(len(density)):
    length = int(density[item]) * '.'
    print('{}\t| {}'.format(right_oh_lst[item], length))

# Without 0 in the range!
# NOTE the 0 ones seem to generally be exact matches! (16108 26M in the sam file)
for item in range(len(density)):
    length = int(12 * density[item]) * '.'
    if right_oh_lst[item] == '0':
        pass
    else:
        print('{}\t| {}'.format(right_oh_lst[item], length))



# Without including 26M cigar strings!
samfile = pysam.AlignmentFile('Cel_M_passenger.bam', 'rb')
alt_oh_dic = defaultdict(lambda:0)
for line in samfile:
    if line.cigarstring != None:
        if line.cigarstring != '26M':
            if ('D' or 'I') not in line.cigarstring:
                ref_pos = line.get_reference_positions(full_length = True)
                try:
                    alt = str(right_overhang(samfile, line, ref_pos))
                    alt_oh_dic[alt] += 1
                except RightOverhangError:
                    continue

# Having a look at the overhang frequencies...
alt_oh_lst = list(alt_oh_dic)
alt_oh_lst.sort()
temp = alt_oh_lst[:8]
temp.reverse()
alt_oh_lst = temp + alt_oh_lst[8:]
alt_oh_lst_values = []
for item in alt_oh_lst:
    alt_oh_lst_values.append(alt_oh_dic.get(item))

for item in range(len(alt_oh_lst)):
    print(alt_oh_lst[item],'has', alt_oh_lst_values[item])


total = 0
for key in alt_oh_lst:
    total += alt_oh_dic[key]


density = []
for item in alt_oh_lst:
    density.append(100 * alt_oh_dic.get(item)/total)

# With 0 in the range!
for item in range(len(density)):
    length = int(density[item]) * '.'
    print('{}\t| {}'.format(alt_oh_lst[item], length))
