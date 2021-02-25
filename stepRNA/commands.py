#Python default modules

def left_overhang(sorted_bam, line, ref_positions):
    '''Get the length of the left overhang; 0 = no overhang, -ve = reference overhang, +ive = query overhang
    
    sorted_bam [PYSAM_AlignFile] - an open pysam.AlignmentFile
    line [PYSAM_AlignFile_Record] - a pysam.AlignmentFile individual record
    ref_positions [PYSAM_Record_RefPos] - a pysam.AlignmentFile record.get_reference_positions()
    
    Returns: The length of the left overhang and the type of overhang'''
    if ref_positions[0] == 0:
        return 0, '5primeExact'
    elif ref_positions[0] == None:
        if line.reference_start != 0:
            raise Exception
        else:
            return line.query_alignment_start, '5primeRead'
    else:
        #If reference_position[0] > 0
        return -ref_positions[0], '5primeRef'

def right_overhang(sorted_bam, line, ref_positions):
    '''Get the length of the right overhang; 0 = no overhang, -ve = reference overhang, +ive = query overhang

    sorted_bam [PYSAM_AlignFile] - an open pysam.AlignmentFile
    line [PYSAM_AlignFile_Record] - a pysam.AlignmentFile individual record
    ref_positions [PYSAM_Record_RefPos] - a pysam.AlignmentFile record.get_reference_positions()
    
    Returns: The length of the left overhang and the type of overhang'''
    ref_length = sorted_bam.lengths[sorted_bam.get_tid(line.reference_name)]
    if ref_positions[-1] == ref_length - 1:
        return 0, '3primeExact'
    elif ref_positions[-1] == None:
        if line.reference_end != ref_length:
           raise Exception
        else:
            return ref_length - line.query_alignment_length , '3primeRead'
            #return ref_length - line.query_length, 'RQ'
    else:
        #If reference_position[-1] < ref_length
        return ref_positions[-1] - (ref_length - 1), '3primeRef'


