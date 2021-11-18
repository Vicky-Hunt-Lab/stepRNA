import pysam

def make_mmap_dict(bamfile):
    alignmentfile = pysam.AlignmentFile(bamfile)
    multi_dict = {}
    for alignment in alignmentfile:
        #Check to see if mapped and not Insertions or Deletions
        if not alignment.is_unmapped and ('D' or 'I') not in alignment.cigarstring:
            #Search for duplicates (using best score)
            #Add counts to dictionary for each reference
            alignmentfile, ref_info = get_mapping_info(alignment, alignmentfile)
            mmap = ref_info[0]
            for ref_name in ref_info[1]:
                try:
                    multi_dict[ref_name][mmap] += 1
                except KeyError:
                    try:
                        multi_dict[ref_name][mmap] = 1 
                    except KeyError:
                        multi_dict[ref_name] = {mmap : 1}
        else:
            continue

    return multi_dict

def get_mapping_info(alignment, alignmentfile):
    '''alignment is a pysam.AlignedSegment object'''
    query_name = alignment.query_name
    count = 0
    names = {}
    #Get all the queries with that name
    while alignment.query_name == query_name:
        q_name = alignment.query_name
        r_name = alignment.reference_name
        alignment_score = get_alignmentscore(alignment)
        if r_name not in names:
            count += 1
            names[r_name] = [q_name, alignment_score]
            alignment = next(alignmentfile)
        else:
            names = compare_alignments(names, r_name, q_name, alignment_score)
            alignment = next(alignmentfile)

    ref_names = []
    _ = [ref_names.append(ref_name) for ref_name in names.keys()]
    if len(ref_names) > 1:
        print(ref_names)

    return alignmentfile, (count, ref_names)

def compare_alignments(names, q_name, r_name, alignment_score):
    comparison_score = names[q_name][1]
    if alignment_score > comparison_score:
        names[q_name] = [r_name, alignment_score]

    return names

def get_alignmentscore(alignment):
    for tag in alignment.tags:
        if 'AS' in tag:
            return tag[1]
