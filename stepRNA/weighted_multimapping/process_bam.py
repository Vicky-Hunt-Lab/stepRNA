import pysam

def check_unmapped(alignment):
    return (not alignment.is_unmapped) and (('D' and 'I') not in alignment.cigarstring)


def make_mmap_dict(bamfile):
    alignmentfile = pysam.AlignmentFile(bamfile)
    multi_dict = {}
    for alignment in alignmentfile:
        #Check to see if mapped and not Insertions or Deletions
        if check_unmapped(alignment):
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

    return alignmentfile, (count, ref_names)

def check_secondary(alignment):
    return (alignment.is_secondary) and (('D' and 'I') not in alignment.cigarstring)

def get_multimappers(bamfile):
    alignmentfile = pysam.AlignmentFile(bamfile)
    ref_names = None 
    print('Before loop')
    previous_ref_name = None
    for alignment in alignmentfile:
        print('###')
        print(check_secondary(alignment))
        print(previous_ref_name)
        if check_secondary(alignment):
            ref_names = [previous_ref_name]
            print(ref_names)
            while alignment.is_secondary:
                ref_names.append(alignment.reference_name)
                alignment = next(alignmentfile)
            yield ref_names
        else:
            previous_ref_name = alignment.reference_name
            print(previous_ref_name)

def compare_alignments(names, q_name, r_name, alignment_score):
    comparison_score = names[q_name][1]
    if alignment_score > comparison_score:
        names[q_name] = [r_name, alignment_score]

    return names

def get_alignmentscore(alignment):
    for tag in alignment.tags:
        if 'AS' in tag:
            return tag[1]

def main(bamfile):
    multi_dict = make_mmap_dict(bamfile)
    print(multi_dict)
    a = input()
    if a is not 'y':
        raise Exception()
        return multi_dict
    multimap_gen = get_multimappers(bamfile)
    err_count = 0
    for ref_names in multimap_gen:
        print(ref_names)
        fractional_read_dict = {}
        for ref_name in ref_names:
            try:
                fractional_read_dict[ref_name] = multi_dict[ref_name]
            except KeyError:
                print(ref_names)
                err_count += 1
                print('#############')
                print(err_count)
                a = input()
                if a is not 'y':
                    raise Exception()
                    return multi_dict

        yield fractional_read_dict
