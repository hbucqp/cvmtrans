#!/usr/bin/python3

import os
import re
import pysam
import pandas as pd
import numpy as np
import tempfile


bamfile = pysam.AlignmentFile("1-18W_cuttags_2_sorted.bam", "rb")


# def cigar_parse(cigar, coordinate):
#     cigar_results = {'start': 0, 'end': 0}
#     cigar_parts = re.findall(r'(\d+[MIDNSHP=X])', cigar)
#     for cigar_item in cigar_parts:
#         match = re.match(r'(\d+)([MIDNSHP=X])', cigar_item)
#         if match:
#             number, action = match.groups()
#             if action in ['M', 'X', '=']:
#                 if cigar_results['start'] == 0:
#                     cigar_results['start'] = coordinate
#                 coordinate += int(number)
#                 # print(current_coordinate)
#                 cigar_results['end'] = coordinate - \
#                     1 if cigar_results['end'] < coordinate else cigar_results['end']
#             elif action in ['S', 'D', 'N']:
#                 if cigar_results['start'] == 0:
#                     coordinate -= int(number)
#                 cigar_results['start'] = coordinate if cigar_results['start'] == 0 else cigar_results['start']
#                 # print(results['start'])
#                 # print(current_coordinate)
#                 coordinate += int(number)
#                 cigar_results['end'] = coordinate - \
#                     1 if cigar_results['end'] < coordinate else cigar_results['end']
#     return cigar_results


with open('bam_parse_result_2.txt', 'w') as output:
    read_counts = {}
    for read in bamfile:
        if not read.is_unmapped and read.mapping_quality >= 30:
            # print(dir(read))
            if read.is_reverse:
                strand = -1
            else:
                strand = 1
            # get the aligned information, pysam already including the following start and end postion, so we skip calculate the start or end postion using cigar string and position field in bam file

            # length = read.reference_length
            # cigar = read.cigarstring
            # cigar_results = cigar_parse(cigar, pos)
            aligned_start = read.pos
            aligned_end = read.reference_end - 1

            # get the reference sequence name
            ref_name = read.reference_name

            # decision the insert site based on the strand
            if strand == 1:
                ins_pos = aligned_start
            else:
                ins_pos = aligned_end

            read_counts.setdefault(ref_name, {}).setdefault(ins_pos, 0)
            read_counts[ref_name][ins_pos] += 1

            print(
                f'{aligned_start}\t{aligned_end}\t{ins_pos}\t{ref_name}\t{strand}')
            output.write(
                f'{aligned_start}\t{aligned_end}\t{ins_pos}\t{ref_name}\t{strand}\n')


fd, path = tempfile.mkstemp()
try:
    with os.fdopen(fd, 'w') as tmp:
        # do stuff with temp file
        tmp.write('stuff')
finally:
    os.remove(path)

print(read_counts)
df = pd.DataFrame.from_dict(read_counts)
print(df)
# df.reset_index(inplace=True)
# df.columns = ['ins_position', 'Count']
# df.to_excel('forward_insert_stat.xlsx', index=False)
