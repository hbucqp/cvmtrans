#!/usr/bin/python3

#!/usr/bin/python3

from Bio import SeqIO
import pandas as pd
import numpy as np
import pysam


def is_gene_within_cds(cds_coordinates, gene_feature):
    """检查基因是否在CDS范围内"""
    for current_coords in cds_coordinates:
        if current_coords[0] <= gene_feature.location.start.position \
                and current_coords[1] >= gene_feature.location.end.position:
            return True
    return False


def cds_locations(embl_file):
    """提取CDS的位置"""
    cds_coordinates = []
    for record in SeqIO.parse(embl_file, "embl"):
        for feature in record.features:
            if feature.type == "CDS":
                cds_coordinates.append(
                    (feature.location.start.position, feature.location.end.position))
    return cds_coordinates


def get_feature_id(feature):
    """获取特征ID"""
    if "locus_tag" in feature.qualifiers:
        return feature.qualifiers["locus_tag"][0]
    elif "ID" in feature.qualifiers:
        return feature.qualifiers["ID"][0]
    elif "systematic_id" in feature.qualifiers:
        return feature.qualifiers["systematic_id"][0]
    else:
        return f"{feature.id}_{feature.location.strand}_{feature.location.start.position}_{feature.location.end.position}"


def get_gene_name(feature):
    """获取基因名称"""
    if "gene" in feature.qualifiers:
        return feature.qualifiers["gene"][0]
    else:
        return get_feature_id(feature)


def get_product_value(feature):
    """获取产品信息"""
    if "product" in feature.qualifiers:
        return feature.qualifiers["product"][0]
    elif "pseudo" in feature.qualifiers:
        return "pseudogene"
    else:
        return ""


# df = pd.read_excel('instert_pos.xlsx')

# insert_sites = dict(zip(df['ins_position'], df['Count']))


# Replace 'your_embl_file.embl' with the actual path to your EMBL file
embl_file = 'H16SD2BY4.embl'

cds_coordinates = cds_locations(embl_file)
# Parse the EMBL file into SeqRecord objects
records = SeqIO.parse(embl_file, 'embl')
with open('embl_test.txt', 'w') as output:
    output.write(
        'locus_tag\tgene_name\tncrna\tstart\tend\tstrand\tread_count\tins_count\tins_index\tproduct')
    for record in records:
        # print(dir(record))
        # Access information from the SeqRecord
        print(f"Accession: {record.name}")
        print(f"Description: {record.description}")
        print(f"Sequence length: {len(record.seq)}")
        # ... other relevant information

        # Access features (e.g., CDS, gene, etc.)
        for feature in record.features:
            if feature.type == "gene" and not is_gene_within_cds(cds_coordinates, feature):
                continue
            if feature.type in ["CDS", "polypeptide"]:
                feature_id = get_feature_id(feature)
                gene_name = get_gene_name(feature)
                product_value = get_product_value(feature)
                rna_value = "1" if "ncRNA" in feature.qualifiers else "0"
                strand = feature.location.strand
                start = feature.location.start
                end = feature.location.end - 1
                feature_start = start
                feature_end = end + 1
                # print(f'{feature_start}\t{feature_end}')
                print(
                    f'{feature.type}\t{feature_id}\t{gene_name}\t{rna_value}\t{start}\t{end}\t{strand}\t{product_value}')


# read bam file
# bamfile = pysam.AlignmentFile("SC2021ZGS032_sorted.bam", "rb")

# read_counts = {}
# for read in bamfile:
#     if not read.is_unmapped and read.mapping_quality >= 30:
#         # print(dir(read))
#         if read.is_reverse:
#             strand = -1
#         else:
#             strand = 1
#         # get the aligned information, pysam already including the following start and end postion, so we skip calculate the start or end postion using cigar string and position field in bam file

#         # length = read.reference_length
#         # cigar = read.cigarstring
#         # cigar_results = cigar_parse(cigar, pos)
#         aligned_start = read.pos
#         aligned_end = read.reference_end - 1

#         # get the reference sequence name
#         ref_name = read.reference_name

#         # decision the insert site based on the strand
#         if strand == 1:
#             ins_pos = aligned_start
#         else:
#             ins_pos = aligned_end

#         read_counts.setdefault(ref_name, {}).setdefault(ins_pos, 0)
#         read_counts[ref_name][ins_pos] += 1

#         print(
#             f'{aligned_start}\t{aligned_end}\t{ins_pos}\t{ref_name}\t{strand}')
# df = pd.DataFrame.from_dict(read_counts)
# print(df)
