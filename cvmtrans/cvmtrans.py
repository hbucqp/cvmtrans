#!/usr/bin/python3
import os
import re
import sys
import pysam
import pandas as pd
import numpy as np
import tempfile
import argparse
from Bio import SeqIO


def args_parse():
    "Parse the input argument, use '-h' for help."
    parser = argparse.ArgumentParser(
        usage='cvmtrans -bam_file <sorted_bam_files> -embl_file <reference embl file> -output <tab delimited text file> \n\nAuthor: Qingpo Cui(SZQ Lab, China Agricultural University)\n')
    # group = parser.add_mutually_exclusive_group(required=False)
    parser.add_argument(
        "-bam_file", help="<input_fastq>: output filename of the fastq file")
    parser.add_argument(
        "-embl_file", help="<input_file>: the tag sequence")
    parser.add_argument(
        "-output", help="<output_fastq>: output filename of the fastq file")

    # parser.add_argument('-v', '--version', action='version',
    #                     version='Version: ' + get_version("__init__.py"), help='<display version>')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def bam2insdict(bamfile):
    """
    Genreate the insertion dictionary from bam files
    Paramters:
    ----------
    bamfile: str
        The sorted bam file of transposon sequencing data
    Returen:
    ----------
    A dictionary of the insertion site as key, and the number of reads of the insertion site
    {"seq_name":{1:10, 3: 2, 4: 300}}
    """
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
            # print(
            #     f'{aligned_start}\t{aligned_end}\t{ins_pos}\t{ref_name}\t{strand}')
    return read_counts


def is_gene_within_cds(cds_coordinates, gene_feature):
    """
    检查基因是否在CDS范围内
    """
    for current_coords in cds_coordinates:
        if current_coords[0] <= gene_feature.location.start.position \
                and current_coords[1] >= gene_feature.location.end.position:
            return True
    return False


def cds_locations(embl_file):
    """
    Extract the CDS coordinates.

    Parameters
    -----------
    embl_file: str
        The embl file name

    Returns
    -----------
    A list of tuple that contains the CDS start and end position
    """
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


def get_insertion_stat(bam_file, embl_file, output):
    """
    Get insertion information from sorted bam files

    Parameters:
    -----------
    bam_file: str
        The path of bam format file
    embl_file: str
        The path of the embl reference file
    output: str
        The output filename
    """

    # parse bamfile
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    insert_read_counts = bam2insdict(bamfile)

    # parse the embl file
    embl_file = embl_file
    cds_coordinates = cds_locations(embl_file)

    # process the insertion sites
    records = SeqIO.parse(embl_file, 'embl')
    for record in records:
        seq_id = record.id
        if seq_id == '':
            sys.exit(
                "Could not find sequence name from your given embl file. Exiting ...")

        # print(dir(record))
        # Access information from the SeqRecord
        # print(f"Accession: {record.id}")
        # print(f"Description: {record.description}")
        # print(f"Sequence length: {len(record.seq)}")
        # ... other relevant information

        # check if embl seq id in the insertion site dictionary

        try:
            insert_sites = insert_read_counts[record.id]
        except:
            sys.exit("Could not found the corresponding seq name %s in insertion dictionary, please check your embl file and make sure your mapping reference file has the same header with your embl file." % record.id)
        with open(output, 'w') as out:
            out.write(
                'locus_tag\tgene_name\tncrna\tstart\tend\tstrand\tread_count\tins_count\tins_index\tproduct\n')
            # Access features (e.g., CDS, gene, etc.)
            for feature in record.features:
                if feature.type == "gene" and not is_gene_within_cds(cds_coordinates, feature):
                    continue
                # if feature.type in ["CDS", "polypeptide", "gene"]: remove gene in order to delete duplicates annotation
                if feature.type in ["CDS", "polypeptide"]:
                    feature_id = get_feature_id(feature)
                    gene_name = get_gene_name(feature)
                    product_value = get_product_value(feature)
                    rna_value = "1" if "ncRNA" in feature.qualifiers else "0"
                    strand = feature.location.strand

                    start = feature.location.start
                    end = feature.location.end - 1
                    feature_start = start + 1
                    feature_end = end + 1
                    # print(
                    # f'{feature_id}\t{gene_name}\t{rna_value}\t{start}\t{end}\t{strand}\t{product_value}')

                    # calculate the insert_count and the total reads in the corresponding insertion gene
                    count = 0
                    inserts = 0

                    for j in np.arange(start, end + 1):
                        if j in insert_sites.keys():
                            # print(j)
                            count += insert_sites[j]
                            inserts += 1
                    ins_index = inserts / \
                        (end - start + 1) if end != start else 0
                    #         count += sum(insert_sites[j])
                    #         inserts += 1 if sum(insert_sites[j]) > 0 else 0
                    # ins_index = inserts / \
                    #     (end - start +
                    #      1) if end != start else 0
                    print(
                        f'{feature_id}\t{gene_name}\t{rna_value}\t{feature_start}\t{feature_end}\t{strand}\t{count}\t{inserts}\t{ins_index}\t{product_value}')
                    out.write(
                        f'{feature_id}\t{gene_name}\t{rna_value}\t{feature_start}\t{feature_end}\t{strand}\t{count}\t{inserts}\t{ins_index}\t{product_value}\n')


def main():
    args = args_parse()
    bamfile = os.path.abspath(args.bam_file)
    emblfile = os.path.abspath(args.embl_file)
    current_path = os.path.abspath(os.path.dirname(__file__))
    output = args.output

    if re.search(r'[\/]', output):
        sys.exit('Please check your output file name, do not give a PATH string')
    else:
        output = os.path.join(current_path, args.output)

    get_insertion_stat(bamfile, emblfile, output)


if __name__ == '__main__':
    main()
