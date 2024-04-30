#!/usr/bin/python3

from Bio import SeqIO
import pandas as pd
import numpy as np

# Replace 'your_embl_file.embl' with the actual path to your EMBL file
embl_file = 'H16SD2BY4.embl'


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


df = pd.read_excel('instert_pos.xlsx')

insert_sites = dict(zip(df['ins_position'], df['Count']))


# Parse the EMBL file into SeqRecord objects
records = SeqIO.parse(embl_file, 'embl')


cds_coordinates = cds_locations(embl_file)


with open('insert_stat.txt', 'w') as output:
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
            if feature.type in ["CDS", "polypeptide", "gene"]:
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
                ins_index = inserts / (end - start + 1) if end != start else 0
                #         count += sum(insert_sites[j])
                #         inserts += 1 if sum(insert_sites[j]) > 0 else 0
                # ins_index = inserts / \
                #     (end - start +
                #      1) if end != start else 0
                print(
                    f'{feature_id}\t{gene_name}\t{rna_value}\t{feature_start}\t{feature_end}\t{strand}\t{count}\t{inserts}\t{ins_index}\t{product_value}')
                output.write(
                    f'{feature_id}\t{gene_name}\t{rna_value}\t{feature_start}\t{feature_end}\t{strand}\t{count}\t{inserts}\t{ins_index}\t{product_value}\n')
