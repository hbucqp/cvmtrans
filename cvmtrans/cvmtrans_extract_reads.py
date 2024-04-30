#!/usr/bin/python3

import dnaio
import os
import re
import sys
from pathlib import Path
import argparse


def args_parse():
    "Parse the input argument, use '-h' for help."
    parser = argparse.ArgumentParser(
        usage='cvmtrans_extract_reads -i <input_fastq> -tag <sequence should inclued in reads> -o <output_fastq> \n\nAuthor: Qingpo Cui(SZQ Lab, China Agricultural University)\n')
    # group = parser.add_mutually_exclusive_group(required=False)
    parser.add_argument(
        "-i", help="<input_fastq>: output filename of the fastq file")
    parser.add_argument(
        "-tag", help="<input_file>: the tag sequence")
    parser.add_argument(
        "-o", help="<output_fastq>: output filename of the fastq file")

    # parser.add_argument('-v', '--version', action='version',
    #                     version='Version: ' + get_version("__init__.py"), help='<display version>')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def get_version(rel_path: str) -> str:
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")


def read(rel_path: str) -> str:
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path)) as fp:
        return fp.read()


def split_extensions(path: str) -> tuple[str, str]:
    """
    split the extension of file name and extension as tuple
    T3.fastq.gz -> (T3, fastq.gz)
    """
    ext = "".join(Path(path).suffixes)
    return str(path).rstrip(ext), ext


def filter_reads(input_file, output_file, tag):
    """
    Extract reads contains the given tag in input file and write to output files

    Parameters
    ----------
    input : str
        The input filename
    output : str
        The output filename
    tag : string
        The seuqence of tags

    Returns
    -------
    A fastq file that make sure all reads containing the tag sequence

    """
    tag = str(tag).upper()

    with dnaio.open(input_file) as reader, dnaio.open(output_file, mode='w') as writer:
        total_reads = 0
        saved_reads = 0
        for read in reader:
            total_reads += 1
            if re.findall(tag, read.sequence):
                saved_reads += 1
                writer.write(read)
            else:
                next
    # print(f'Total reads: {total_reads}')
    # print(f'Saved reads: {saved_reads}')
    return total_reads, saved_reads


def main():
    args = args_parse()
    input_fq = os.path.abspath(args.i)
    current_path = os.path.abspath(os.path.dirname(__file__))
    output_fq = os.path.join(current_path, args.o)
    tag = args.tag
    # print(tag)
    # print(input_fq)
    # print(output_fq)
    # print(tag)
    print(f'Start processing {args.i} ...')
    total_reads, saved_reads = filter_reads(
        input_file=input_fq, output_file=output_fq, tag=tag)
    # print(f'Total reads: {total_reads}')
    # print(f'Saved reads: {saved_reads}')

    fraction = saved_reads * 100 / total_reads
    percent = '%.2f' % fraction
    print(
        f'A total number of {total_reads} reads in {args.i} and {saved_reads} reads contain {tag.upper()}.\n{percent}% were saved to {args.o}')


if __name__ == '__main__':
    main()
