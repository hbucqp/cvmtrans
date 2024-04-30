#!/usr/bin/python3

# import dnaio
import os
import re
import sys
import subprocess
from pathlib import Path
import argparse


def args_parse():
    "Parse the input arguments, use '-h' for help."
    parser = argparse.ArgumentParser(
        usage="cvmtrans_cut_tags -i <input_fastq> -tag_5 <5' end tag> -tag_3 <3' end tag> -o <output_fastq> \n\nAuthor: Qingpo Cui(SZQ Lab, China Agricultural University)\n")
    # group = parser.add_mutually_exclusive_group(required=False)
    parser.add_argument(
        "-i", help="<input_fastq>: output filename of the fastq file")
    parser.add_argument(
        "-tag_5", help="<5‘ end tags or barcode>: The tag sequence at the 5’ end will be removed")
    parser.add_argument(
        "-tag_3", help="<3‘ end tags or barcode>: The tag sequence at the 3’ end will be removed")
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


def run_shell(shell_command):
    """
    Run the shell command
    """
    try:
        result = subprocess.run(shell_command, shell=True,
                                check=True, capture_output=True, text=True)
        print("Command executed successfully. Output:")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error executing the command. Exit code: {e.returncode}")
        print(f"Error message: {e.stderr}")


def check_tags(tag_5, tag_3):
    if tag_5 is None and tag_3 is None:
        # return False
        print("The tag sequence is needed!")
        sys.exit(1)
    else:
        return True


def main():
    args = args_parse()
    input_fq = os.path.abspath(args.i)
    current_path = os.path.abspath(os.path.dirname(__file__))
    output_fq = os.path.join(current_path, args.o)
    tag5 = args.tag_5
    tag3 = args.tag_3

    if tag5 is None and tag3 is None:
        print("The tag sequence must provide!")
        sys.exit(1)
    elif tag5 is None:
        command = f'cutadapt -a {tag3.upper()}  -j 0 -n 5 -o {output_fq} {input_fq}'
        print(command)
        print(f'Running {command} ...')
    elif tag3 is None:
        command = f'cutadapt -g {tag5.upper()}  -j 0 -n 5 -o {output_fq} {input_fq}'
        print(command)
        print(f'Running {command} ...')
    else:
        command = f'cutadapt -g {tag5.upper()} -a {tag3.upper()}  -j 0 -n 5 -o {output_fq} {input_fq}'
        print(command)
        print(f'Running {command} ...')

    # run_shell(command)


if __name__ == '__main__':
    main()
