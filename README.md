# cvmtrans

![PYPI](https://img.shields.io/pypi/v/cvmtrans)
![Static Badge](https://img.shields.io/badge/OS-_Mac_%7C_Linux-steelblue)


## Installation

pip3 install cvmtrans


## Usage

### 1. Extract reads contains specific tags

```
Usage: cvmtrans_extract_reads -i <input_fastq> -tag <sequence should inclued in reads> -o <output_fastq>

Author: Qingpo Cui(SZQ Lab, China Agricultural University)

optional arguments:
  -h, --help  show this help message and exit
  -i I        <input_fastq>: output filename of the fastq file
  -tag TAG    <input_file>: the tag sequence
  -o O        <output_fastq>: output filename of the fastq file
```

### 2. Cut 5' end or 3' end tags
```
Usage: cvmtrans_cut_tags -i <input_fastq> -tag <sequence should inclued in reads> -o <output_fastq>

Author: Qingpo Cui(SZQ Lab, China Agricultural University)

optional arguments:
  -h, --help    show this help message and exit
  -i I          <input_fastq>: output filename of the fastq file
  -tag_5 TAG_5  <5‘ end tags or barcode>: The tag sequence at the 5’ end will be removed
  -tag_3 TAG_3  <3‘ end tags or barcode>: The tag sequence at the 3’ end will be removed
  -o O          <output_fastq>: output filename of the fastq file
```


### 3. Get the insertion site statistics table
```
Usage: cvmtrans -bam_file <sorted_bam_files> -embl_file <reference embl file> -output <tab delimited text file>

Author: Qingpo Cui(SZQ Lab, China Agricultural University)

optional arguments:
  -h, --help            show this help message and exit
  -bam_file BAM_FILE    <input_fastq>: output filename of the fastq file
  -embl_file EMBL_FILE  <input_file>: the tag sequence
  -output OUTPUT        <output_fastq>: output filename of the fastq file
```


