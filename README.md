# SNPfilter
A handy little tool for filtering SNPs

## Installation
```
conda install -c bioconda samtools bwa
pip install snpfilter
```

## Usage

**prepare**: Use BWA and SAMTools to map the sequencing data, generate BAM files and BCF files, etc.

```
usage: SNPfilter prepare [-h] [-t THREADS] sample_id reference R1 R2

Prepare work environment

positional arguments:
  sample_id             sample id
  reference             reference genome in fasta format, should be indexed by bwa
  R1                    fastq file 1
  R2                    fastq file 2

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        num of threads
```

**qcfilter**: filtering SNPs from BCF file with min depth and min variant frequency

```
usage: SNPfilter qcfilter [-h] [-d MIN_DEPTH] [-v MIN_VARIANT_FREQUENCY] [-b BACKGROUND_SITE] sample_id input_bcf

filtering snp from bcf file with min depth and min variant frequency

positional arguments:
  sample_id             sample id
  input_bcf             input bcf file

optional arguments:
  -h, --help            show this help message and exit
  -d MIN_DEPTH, --min_depth MIN_DEPTH
                        min depth
  -v MIN_VARIANT_FREQUENCY, --min_variant_frequency MIN_VARIANT_FREQUENCY
                        min variant frequency
  -b BACKGROUND_SITE, --background_site BACKGROUND_SITE
                        A comma-separated list of bcf files, loci that appear in these bcf files will be filtered out
```

## Example
1. Indexing of the reference genome using BWA
```
bwa index reference.fasta
```

2. Use **prepare** to generate bam files and BCF files
```
SNPfilter prepare -t 10 sample1 reference.fasta sample1_R1.fastq sample1_R2.fastq
SNPfilter prepare -t 10 sample2 reference.fasta sample2_R1.fastq sample2_R2.fastq
SNPfilter prepare -t 10 sample3 reference.fasta sample3_R1.fastq sample3_R2.fastq
```

3. Use **qcfilter** to filter SNPs with a relaxed threshold（d=2, v=0.3)
```
SNPfilter qcfilter -d 2 -v 0.3 sample1 sample1.bcf
SNPfilter qcfilter -d 2 -v 0.3 sample2 sample2.bcf
SNPfilter qcfilter -d 2 -v 0.3 sample3 sample3.bcf
```

4. Use **qcfilter** to filter SNPs with a strict threshold（d=5, v=0.5, other sample BCF file as background)
```
SNPfilter qcfilter -d 5 -v 0.5 -b sample2.d2.v0.30.bcf,sample3.d2.v0.30.bcf sample1 sample1.d2.v0.30.bcf
SNPfilter qcfilter -d 5 -v 0.5 -b sample1.d2.v0.30.bcf,sample3.d2.v0.30.bcf sample2 sample2.d2.v0.30.bcf
SNPfilter qcfilter -d 5 -v 0.5 -b sample1.d2.v0.30.bcf,sample2.d2.v0.30.bcf sample3 sample3.d2.v0.30.bcf
```
