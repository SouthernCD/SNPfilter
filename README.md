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
usage: SNPfilter prepare [-h] [-t THREADS] [-q MIN_MQ] [-Q MIN_BQ] sample_id reference R1 R2

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
  -q MIN_MQ, --min-MQ MIN_MQ
                        skip alignments with mapQ smaller than INT
  -Q MIN_BQ, --min-BQ MIN_BQ
                        skip bases with baseQ/BAQ smaller than INT
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

**codefilter**: filtering SNPs based on whether they cause changes in coding amino acids

```
usage: SNPfilter codefilter [-h] sample_id input_bcf reference gff_file

filtering SNPs based on whether they cause changes in coding amino acids

positional arguments:
  sample_id   sample id
  input_bcf   input bcf file
  reference   reference genome in fasta format
  gff_file    gff file for reference genome

optional arguments:
  -h, --help  show this help message and exit
```

## Example
1. Indexing of the reference genome using BWA
```
bwa index reference.fasta
```

2. Trim the sequencing data 
```
zcat sample1_R1.fastq.gz | bioawk -cfastx '{print "@"$name"\n"substr($seq, 16)"\n+"$name"\n"substr($qual, 16)}' | gzip > sample1_R1.trimmed.fastq.gz
zcat sample1_R2.fastq.gz | bioawk -cfastx '{print "@"$name"\n"substr($seq, 16)"\n+"$name"\n"substr($qual, 16)}' | gzip > sample1_R2.trimmed.fastq.gz
zcat sample2_R1.fastq.gz | bioawk -cfastx '{print "@"$name"\n"substr($seq, 16)"\n+"$name"\n"substr($qual, 16)}' | gzip > sample2_R1.trimmed.fastq.gz
zcat sample2_R2.fastq.gz | bioawk -cfastx '{print "@"$name"\n"substr($seq, 16)"\n+"$name"\n"substr($qual, 16)}' | gzip > sample2_R2.trimmed.fastq.gz
zcat sample3_R1.fastq.gz | bioawk -cfastx '{print "@"$name"\n"substr($seq, 16)"\n+"$name"\n"substr($qual, 16)}' | gzip > sample3_R1.trimmed.fastq.gz
zcat sample3_R2.fastq.gz | bioawk -cfastx '{print "@"$name"\n"substr($seq, 16)"\n+"$name"\n"substr($qual, 16)}' | gzip > sample3_R2.trimmed.fastq.gz
```

3. Use **prepare** to generate bam files and BCF files
```
SNPfilter prepare -t 10 sample1 reference.fasta sample1_R1.trimmed.fastq.gz sample1_R2.trimmed.fastq.gz
SNPfilter prepare -t 10 sample2 reference.fasta sample2_R1.trimmed.fastq.gz sample2_R2.trimmed.fastq.gz
SNPfilter prepare -t 10 sample3 reference.fasta sample3_R1.trimmed.fastq.gz sample3_R2.trimmed.fastq.gz
```

4. Use **qcfilter** to filter SNPs with a relaxed threshold（d=2, v=0.3)
```
SNPfilter qcfilter -d 2 -v 0.3 sample1 sample1.bcf
SNPfilter qcfilter -d 2 -v 0.3 sample2 sample2.bcf
SNPfilter qcfilter -d 2 -v 0.3 sample3 sample3.bcf
```

5. Use **qcfilter** to filter SNPs with a strict threshold（d=5, v=0.9, other sample BCF file as background)
```
SNPfilter qcfilter -d 5 -v 0.9 -b sample2.d2.v0.30.bcf,sample3.d2.v0.30.bcf sample1 sample1.d2.v0.30.bcf
SNPfilter qcfilter -d 5 -v 0.9 -b sample1.d2.v0.30.bcf,sample3.d2.v0.30.bcf sample2 sample2.d2.v0.30.bcf
SNPfilter qcfilter -d 5 -v 0.9 -b sample1.d2.v0.30.bcf,sample2.d2.v0.30.bcf sample3 sample3.d2.v0.30.bcf
```

6. Use **codefilter** to filter SNPs that do not cause changes in coding amino acids
```
SNPfilter codefilter sample1 sample1.d5.v0.90.bcf reference.fasta reference.gff
SNPfilter codefilter sample2 sample2.d5.v0.90.bcf reference.fasta reference.gff
SNPfilter codefilter sample3 sample3.d5.v0.90.bcf reference.fasta reference.gff
```
