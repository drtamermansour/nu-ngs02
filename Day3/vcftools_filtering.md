
# VCFtools Filtering Tutorial

## Installation

```
conda activate ngs1
conda install -y vcftools
```

## Data Preparation

```
mkdir -p workdir/vcftools_tut && cd workdir/vcftools_tut

# Download data.zip file
curl -L -o data.zip https://www.dropbox.com/sh/bf9jxviaoq57s5v/AAD2Kv5SPpHlZ7LC7sBz4va8a?dl=1

# Extract
unzip data.zip

```

## Tutorial

### Command #1

```
vcftools --gzvcf raw.vcf.gz --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out raw.g5mac3
```
- We call **vcftools**, feed it a vcf file after the `--vcf` flag.
- The `--max-missing 0.5` tells it to filter genotypes called *below 50%* (across all individuals).
- The `--mac 3` flag tells it to filter SNPs that have a minor allele count *less than 3.* this is relative to genotypes, so - it has to be called in at least **1 homozygote and 1 heterozygote** or **3 heterozygotes**.
- The `--recode` flag tells the program to write a new vcf file with the filters.
- The `--recode-INFO-all` keeps all the INFO flags from the old vcf file in the new one.
- Lastly, `--out` designates the name of the output, The output will scroll through a lot of lines, but should end like:


```
After filtering, kept 40 out of 40 Individuals
Outputting VCF file...
After filtering, kept 78434 out of a possible 147540 Sites
Run Time = 15.00 seconds
```
