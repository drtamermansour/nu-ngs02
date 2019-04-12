# VCF-tools

## Installation

```
conda activate ngs1
conda install vcftools
conda install bcftools
```

## Test Data Download

### Small sample

```bash
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/release/2010_07/exon/snps/YRI.exon.2010_03.sites.vcf.gz -O pilot.vcf.gz

wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/release/2010_07/exon/snps/YRI.exon.2010_03.genotypes.vcf.gz -O gt.vcf.gz
```

### Extract both and keep the compressed file

`gunzip -k *gz`

---

## Get basic file statistics

### Working with the compressed `--gzvcf` flag

`vcftools --gzvcf pilot.vcf`

> **Try to figure out what is the filtered & unfiltered sites!**

### Doing the same from stdin on the uncompressed file

`cat pilot.vcf | vcftools --vcf -`

> Don't forget the last "-" char for receiving from the stdin.

---

## Applying a filter

`vcftools --vcf pilot.vcf --chr 1 --from-bp 1000000 --to-bp 2000000`

> VCFtools options are very clear right?
> **--chr** > For the #CHROM field
> **--from-pb** starting from that position to the **--to--pb** position.

### Let's write the filtered records in a new VCF.

`vcftools --vcf pilot.vcf --chr 1 --from-bp 1000000 --to-bp 2000000 --recode --recode-INFO-all --out F1`

>  **--recode-INFO-all** to include all data from the INFO fields in the output.
> **--recode** To write out the variants that pass through filters.

### piping the output to the stdout

`vcftools --vcf pilot.vcf --chr 1 --from-bp 1000000 --to-bp 2000000 --recode --recode-INFO-all --stdout`

### Write a compressed output (useful for large data)

`vcftools --vcf pilot.vcf --chr 1 --from-bp 1000000 --to-bp 2000000 --recode --recode-INFO-all -c | bgzip -c > subset.vcf.gz`

> **-c** == **--stdout**

> ### You got an error? try to find which conda package is missing and install it then re-run the previous command.

---

## Getting allele frequency

`vcftools --vcf gt.vcf --freq --out freq`

---

## Getting sequencing depth information

`vcftools --vcf gt.vcf --depth -c > gt_depth_summary.txt`

