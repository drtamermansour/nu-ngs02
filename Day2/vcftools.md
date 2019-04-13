# VCF-tools

## Installation

```
conda activate ngs1
conda install vcftools
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

`vcftools --gzvcf pilot.vcf.gz`

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

```
CHROM	POS	N_ALLELES	N_CHR	{ALLELE:FREQ}
1	1105324	2	126	C:0.968254	T:0.031746
```

> The N_ALLELES column indicates the number of possible alleles at that locus.
> In our case it is 2 at both sites, the reference and one alternate allele. 
> The N_CHR column indicates the amount of chromosomes you have data available for at that locus. In our case, we have data for 63 diploid individuals at each site, giving us 126 chromosomes.

---

## Extract specific SNPs

extract SNP: rs725021

`vcftools --vcf gt.vcf --snp rs725021 --recode --recode-INFO-all --stdout > rs725021.vcf`

### Extract list of SNPs

`nano snps.txt`

```
rs725021
rs34506306
rs35595233
rs1382603
rs34752670
```

`vcftools --vcf gt.vcf --snps snps.txt --recode --recode-INFO-all --stdout > snps_list.vcf`

## FILTER FLAG FILTERING

> Removes all sites with a FILTER flag other than PASS.

`vcftools --vcf gt.vcf --remove-filtered-all --recode --recode-INFO-all --stdout > all_passed.vcf`

## Export all genotypes depth

`vcftools --vcf gt.vcf --geno-depth -c > gt_genotypes_depth.vcf`

## Generates a file containing the mean depth per individual. This file has the suffix.

`vcftools --vcf gt.vcf --depth -c > mean_depth.vcf`

## the depth  for each site summed across individuals.

`vcftools --vcf gt.vcf --site-depth -c > site_depth.vcf`

## Alleles Count
`vcftools --vcf gt.vcf --counts -c > gt_counts.vcf`

---

## Convert to PLINK

`vcftools --vcf gt.vcf --plink  --out gt_plink`

