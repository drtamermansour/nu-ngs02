Genome Wide Association Analysis (GWAS)
=======================================

Install [Plink](https://www.cog-genomics.org/plink)
```
conda activate ngs1
conda install -c bioconda plink 
```

Download the sample VCF file and phenotype data: Genotyping of 476840 SNPs in 53 dogs (24 yellow coat and 29 dark coat)
```
mkdir ~/GWAS && cd ~/GWAS
wget https://de.cyverse.org/dl/d/E0A502CC-F806-4857-9C3A-BAEAA0CCC694/pruned_coatColor_maf_geno.vcf.gz
gunzip pruned_coatColor_maf_geno.vcf.gz
wget https://de.cyverse.org/dl/d/3B5C1853-C092-488C-8C2F-CE6E8526E96B/coatColor.pheno
## explore the headers line of the VCF
grep "^#CHROM" pruned_coatColor_maf_geno.vcf
## explore the phenotype file
less coatColor.pheno
```

convert VCF into Plink readable format (map,ped) then Plink binary format (fam,bed,bim)
```
#conda install vcftools
vcftools --vcf pruned_coatColor_maf_geno.vcf --plink --out coatColor
plink --file coatColor --allow-no-sex --dog --make-bed --noweb --out coatColor.binary
```

**Run a simple association analysis**

  *  **Note1:** –assoc performs a standard case/control association analysis which is a chi-square test of allele frequency.

  *  **Note2:** –reference-allele allow you to use your list of A1 alleles. By default, the minor allele is coded A1 and tested for being the risk allele. This will be confusing if the reference allele happens to be the minor allele. So we generate a list of alternative alleles and used it with –reference-allele option to be used as A1  

  *  **Note3:** –adjust enables correction for multiple analysis and automatically calculates the genomic inflation factor

```
#create list of alternative alleles
cat pruned_coatColor_maf_geno.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}'  > alt_alleles
plink --bfile coatColor.binary --make-pheno coatColor.pheno "yellow" --assoc --reference-allele alt_alleles --allow-no-sex --adjust --dog --noweb --out coatColor
```
