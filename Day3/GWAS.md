Genome Wide Association analysis (GWAS)
=======================================

Install [Plink](https://www.cog-genomics.org/plink)
```
conda activate ngs1
conda install -c bioconda plink 
```

Download the sample VCF file and phenotype data¶
```
mkdir ~/GWAS && cd ~/GWAS
wget https://de.cyverse.org/dl/d/E0A502CC-F806-4857-9C3A-BAEAA0CCC694/pruned_coatColor_maf_geno.vcf.gz
gunzip pruned_coatColor_maf_geno.vcf.gz
wget https://de.cyverse.org/dl/d/3B5C1853-C092-488C-8C2F-CE6E8526E96B/coatColor.pheno
```

Convert your VCF into Plink format
```
#conda install vcftools

```

