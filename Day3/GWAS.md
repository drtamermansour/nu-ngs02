Genome Wide Association analysis (GWAS)
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
```

convert VCF into Plink readable format (map,ped) then Plink binary format (fam,bed,bim)
```
#conda install vcftools
vcftools --vcf pruned_coatColor_maf_geno.vcf --plink --out coatColor
plink --file coatColor --allow-no-sex --dog --make-bed --noweb --out coatColor.binary
```

