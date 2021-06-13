Genome Wide Association Analysis (GWAS)
=======================================

Install [Plink](https://www.cog-genomics.org/plink). 

**Notes:**
1.  Shouldn't we use PLINK2 instead of the old PLINK? Actually, things are not that simple. You can read more [here](https://www.biostars.org/p/299855/) about this
2.  Unfortunately, Pink documentation is kind of cumulative. You need always to start from the [documentation of Plink 1.07](https://zzz.bwh.harvard.edu/plink/)    
```
conda activate ngs1
conda install -c bioconda plink 
```

Download the sample VCF file and phenotype data: Genotyping of 476840 SNPs in 53 dogs (24 yellow coat and 29 dark coat)
```
mkdir -p ~/GWAS && cd ~/GWAS
wget https://de.cyverse.org/dl/d/E0A502CC-F806-4857-9C3A-BAEAA0CCC694/pruned_coatColor_maf_geno.vcf.gz
gunzip pruned_coatColor_maf_geno.vcf.gz
wget https://de.cyverse.org/dl/d/3B5C1853-C092-488C-8C2F-CE6E8526E96B/coatColor.pheno
## explore the headers line of the VCF
grep "^#CHROM" pruned_coatColor_maf_geno.vcf
## explore the phenotype file
less coatColor.pheno  ## phenotype file: family IDs, within-family IDs, and phenotype
```

convert VCF into Plink readable format (map,ped) then Plink binary format (fam,bed,bim)
```
#conda install vcftools
vcftools --vcf pruned_coatColor_maf_geno.vcf --plink --out coatColor
plink --file coatColor --allow-no-sex --dog --make-bed --out coatColor.binary
```

[map file format](https://www.cog-genomics.org/plink/1.9/formats#map)
A text file with no header file, and one line per variant with the following fields:

*  Chromosome code. 
*  Variant identifier
*  Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
*  Base-pair coordinate

[ped file format](https://www.cog-genomics.org/plink/1.9/formats#ped)
Contains no header line, and one line per sample with 2V+6 fields where V is the number of variants. The first six fields are:

*  Family ID ('FID')
*  Within-family ID ('IID'; cannot be '0')
*  Within-family ID of father ('0' if father isn't in dataset)
*  Within-family ID of mother ('0' if mother isn't in dataset)
*  Sex code ('1' = male, '2' = female, '0' = unknown)
*  Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)


**Run a simple association analysis**

  *  **Note1:** –-assoc performs a standard case/control association analysis which is a chi-square test of allele frequency.
  
  *  **Note2:** If the ped file did not have informative phenotype values or we want to replace these phenotype values, --make-pheno <phenotype_file> <phenotype> can be used 

  *  **Note3:** By default, the minor allele is coded A1 and tested for being the risk allele. This will be confusing if the reference allele happens to be the minor allele. --a2-allele <filename> [a2col] [IDcol] [skip] : Force alleles in the "a2col" column of the file to A2.

  *  **Note4:** –-adjust enables correction for multiple analysis and automatically calculates the genomic inflation factor

```
plink --bfile coatColor.binary --assoc --make-pheno coatColor.pheno "yellow" --a2-allele pruned_coatColor_maf_geno.vcf 4 3 '#' --allow-no-sex --adjust --dog --out coatColor
```

check the output files!

**Create Manhattan plot**

Install qqman package
```
Rscript -e "install.packages('qqman',  contriburl=contrib.url('http://cran.r-project.org/'))"
```

Identify statistical cutoffs
```
unad_cutoff_sug=$(tail -n+2 coatColor.assoc.adjusted | awk '$5>=0.05' | head -n1 | awk '{print $3}')
unad_cutoff_conf=$(tail -n+2 coatColor.assoc.adjusted | awk '$5>=0.01' | head -n1 | awk '{print $3}')
```

Run the plotting function
```
Rscript -e 'args=(commandArgs(TRUE));library(qqman);'\
'data=read.table("coatColor.assoc", header=TRUE); data=data[!is.na(data$P),];'\
'pdf("coatColor_man.pdf", width=20, height=10);'\
'manhattan(data, p = "P", col = c("blue4", "orange3"),'\
'suggestiveline = -log(as.numeric(args[1]),10),'\
'genomewideline = -log(as.numeric(args[2]),10),'\
'chrlabs = c(1:38, "X"), annotateTop=TRUE, cex = 1.2);'\
'graphics.off();' $unad_cutoff_sug $unad_cutoff_conf
```
Open the output file in the browser
The top associated mutation is a nonsense SNP in MC1R (c.916C>T) known to control pigment production. 

