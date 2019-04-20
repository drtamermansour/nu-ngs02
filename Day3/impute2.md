# [Impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)


## Download

```bash
mkdir -p ~/impute && cd ~/impute
wget -c https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz
tar -zxvf impute_v2.3.2_x86_64_static.tgz
```


## Example (IMPUTATION WITH ONE PHASED REFERENCE PANEL)

> *impute untyped SNPs in a study dataset from a panel of reference haplotypes.*

| option                | description                                                                                                                            | 
|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------| 
| [-m](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-m) <file>             | Fine-scale recombination map for the region to be analyzed. This file should have three columns: physical position (in base pairs), recombination rate between current position and next position in map (in cM/Mb), and genetic map position (in cM). The file should also have a header line with an unbroken character string for each column (e.g., "position COMBINED_rate(cM/Mb) Genetic_Map(cM)").       | 
| [-h](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-h)  <file_1>  <file_2>  | File of known haplotypes, with one row per SNP and one column per haplotype. All alleles must be coded as 0 or 1, and each -h file must be provided with a corresponding legend file. impute2 provides formatted haplotypes from the HapMap Project and the 1,000 Genomes Project. Also it has a script to format custom phased VCF into reference panel format                      | 
| [-l](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-l) <file_1> <file_2>  | Legend file(s) with information about the SNPs in the -h file(s).                                                                      | 
| [-g](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-g) <file>             | File containing genotypes for a study cohort that you want to impute or phase. The gen file has the format of ID1,ID2,POS,A,B followed by three genotype probabilities P(AA), P(AB), P(BB) for each sample. Impute2 has a script to convert a phased or unphased VCF file into genotype file format                                                        | 
| -[strand_g](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-strand_g) <file>      | File showing the strand orientation of the SNP allele codings in the -g file                                                           | 
| [-int](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-int) <lower> <upper>  | Genomic interval to use for inference, as specified by <lower> and <upper> boundaries in base pair position                            | 
| [-Ne](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-ne) <int>             | "Effective size" of the population (commonly denoted as Ne in the population genetics literature) from which your dataset was sampled. | 
| [-o](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-o)                    | Name of main output file.                                                                                                              | 



```bash
cd impute_v2.3.2_x86_64_static/

./impute2 \
 -m ./Example/example.chr22.map \
 -h ./Example/example.chr22.1kG.haps \
 -l ./Example/example.chr22.1kG.legend \
 -g ./Example/example.chr22.study.gens \
 -strand_g ./Example/example.chr22.study.strand \
 -int 20.4e6 20.5e6 \
 -Ne 20000 \
 -o ./Example/example.chr22.one.phased.impute2
```


### [Output Description](http://mathgen.stats.ox.ac.uk/impute/concordance_table_description.html)

Let us check what happened

a) Recover the original VCF file before imputation 
```
plink --gen Example/example.chr22.study.gens --sample Example/example.study.samples --oxford-single-chr 1 --recode vcf --out Example/example.chr22.study.gens.plink
```
b) Convert the final imput2 output to VCF as well
```
plink --gen Example/example.chr22.one.phased.impute2 --sample Example/example.study.samples --oxford-single-chr 1 --recode vcf --out Example/example.chr22.one.phased.impute2.plink
```

Can you check the log files of the Plink processes?
Check for the no of variants and the genotyping rate?
