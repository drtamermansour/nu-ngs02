# Beagle software

> Beagle is a software package for phasing genotypes and for imputing ungenotyped markers (missing genotypes.)

## Check Java installation

`java -version`

if not installed type:

```bash
conda activate ngs1
conda install java-jdk
```

---

## Software Parameters

### Data Parameters

|      param     	|               Description               	|
|:--------------:	|:---------------------------------------:	|
|       gt       	|        Input VCF file (required)        	|
|       ref      	| bref3 or VCF file with phased genotypes 	|
|       out      	|      output file prefix (required)      	|
|       map      	|       PLINK map file with cM units      	|
|      chrom     	|     [chrom] or [chrom]:[start]-[end]    	|
| excludesamples 	|      file with 1 sample ID per line     	|
| excludemarkers 	|      file with 1 marker ID per line     	|







