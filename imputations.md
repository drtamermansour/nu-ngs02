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

|      param     |               Description               |
|:--------------:|:---------------------------------------:|
|       gt       |        Input VCF file (required)        |
|       ref      | bref3 or VCF file with phased genotypes |
|       out      |      output file prefix (required)      |
|       map      |       PLINK map file with cM units      |
|      chrom     |     [chrom] or [chrom]:[start]-[end]    |
| excludesamples |      file with 1 sample ID per line     |
| excludemarkers |      file with 1 marker ID per line     |

### Imputation parameters

|     Param    |               Description               | Default |
|:------------:|:---------------------------------------:|:-------:|
|    impute    | impute ungenotyped markers (true/false) |   True  |
|  imp-states  |       model states for imputation       |   1600  |
|  imp-segment |    min haplotype segment length (cM)    |   6.0   |
| imp-clusters |        max cM in a marker cluster       |  0.005  |
|    imp-ap    |   print posterior allele probabilities  |  false  |
|    imp-gp    |  print posterior genotype probabilities |  false  |

### General parameters

| param    | description               | default           | 
|----------|---------------------------|-------------------| 
| ne       | effective population size | 1000000           | 
| err      | allele mismatch rate      | 1.0E-4            | 
| window   | window length in cM       | 40.0              | 
| overlap  | window overlap in cM      | 4.0               | 
| seed     | random seed               | -99999            | 
| nthreads | number of threads         | machine-dependent | 
| step     | IBS step length (cM       | 0.1               | 
| nsteps   | number of IBS steps       | 7                 | 











