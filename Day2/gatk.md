
GATK Variant calling
====================

> **What is the required data coverage?**

Did you run the [crash variant calling tutorial](https://github.com/drtamermansour/nu-ngs02/blob/master/Crash_variant_calling.md)? Make sure you run it to download data and software needed for this tutorial

Data trimming: 
- Trimming is data loss so be careful.
- Sequence trimming is complementary to variant filtration
- Sources of errors: 
    * The call is suspicious ==> low quality score (variant filtration is better than quality trimming) 
    * Technical problems (e.g. sequencing chemistry or physics) ==> systematic errors (can be removed by careful kmer based trimming. GATK recalibration is an alternative)
    
- Very mild quality trimming: SLIDINGWINDOW:4:2 ==> this means that the Base call accuracy is ~ 40%. Check this [tutorial](https://github.com/drtamermansour/nu-ngs01/blob/master/Day-2/Trimmomatic_tutorial.md) for more info.

## Explore the sample names
```
ls -tral ~/workdir/fqData/*_R*_001.pe.fq.gz
```

## Add [Read group information] and align all reads
```
mkdir -p ~/workdir/GATK_tutorial && cd ~/workdir/GATK_tutorial
for R1 in ~/workdir/fqData/*_R1_001.pe.fq.gz;do
    SM=$(basename $R1 | cut -d"_" -f1)                                          ##sample ID
    LB=$(basename $R1 | cut -d"_" -f1,2)                                        ##library ID
    PL="Illumina"                                                           ##platform (e.g. illumina, solid)
    RGID=$(zcat $R1 | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)       ##read group identifier 
    PU=$RGID.$LB                                                            ##Platform Unit
    echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

    R2=$(echo $R1 | sed 's/_R1_/_R2_/')
    echo $R1 $R2
    index="$HOME/workdir/bwa_align/bwaIndex/dog_chr5.fa"
    bwa mem -t 4 -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" $index $R1 $R2 > $(basename $R1 _R1_001.pe.fq.gz).sam
done
```

Note:
```
a) Sample and library tags
    - SM = biological sample name.
    - LB = name of DNA preparation library tube = {SM}.{library-specific identifier}(Important To identify PCR duplicates in MarkDuplicates step. Ignore in PCR free libraries)

Can be autmatically detected from current sample naming scheme:<Sample.ID><Index.Sequence><Lane.ID><Set.number>.fastq

    - SM = <Sample.ID>
    - LB = <Sample.ID>_<Index.Sequence>

b) ID and PU (to enable merging replictes)
    - ID = Read group identifier = {FLOWCELL_BARCODE}.{LANE}
    - PU = Platform Unit = {FLOWCELL_BARCODE}.{LANE}.{library-specific identifier}. This is the most specific definition for a group of reads.

Also can be identified from the name of a sequence read in the Fastq file:@(instrument id):(run number):(flowcell ID):(lane):(tile):(x_pos):(y_pos) (read):(is filtered):(control number):(index sequence)FLOWCELL_BARCODE = @(instrument id):(run number):(flowcell ID)

special Notes:
- One sample (SM) can have multiple libraries (e.g SE, PE50, and PE100) (LB), can run on multiple lanes and/or multiple flow cells (RGID), and can run on multiple platforms (PL).
- One library can run on multiple lanes or multiple flow cells (PU).
- If we have multiple samples for the same individual e.g. before and after treatment, each sample should have a different SM but unless you expect change of sequence, we can consider them multiple libraries of the same sample)
- Multiple samples can share the same Read group ID (When manuals sya “must be unique”. They mean unique in a BAM file. So it is ok that multiple samples can share the same Read group ID)
- if you have one library for each sample running on one lane of a sequencing machine then you can make SM=LB=RGID=PU

```

## generate & sort BAM file

```
for samfile in *.sam;do
  sample=${samfile%.sam}
  samtools view -hbo $sample.bam $samfile
  samtools sort $sample.bam -o $sample.sorted.bam
done
```
Explore files size!

## Merge replicates (e.g. one library running on two lanes) with [Picard tools](http://broadinstitute.github.io/picard/):
```
# Install Picard tools
conda install -c bioconda picard 
picard_path=$CONDA_PREFIX/share/picard-2.19.0-0


# merge the replicates
java  -Xmx2g -jar $picard_path/picard.jar MergeSamFiles I=BD143_TGACCA_L005.sorted.bam I=BD143_TGACCA_L006.sorted.bam OUTPUT=BD143_TGACCA_merged.sorted.bam

# check for the changes in the header
samtools view -H BD143_TGACCA_L005.sorted.bam
samtools view -H BD143_TGACCA_L006.sorted.bam
samtools view -H BD143_TGACCA_merged.sorted.bam

# remove the individual replicates
rm BD143_TGACCA_L00*.sorted.bam
```

**Note**: Duplicate marking should NOT be applied to amplicon sequencing data or other data types where reads start and stop at the same positions by design.


## mapping QC
```
for bamFile in *.sorted.bam;do
  output=${bamFile%.sorted.bam}
  samtools depth $bamFile | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' > $output.cov
  samtools flagstat $bamFile > $output.stat
done
```

## Mark duplicate
```
for sample in *.sorted.bam;do
  name=${sample%.sorted.bam}
  java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
done
```

## Install GATK
```
conda install -c bioconda gatk4 
```
check for best practice [here](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145)

## indexing
```
# samples
for sample in *.dedup.bam;do
  #name=${sample%.dedup.bam}
  java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=$sample
done

# Reference
ln -s ~/workdir/sample_data/dog_chr5.fa .
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=dog_chr5.fa O=dog_chr5.dict
samtools faidx dog_chr5.fa
```

## Download known varinats

```
# Download known polymorphic sites
wget 'ftp://ftp.ensembl.org/pub/release-89/variation/vcf/canis_familiaris/Canis_familiaris.vcf.gz' -O canis_familiaris.vcf.gz

# Select variants on chr5 and correct chr name
gunzip canis_familiaris.vcf.gz
grep "^#" canis_familiaris.vcf > canis_fam_chr5.vcf
grep "^5" canis_familiaris.vcf | sed 's/^5/chr5/' >> canis_fam_chr5.vcf
gatk IndexFeatureFile -F canis_fam_chr5.vcf
```

Note the differences between genome annotation databases. Not only chromosome names but more importantly the coordinate system [interesting post](https://www.biostars.org/p/84686/)

## Recalibrate Bases [BQSR](https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr)

- Data pre-processing step that detects systematic errors made by the sequencer when it estimates the quality score of each base call.
- Various sources of systematic (non-random) technical error, leading to over- or under-estimated base quality scores in the data. Some of these errors are due to the physics or the chemistry of how the sequencing reaction works, and some are probably due to manufacturing flaws in the equipment.
- We apply machine learning to model these errors empirically and adjust the quality scores accordingly. For example we can identify that, for a given run, whenever we called two A nucleotides in a row, the next base we called had a 1% higher rate of error. So any base call that comes after AA in a read should have its quality score reduced by 1%.
- We do that over several different covariates (mainly sequence context and position in read, or cycle) in a way that is additive. So the same base may have its quality score increased for one reason and decreased for another
- The base recalibration process involves two key steps: first one tool [(BaseRecalibrator)](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php#--use-original-qualities) builds a model of covariation based on the data and a set of known variants, then another tool [(ApplyBQSR)](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php) adjusts the base quality scores in the data based on the model.

```
for sample in *.dedup.bam;do
  name=${sample%.dedup.bam}

  gatk --java-options "-Xmx2G" BaseRecalibrator \
-R dog_chr5.fa -I $sample --known-sites canis_fam_chr5.vcf \
-O $name.report

  gatk --java-options "-Xmx2G" ApplyBQSR \
-R dog_chr5.fa -I $sample -bqsr $name.report \
-O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals
done
```

## Joint variant calling using [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php)

Call germline SNPs and indels via **local re-assembly** of haplotypes

```
## assess genotype likelihood per-sample
for sample in *.bqsr.bam;do
  name=${sample%.bqsr.bam}

  gatk --java-options "-Xmx2G" HaplotypeCaller \
  -R dog_chr5.fa -I $sample \
  --emit-ref-confidence GVCF \
  --pcr-indel-model NONE \
  -O $name.gvcf
done

## combine samples
gatk --java-options "-Xmx2G" CombineGVCFs \
-R dog_chr5.fa \
-V BD143_TGACCA_merged.gvcf \
-V BD174_CAGATC_L005.gvcf \
-V BD225_TAGCTT_L007.gvcf \
-O raw_variants.gvcf

## Joint Genotyping
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R dog_chr5.fa \
-V raw_variants.gvcf \
--max-alternate-alleles 6 \
-O raw_variants.vcf

## annotated output
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R dog_chr5.fa \
-V raw_variants.gvcf \
--max-alternate-alleles 6 \
--dbsnp canis_fam_chr5.vcf \
-O raw_variants_ann.vcf

## check how many variant got annotated
grep -v "^#" raw_variants_ann.vcf | awk '{print $3}' | grep "^rs" | wc -l
```

## VCF statitics

First letus index the VCF file
```
conda install -c bioconda tabix
bgzip -c raw_variants_ann.vcf > raw_variants_ann.vcf.gz
tabix -p vcf raw_variants_ann.vcf.gz
```

Calc some stats about your vcf
```
conda install -c bioconda rtg-tools
rtg vcfstats raw_variants_ann.vcf.gz > stats.txt
```

read more about [RTG tools](https://www.realtimegenomics.com/products/rtg-tools) and explore there [manual](https://cdn.rawgit.com/RealTimeGenomics/rtg-tools/master/installer/resources/tools/RTGOperationsManual/index.html) for the possible commands

## Split SNPs and indels

```
gatk --java-options "-Xmx2G" SelectVariants \
-R dog_chr5.fa \
-V raw_variants_ann.vcf \
--select-type-to-include SNP \
-O raw_variants_ann_SNP.vcf

gatk --java-options "-Xmx2G" SelectVariants \
-R dog_chr5.fa \
-V raw_variants_ann.vcf \
--select-type-to-include INDEL \
-O raw_variants_ann_INDEL.vcf
```

## Assess the different filters in both known and novel
```
for var in "SNP" "INDEL";do
 input="raw_variants_ann_"$var".vcf"
 for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AN" "DP" "InbreedingCoeff";do
  filterValues=$var.$filter
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > $filterValues
  grep -v "^\." $filterValues > known.$var.$filter
  grep "^\." $filterValues > novel.$var.$filter
done; done
```

let us make things a little bit more organized
```
mkdir filters && cd filters
mv ../{*.SNP.*,SNP.*,*.INDEL.*,INDEL.*} .
```

Figures!!!
```
wget https://raw.githubusercontent.com/dib-lab/dogSeq/master/scripts/densityCurves.R
sudo Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"
for f in SNP.* INDEL.*;do
  Rscript densityCurves.R "$f"
done
```

Calc the DP threathols
```
cat SNP.DP INDEL.DP | awk '{sum+= $2; sumsq+= ($2)^2} END { print sum/NR, sqrt((sumsq-sum^2/NR)/NR), sum/NR + 5*sqrt((sumsq-sum^2/NR)/NR) }' 
```

## SNP Variant filteration
```
cd ~/workdir/GATK_tutorial
gatk --java-options "-Xmx2G" VariantFiltration \
-R dog_chr5.fa \
-V raw_variants_ann_SNP.vcf \
--filter-name "snpQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "snpMQ" \
--filter-expression "vc.hasAttribute('MQ') && MQ < 40.0" \
--filter-name "snpMQRankSum" \
--filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
--filter-name "snpFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 60.0" \
--filter-name "snpSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 4.0" \
--filter-name "snpReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
--filter-name "snpDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O raw_variants_ann_SNP_clean.vcf
```

check the filtered recored 
```
grep -v "^#" raw_variants_ann_SNP_clean.vcf | awk '{if($7!="PASS")print $0}'
```

## INDEL Variant filteration
```
cd ~/workdir/GATK_tutorial
gatk --java-options "-Xmx2G" VariantFiltration \
-R dog_chr5.fa \
-V raw_variants_ann_SNP.vcf \
--filter-name "indelQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "indelFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 200.0" \
--filter-name "indelSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 10.0" \
--filter-name "indelReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
--filter-name "indelInbreedingCoeff" \
--filter-expression "vc.hasAttribute('InbreedingCoeff') && InbreedingCoeff < -0.8" \
--filter-name "indelDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O raw_variants_ann_INDEL_clean.vcf
```

## How RNA variant calling is different?

https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq

## How somatic variant calling is different?

https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146
