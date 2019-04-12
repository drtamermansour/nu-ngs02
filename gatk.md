
GATK Variant calling
====================

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
picard_path="/home/ngs-01/miniconda3/envs/ngs1/share/picard-2.19.0-0"

# merge the replicates
java  -Xmx2g -jar $picard_path/picard.jar MergeSamFiles I=BD143_TGACCA_L005.sorted.bam I=BD143_TGACCA_L006.sorted.bam OUTPUT=BD143_TGACCA_merged.sorted.bam
```

Note: Duplicate marking should NOT be applied to amplicon sequencing data or other data types where reads start and stop at the same positions by design.
