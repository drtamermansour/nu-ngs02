# ChIP-Seq

## Get some sample data
```
conda activate ngs1
cd ~/workdir
wget https://www.ebi.ac.uk/~emily/Online%20courses/NGS/ChIP-seq.zip
unzip ChIP-seq.zip
```
---
## Install prerequisites
For this workflow, we are going to need the following tools:
* Bowtie
* SAMtools
* BEDTools
* UCSC Tools
* MACS2
* deepTools
* bedGraphToBigWig
```
conda install -y bowtie2
conda install -y samtools
conda install -y bedtools
conda install -y deeptools
conda install -y macs2
conda install -y bioconductor-deseq2
conda install -y bioconductor-edger
conda install -y ucsc-bedgraphtobigwig
conda install -y -c r r-essentials
```
---
## Mapping
```
cd ~/workdir/ChIP-seq
bowtie2-build bowtie_index/mm10.fa bowtie_index/mm10
bowtie2 -x bowtie_index/mm10 -U Oct4.fastq -S Oct4.sam -p 1
```
---
## Manipulate SAM output
```
samtools view -bo Oct4.bam Oct4.sam
samtools sort Oct4.bam -o Oct4.sorted.bam
samtools index Oct4.sorted.bam
```
---
## Visualization
### Viewing with tview
```
samtools tview Oct4.sorted.bam bowtie_index/mm10.fa
# or alternatively specify the region you want to see
samtools tview -p chr1:173389928 Oct4.sorted.bam bowtie_index/mm10.fa
```
### Viewing with Online Browsers
```
bamCoverage -b Oct4.sorted.bam --normalizeUsing RPKM -p 5 --extendReads 200 -o oct4.bw
```
---
## Aligning the control sample
```
bowtie2 -x bowtie_index/mm10 -U gfp.fastq -S gfp.sam -p 6
samtools view -bSo gfp.bam gfp.sam
samtools sort gfp.bam -T gfp.temp -o gfp.sorted.bam
samtools index gfp.sorted.bam
bamCoverage -b gfp.sorted.bam --normalizeUsing RPKM -p 5 --extendReads 200 -o gfp.bw
```
---
## macs2 peak calling
MACS2 stands for Model based analysis of ChIP-seq. It was designed for identifying transcription factor binding sites.
MACS2 captures the influence of genome complexity to evaluate the significance of enriched ChIP regions, and improves the 
spatial resolution of binding sites through combining the information of both sequencing tag position and orientation. 
MACS2 can be easily used for ChIP-Seq data alone, or with a control sample to increase specificity.
```
macs2 callpeak -t Oct4.bam -c gfp.bam --format=BAM --name=Oct4 --gsize=138000000 --tsize=26
# cut out only the first 4 columns of the narrowPeak file
cat Oct4_peaks.narrowPeak | cut -f1-4 > Oct_peaks.bed
```
---
## Use bedtools to intersect bed files
It is very common to do genomic region overlap analysis.
For this exercise, we want to find out which promoter regions are bound by Oct4.
To do this we will first need the genomic coordinates of the promoters.
```
curl -LO http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz

# unzip the file
gunzip refGene.txt.gz

# have a look of the file
less -S refGene.txt
```

```
1340    NM_001162506    chr15   +       99074972        99083409        99075354
1521    NM_001162503    chr2    +       122765358       122809551       12278496
901     NM_001271709    chr12   -       41451667        41486057        41452192
103     NM_001357281    chr16   +       31695547        31873356        31696063
1671    NM_001242345    chr6    +       142413430       142419740       14241344
1355    NM_001214911    chr15   -       101050191       101054399       10105246
1387    NM_001309484    chr12   -       105216754       105222737       10521741
1523    NM_001347456    chr5    -       123037126       123047016       12303731
1971    NM_001291205    chr2    -       181688418       181691817       18168910
1969    NM_001291204    chr2    +       181497141       181517962       18149731
1969    NM_001291203    chr2    +       181497141       181517962       18149731
1969    NM_001291202    chr2    +       181497141       181517962       18149731
1969    NM_001291201    chr2    +       181497141       181517962       18149731
1969    NM_001291200    chr2    +       181497141       181517962       18149731
1969    NM_001291197    chr2    +       181497141       181517962       18149731
137     NM_001291196    chr7    -       67231162        67372858        67234867
1096    NM_001166218    chr13   -       67038593        67061170        67040343
1131    NM_001166213    chr8    +       71597645        71608149        71597768
122     NM_001330608    chr7    +       51862014        51994461        51897276
152     NM_001347566    chr13   +       83524742        83667079        83575676
1393    NM_029604       chr3    -       105996956       106001497       10599713
1242    NM_147017       chr2    -       86197802        86198750        86197802
1243    NM_147016       chr2    -       86254764        86255691        86254764
```
<img src="https://angus.readthedocs.io/en/stable/_images/mm10_refgene.2400x2400.jpeg"/>
We will need field 3,4,5,6,13 and rearrange the columns to a bed format.

```
cat refGene.txt | cut -f3-6,13 | awk -v OFS='\t' '{print $1,$3,$4,$5,".",$2}' > mm10_refgene.bed
```
Now, we want to get the promoter regions of the genes. The promoter definition is still arbitrary.
We will use 5kb upstream of the transcription start site (TSS) as a promoter region.
```
# download the genome info file
curl -LO http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# get the 5kb promoter region
bedtools flank -i mm10_refgene.bed -g mm10.chrom.sizes -l 5000 -r 0 -s > mm10.genes.5kb.promoters.bed

# have a look of the file
head mm10.genes.5kb.promoters.bed

# compare to the original file
head mm10_refgene.bed
```

```
chr15	99069972	99074972	Troap	.	+
chr2	122760358	122765358	Sqor	.	+
chr12	41486057	41491057	Lrrn3	.	-
chr16	31690547	31695547	Dlg1	.	+
chr6	142408430	142413430	Spx	.	+
chr15	101054399	101059399	Fignl2	.	-
chr12	105222737	105227737	Tcl1	.	-
chr5	123047016	123052016	Morn3	.	-
chr2	181691817	181696817	Rgs19	.	-
chr2	181492141	181497141	Tpd52l2	.	+
```
Peaks overlap with the promoter regions
```
bedtools intersect -a Oct_peaks.bed -b mm10.genes.5kb.promoters.bed -wa -wb > Oct4_peaks_overlap_promoter.txt
# select genes
cat Oct4_peaks_overlap_promoter.txt | cut -f8 | sort | uniq > genes_with_Oct4_binding.txt
```
## Alternatively you can use HOMER
```
conda install homer
perl $CONDA_PREFIX/share/homer-4.9.1-6/configureHomer.pl -install mm10
annotatePeaks.pl Oct4_summits.bed mm10 > genes_with_Oct4_binding_homer.txt
less -S genes_with_Oct4_binding_homer.txt
```
