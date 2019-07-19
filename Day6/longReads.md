## Oxford Nanopore MinION

Oxford Nanopore MinION data from a new bacterial species
```
mkdir -p nano && cd nano
wget https://s3.amazonaws.com/ngs2016/ectocooler_subset.zip
unzip ectocooler_subset.zip
ls ectocooler_subset/
```
You should see a bunch of .fast5 files.
This is only a subset of the reads from the whole run. (Click [here](https://github.com/ljcohen/dib_ONP_MinION/blob/master/Ectocooler/Ectocooler_read_stats_all3runs.ipynb) for stats from the full data set.)


Explore the [hdf5 files](https://support.hdfgroup.org/HDF5/whatishdf5.html)
```
source activate ngs1
conda install -c conda-forge hdf5
h5ls FAST5FILE
h5dump FAST5FILE | less
```

To open a GUI to browse the contents of a FAST5 file type the following, then left­ or right­click on on
groups and datasets to display them:
```
conda install -c eumetsat hdfview
```

[poretools](https://poretools.readthedocs.io/en/latest/)

```
conda install -c bioconda poretools
directory="ectocooler_subset/"
# Convert your .fast5 to .fastq and .fasta files:
poretools fastq $directory > ectocooler_subset.fastq
poretools fasta $directory > ectocooler_subset.fasta
# Check stats
poretools stats -q $directory
```

Install [assembly-stats](https://github.com/sanger-pathogens/assembly-stats):
```
git clone https://github.com/sanger-pathogens/assembly-stats.git
cd assembly-stats/
mkdir build
cd build
cmake ..
make
make test
make install  #sudo make install
assembly-stats/build/assembly-stats ectocooler_subset.fastq
```

Assemble the data with [Canu](https://canu.readthedocs.io/en/latest/)
```
conda install -c bioconda canu  ##workEnv2
```

nanopolish
```
conda install -c bioconda nanopolish
```

