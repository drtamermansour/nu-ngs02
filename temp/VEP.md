# [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html)
## Install

```
source activate ngs1
conda install ensembl-vep
## VEP uses perl. The current version of perl in our conda env is 5.26.2. However, for some reason conda install all perl dependencies for 5.22.0
## Therefor, vep will fail to run because to search for the perl libraries at: /home/ngs/miniconda3/envs/ngs1/lib/site_perl/5.26.2
## We can change this behaviour by setting the environmental variable PERL5LIB to the correct path
## This solution is based on suggestions here: https://stackoverflow.com/questions/58290190/how-to-fix-perl-from-anaconda-not-installing-bioperl-bailing-out-the-installat
## and here: https://github.com/simroux/VirSorter/#note-for-conda-installation
export PERL5LIB=/home/ngs/miniconda3/envs/ngs1/lib/perl5/site_perl/5.22.0/
```

## Download Example VCF
```
mkdir -p ~/workdir/VEP && cd ~/workdir/VEP
wget https://data.cyverse.org/dav-anon/iplant/home/drtamermansour/DIBSI/ensembl-vep.zip
unzip ensembl-vep.zip
```
VEP can be used to run against online databases and offline caches.
Using cache is faster than online databases, however cache is very large to download.

## Run VEP

```
vep -i examples/homo_sapiens_GRCh38.vcf --database
less variant_effect_output.txt
```

## View results in browser
```
firefox variant_effect_output.txt_summary.html
```

## To install cache files
```
vep_install --NO_HTSLIB
# The current default conda installation of VEP is old version and thus this command will suggest updating the VEP API. We will ignore an continue.
# It will ask if you want to install cache file, you should say "yes"
# You'll be asked if you want to create the cache folder if it doesn't exist, yes
# The following species/files are available; which do you want (can specify multiple separated by spaces or 0 for all):
# choose 27 for drosophila_melanogaster
# We do not need the FASTA files or more plugins
```
## Run against cached files
```
vep -i examples/drosophila_melanogaster_BDGP6.vcf --cache --species drosophila_melanogaster -o dm_vep.txt
less dm_vep.txt
```
## View results in browser
```
firefox dm_vep.txt_summary.html
```

For reference and details, check for the [VEP tutorial](http://uswest.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html) 



