# [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html)
## Install

```
source activate ngs1
conda install ensembl-vep
```
## Download Example VCF
```
mkdir ~/workdir/VEP && cd ~/workdir/VEP
wget https://wsend.net/a9fb8555134f004dba3030100c090925/ensembl-vep.zip
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
# You'll be asked if you want to create the cache folder if it doesn't exist, yes
# The following species/files are available; which do you want (can specify multiple separated by spaces or 0 for all):
# choose 82 for drosophila_melanogaster
```
## Run against cached files
```
vep -i examples/drosophila_melanogaster_BDGP6.vcf --cache --species drosophila_melanogaster -o vep.txt
less vep.txt
```
## View results in browser
```
firefox vep.txt_summary.html
```
