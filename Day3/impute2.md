# [Impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)


## Download

```bash
wget -c https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz
tar -zxvf impute_v2.3.2_x86_64_static.tgz
```


## Run Example

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
