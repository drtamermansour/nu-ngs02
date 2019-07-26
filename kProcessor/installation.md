## Install kProcessor
```
conda activate ngs1
mkdir ~/workdir/kProcessor && cd ~/workdir/kProcessor

python --version

# if your version is 3.6.x
wget https://github.com/drtamermansour/nu-ngs02/raw/master/kProcessor/dist/kProcessor-0.1-cp36-cp36m-linux_x86_64.whl
pip install kProcessor-0.3-cp36-cp36m-linux_x86_64.whl

# else if your version is 3.7
wget https://github.com/drtamermansour/nu-ngs02/raw/master/kProcessor/dist/kProcessor-0.1-cp37-cp37m-linux_x86_64.whl
pip install kProcessor-0.3-cp37-cp37m-linux_x86_64.whl
```
---
## Install Jupyter
```
conda install jupyter
```
---
## Download kProcessor jupyter notebook
```
wget https://raw.githubusercontent.com/drtamermansour/nu-ngs02/master/kProcessor/Tutorial.ipynb
```
---
## Start Jupyter
```
jupyter notebook
```
