To setup the `tauso` package locally we first need to install conda / mamba environment. 


```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
```

After setting up mamba / conda, you may install the environment with 

```
mamba create -n tauso -f environment.yml
```

Then you may run `pip install -e .` to install locally. 

To also run the tests / notebooks, you may further install the development libraries

```
mamba install -n tauso -f environment-dev.yml
```

To run the tests or train the models, we need to assign canonical genes to the original OligoAI and create an index:

```
python -m notebooks.data.OligoAI.assign_canonical_gene
python -m notebooks.utils.data
```
