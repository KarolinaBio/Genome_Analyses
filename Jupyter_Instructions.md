# Jupyter Instructions for Insulation Plots

First create a conda environment for Jupyter

```
conda create -n hic_cooltools -c conda-forge python=3.9 numpy=1.26 scipy pandas matplotlib cython ipykernel pip

conda activate hic_cooltools

conda install -c conda-forge cooler bioframe multiprocess pybigwig

conda install bioconda::cooltools==0.4.0

python -m ipykernel install --user --name hic_cooltools --display-name "Hi-C (cooltools 0.4.0)"
```

Then start an interactive Jupyter Notebook session: For example: 2 hours, 4 cores

Select the kernel: Hi-C (cooltools 0.4.0)
