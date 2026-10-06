# Jupyter Instructions for Insulation Plots

First create a conda environment for Jupyter

```
conda create -n hic_cooltools -c conda-forge python=3.9 numpy=1.26 scipy pandas matplotlib=3.6.3 cython ipykernel pip

conda activate hic_cooltools

conda install -c conda-forge cooler bioframe multiprocess pybigwig

conda install bioconda::cooltools==0.4.0

python -m ipykernel install --user --name hic_cooltools --display-name "Cool Tools Use This One"
```

Then start an interactive Jupyter Notebook session: For example: 2 hours, 4 cores

Select the kernel: Cool Tools Use This One

Make sure you delete your old ~/.local/matplotlib files because they will interfere with your new installation in ~/.conda/

