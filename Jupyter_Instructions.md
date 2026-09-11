Jupyter Instructions for Insulation Plots

conda create -n cool_tools -c conda-forge -c bioconda python=3.9 numpy=1.26.4 pandas=1.5.3 scipy=1.10.1 matplotlib=3.7.5 cython cooler bioframe multiprocess pybigwig jupyter jpykernel

conda activate cool_tools
pip install cooltools==0.4.0

python -m ipykernel install --user --name cool_tools --display-name "Hi-C (cooltools 0.4.0)"

Then start an interactive Jupyter Notebook session 2 hours, 4 cores
Select the kernel: Hi-C (cooltools 0.4.0)
