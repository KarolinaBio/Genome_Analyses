# Converting to Cool format for insulation plots

Edited from https://github.com/ercanlab/2024_Aharonoff_et_al/blob/main/scripts/Hi-C/hicTocool.sh

Install the environment:

`conda create -n hic_explorer -c conda-forge -c bioconda python=3.10 hicexplorer=3.7.6`


```
#!/bin/bash
#SBATCH --job-name=hicConvert
#SBATCH --account=acc_jfierst
#SBATCH --qos=highmem1
#SBATCH --partition=highmem1-sapphirerapids
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=hicConvert.log

conda activate hic_explorer

# Step 1 - hic to cool
#hicConvertFormat -m out_JBAT.hic --inputFormat hic --outputFormat cool --resolutions 5000 -o JU4113.cool #tacks on the resolution

# Step 2 - remove normalization
hicConvertFormat -m JU4113_5000.cool --inputFormat cool --outputFormat cool -o JU4113_5000_raw.cool --load_raw_values
conda deactivate

conda activate cool_env

# Step 3 - matrix balancing
cooler balance JU4113_5000_raw.cool --max-iters 500 --mad-max 5 --ignore-diags 2 #no output file
```
