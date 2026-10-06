# Converting to Cool format for insulation plots

Edited from https://github.com/ercanlab/2024_Aharonoff_et_al/blob/main/scripts/Hi-C/hicTocool.sh

Install the environment:

`conda create -n hic_explorer -c conda-forge -c bioconda python=3.10 hicexplorer=3.7.6`


```
#!/bin/bash
#SBATCH --job-name=hicConvert
#SBATCH --account=acc_group
#SBATCH --qos=highmem1
#SBATCH --partition=highmem1-sapphirerapids
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=hicConvert.log

conda activate hic_explorer

# hic to cool
#hicConvertFormat -m out_JBAT.hic --inputFormat hic --outputFormat cool --resolutions 5000 -o JUXXXX.cool #tacks on the resolution

# remove normalization
hicConvertFormat -m JUXXXX_5000.cool --inputFormat cool --outputFormat cool --load_raw_values -o JUXXXX_raw.cool

# matrix balancing
cooler balance JUXXXX_raw.cool --max-iters 500 --mad-max 5 --ignore-diags 2 #no output file
```
