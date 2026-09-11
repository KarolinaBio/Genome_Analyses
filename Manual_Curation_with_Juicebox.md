# Manual curation with Juicebox (JBAT) for the Roary Cluster 

Edited from https://github.com/c-zhou/yahs

```
#!/bin/bash
#SBATCH --partition highmem1-sapphirerapids
#SBATCH --qos highmem1
#SBATCH --account acc_group
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --job-name=hic
#SBATCH --output=out_jbat_hic.log

conda activate hic_analysis

juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final_edited.agp JUXXXX_p_purged.fa.fai >out_JBAT.log 2>&1

java -jar -Xmx64G /home/data/group/user/hic_tools/juicer/scripts/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)

```

Note: you need to download a different version of juicer tools into your scripts folder otherwise you'll get a java error

`wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar`


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

