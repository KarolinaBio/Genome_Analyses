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

Also, instead of `<(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')` try simply doing `yahs.out_scaffolds_final.chrom.sizes` because then you're actually able to see the chromosomes instead of just `assembly` which really makes no sense
