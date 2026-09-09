# Genome_Analyses

Manual curation with Juicebox (JBAT) for the Roary Cluster


```
#!/bin/bash
#SBATCH --partition highmem1-sapphirerapids
#SBATCH --qos highmem1
#SBATCH --account acc_group
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --job-name=3778
#SBATCH --output=out_jbat.log

conda activate hic_analysis

juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final_edited.agp JU3778_p_purged.fa.fai >out_JBAT.log 2>&1

conda deactivate
conda activate java_new

(java -jar -Xmx32G /home/data/group/user/hic_tools/juicer/scripts/common/juicer_tools.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)

```

You can replace `<(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}'))` with `<(echo "assembly 73472353"))` found in the `out_JBAT.log` file.
