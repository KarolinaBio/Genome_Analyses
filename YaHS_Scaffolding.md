# Standard Scaffolding with YaHS for the Roary Cluster

```
#!/bin/bash
#SBATCH --partition highmem1-sapphirerapids
#SBATCH --qos highmem1
#SBATCH --account acc_group
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --job-name=3778
#SBATCH --output=out_hic.log

conda activate hic_analysis
(juicer pre yahs.out.bin yahs.out_scaffolds_final_edited.agp JU3778_p_purged.fa.fai | sort -k2,2d -k6,6d -T ./ | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

samtools faidx yahs.out_scaffolds_final.fa
cut -f1,2 yahs.out_scaffolds_final.fa.fai > yahs.out_scaffolds_final.chrom.sizes

(java -Xmx32G -jar /home/data/group/user/hic_tools/juicer/scripts/common/juicer_tools.jar pre alignments_sorted.txt out.hic.part yahs.out_scaffolds_final.chrom.sizes) && (mv out.hic.part out.hic)
```

Note: a better version of `juicer_tools` is `juicer_tools.1.9.9_jcuda.0.8.jar` therefore you should change the above script to `/home/data/group/user/hic_tools/juicer/scripts/juicer_tools.1.9.9_jcuda.0.8.jar`
