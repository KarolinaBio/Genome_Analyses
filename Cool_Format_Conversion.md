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
#hicConvertFormat -m out_JBAT.hic --inputFormat hic --outputFormat cool --resolutions 5000 -o JUXXXX.cool #tacks on the resolution

# Step 2 - remove normalization
# See Jupyter instructions below. Then, come back here for step 3

# Step 3 - matrix balancing
cooler balance JUXXXX_raw.cool --max-iters 500 --mad-max 5 --ignore-diags 2 #no output file
```

Jupyter Instruction - Inside the kernel, once all imports have been made

```
input_file = "/home/data/jfierst/Karolina/JU4110/YaHS/JU4110_test_5000.cool"
output_file = "/home/data/jfierst/Karolina/JU4110/YaHS/JU4110_test_5000_raw_fixed.cool"

clr = cooler.Cooler(input_file)

# Keep only the structural bin columns
bins = clr.bins()[:][["chrom", "start", "end"]].copy()

# Keep the raw contact counts
pixels = clr.pixels()[:][["bin1_id", "bin2_id", "count"]].copy()

print("Bins:", len(bins))
print("Pixels:", len(pixels))
print("Total contacts:", pixels["count"].sum())

cooler.create_cooler(
    output_file,
    bins=bins,
    pixels=pixels,
    ordered=True,
    symmetric_upper=True,
    metadata={
        "genome-assembly": "JU4110_yahs.out_scaffolds_final.chrom.sizes",
        "generated-by": "cooler 0.10.2"
    }
)

print("Created:", output_file)
```

