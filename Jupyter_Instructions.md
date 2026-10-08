# Jupyter Instructions for Insulation Plots
Edited from https://github.com/ercanlab/2024_Aharonoff_et_al/blob/main/scripts/Hi-C/Figure_2.ipynb

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

# Inside the kernel

```
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import bioframe
import cooltools
import cooltools.expected
from cooltools import snipping
import cooltools.lib.plotting
import multiprocess
from matplotlib.colors import LogNorm
from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyBigWig
import csv
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan
from cooltools.insulation import calculate_insulation_score, find_boundaries
```

Load the cool format files into the kernel

```
#bin size
resolution = "5000" 
JU4110_clr = cooler.Cooler(f"/home/data/jfierst/Karolina/JU4110/YaHS/JU4110_test_5000_raw_fixed.cool")
JU4112_clr = cooler.Cooler(f"/home/data/jfierst/Karolina/JU4112/YaHS/JU4112_test_new_5000_raw.cool")
JU4113_clr = cooler.Cooler(f"/home/data/jfierst/Karolina/JU4113/YaHS/JU4113_test_new_5000_raw.cool")
JU4118_clr = cooler.Cooler(f"/home/data/jfierst/Karolina/JU4118/YaHS/JU4118_test_new_5000_raw.cool")
JU4121_clr = cooler.Cooler(f"/home/data/jfierst/Karolina/JU4121/YaHS/JU4121_test_new_5000_raw.cool")

# resolution of your cool file is step size, this is bin size
windows = [50000,100000,150000,200000]

#This part takes a while
JU4110_ins = calculate_insulation_score(JU4110_clr, windows, verbose=True)
JU4112_ins = calculate_insulation_score(JU4112_clr, windows, verbose=True)
JU4113_ins = calculate_insulation_score(JU4113_clr, windows, verbose=True)
JU4118_ins = calculate_insulation_score(JU4118_clr, windows, verbose=True)
JU4121_ins = calculate_insulation_score(JU4121_clr, windows, verbose=True)

```

Once the above step is done, continue in the kernel

```
# home-brew wrapper around smoothing and filling of nan bins, see cooltools for details
def cgi_region(clr,region):
    cg = adaptive_coarsegrain(clr.matrix(balance=True).fetch(region),
                              clr.matrix(balance=False).fetch(region),
                              cutoff=3, max_levels=8)
    cgi = interp_nan(cg)
    return(cg) # chose to not interpolate...the way it fills in bins...makes it seem like stuff is there that isn't

# format x-axis, replace exponentials with intuitive abbreviations
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
```

Next steps, still in the kernel

```
# Feel free to change this region to whatever you want
region = 'X:21,750,000-23,250,000'
start, end = 21_750_000, 23_250_000

# Get the matrix for this region
clr_region = cgi_region(JU4110_clr, region)

# Get insulation values for the same region
ins_region = bioframe.select(JU4110_ins, region)

# Create figure
fig, ax = plt.subplots(figsize=(8, 8))

# Hi-C matrix
im = ax.matshow(
    clr_region,
    cmap='fall',
    norm=LogNorm(vmin=0.00015, vmax=0.01),
    extent=(start, end, end, start),
    aspect='equal'
)
# Show the Hi-C map x- and y-axis ticks
ax.xaxis.set_visible(True)
ax.yaxis.set_visible(True)

# Label the axes
ax.set_xlabel("Chr X position (Mb)")
ax.set_ylabel("Chr X position (Mb)")

# Convert Hi-C map tick labels from bp to Mb
xticks = ax.get_xticks()
ax.set_xticks(xticks)
ax.set_xticklabels(
    [f"{t/1e6:.2f}" for t in xticks]
)
yticks = ax.get_yticks()
ax.set_yticks(yticks)
ax.set_yticklabels(
    [f"{t/1e6:.2f}" for t in yticks]
)
# Colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)

plt.colorbar(im, cax=cax)

# Insulation plot
ax_ins = divider.append_axes(
    "bottom",
    size="30%",
    pad=0.3,
    sharex=ax
)

x = ins_region[['start', 'end']].mean(axis=1)
y = ins_region['log2_insulation_score_150000']

ax_ins.plot(
    x,
    y,
    color='black',
    linewidth=1
)

# Insulation limits - you can set these to your full region for example -2.2,1.8 don't have to constrain it to just -1,1
ax_ins.set_ylim(-1, 1)

# Explicit x limits
ax_ins.set_xlim(start, end)

# Reference line
ax_ins.axhline(
    0,
    color='gray',
    linewidth=0.5,
    linestyle='--'
)
# Formatting
ax.xaxis.set_visible(False)

ax_ins.xaxis.set_visible(True)

ax_ins.set_xlabel("Chromosome X position (Mb)")
ax_ins.set_ylabel("Insulation")

# Convert x-axis labels from bp to Mb
ticks = ax_ins.get_xticks()
ax_ins.set_xticks(ticks)
ax_ins.set_xticklabels(
    [f"{t/1e6:.2f}" for t in ticks],
    rotation=45
)
plt.show()

```

