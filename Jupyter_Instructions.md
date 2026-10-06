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
#For reference, these re the species names, don't paste these into the kernel
#JU4110 = "C. sp. 61"
#JU4112 = "C. sp. 62"
#JU4113 = "C. sp. 63"
#JU4118 = "C. sp. 65"
#JU4121 = "C. sp. 66"
#JU3778 = "C. sp. 59"
#JU3779 = "C. sp. 58"

resolution = "5000" #bin size
JU4110_clr = cooler.Cooler(f"/home/data/jfierst/Karolina/JU4110/YaHS/JU4110_5000_new_raw.cool")

# resolution of your cool file is step size, this is bin size
windows = [50000,100000,150000,200000]

#This part takes a while
JU4110_ins = calculate_insulation_score(JU4110_clr, windows, verbose=True)

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

f, axs = plt.subplots(
    figsize=(25, 10),
    ncols=7,
    nrows=1
)
```

Next steps, still in the kernel

```
JU4110_norm = LogNorm(vmin=0.00015,vmax=0.01)

#Insulation parameter
ins_min,ins_max = [-1.0,1.0]
region = 'X:6,500,000-8,000,000'
start, end =6_500_000, 8_000_000
extents = (start, end, end, start)

#Encountered an error, trying to remake the cool files
ax = axs[6]
clr_region = cgi_region(JU4110_clr,region)
im = ax.matshow(
    clr_region,
    cmap='fall',
    norm=JU4110_norm,
    extent=extents
)


divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1) # color axis for color bar
plt.colorbar(im, cax=cax)

#Insulation
ax_ins = divider.append_axes("bottom", size="30%", pad=0.2, sharex=ax) # axis for insulation score
ins_region = bioframe.select(oonins, region)
ax_ins.plot(ins_region[['start', 'end']].mean(axis=1), 
            ins_region['log2_insulation_score_150000']) # where you input the window size you used

# force min/max values
ax_ins.set_ylim([ins_min,ins_max])

#Formatting
ax.xaxis.set_visible(False) # hide axis labels
ax_ins.xaxis.set_visible(True) # hide axis labels
format_ticks(ax,x=True,y=True,rotate=True) # format y-axis for matrix
plt.xticks(rotation=45)
```

