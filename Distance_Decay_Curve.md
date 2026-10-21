# Jupyter Instructions
Edited from https://github.com/ercanlab/2024_Aharonoff_et_al/blob/main/scripts/Hi-C/Figures_3_S8.ipynb

```
# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
import cooler
import bioframe
import cooltools
import cooltools.expected
from cooltools import snipping
import cooltools.lib.plotting
# Set up parallelization
import multiprocess
from matplotlib.colors import LogNorm
from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyBigWig
import csv
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan
# Import python package for working with cooler files and tools for analysis
from cooltools.insulation import calculate_insulation_score, find_boundaries
```

Go to the terminal and do this: ` awk 'BEGIN {OFS="\t"; print "chrom", "start", "end", "name"} {print $1, 0, $2, $1":0-"$2}' yahs.out_scaffolds_final.chrom.sizes > JU4118_cool.chrom.sizes ` Then, come back to Jupyter

```
resolution = "5000" 

# load cool file
JU4118_clr = cooler.Cooler(f"/home/data/jfierst/Karolina/JU4118/YaHS/JU4118_test_new_5000_raw.cool")

# Load chrom.sizes files
df = pd.read_csv('/home/data/jfierst/Karolina/JU4118/YaHS/JU4118_cool.chrom.sizes', sep="\t")
view = bioframe.make_viewframe(df) # view = bioframe.parse_regions(df)
print(view)


```
Define functions to measure P(s) and its derivative (taken from cooltools)

```
# a home-made wrapper around scaling distance decay curve and its derivative. See cooltools for more details
# returns 2 outputs, P(s) and its derivative

def scaling_agg(input_clr,input_regions):
    cvd = cooltools.expected.diagsum(
        clr=input_clr,
        regions=input_regions,
        transforms={'balanced': lambda p: p['count']*p['weight1']*p['weight2']}
    )
    lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(cvd) #bins_per_order_magnitude=20 to get more dots
    lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(
        lb_cvd,
        binned_exp_slope=lb_slopes
    )
    lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * input_clr.binsize
    lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * input_clr.binsize
    return(lb_cvd_agg,lb_slopes_agg)
```

Define parse_regions

```
def parse_regions_compat(regions, chromsizes):
    if isinstance(regions, pd.DataFrame):
        result = regions.copy()
    else:
        result = pd.DataFrame(
            [bioframe.parse_region(r, chromsizes) for r in regions],
            columns=['chrom', 'start', 'end']
        )

    if 'name' not in result.columns:
        result['name'] = result['chrom'].astype(str)

    return result[['chrom', 'start', 'end', 'name']]

bioframe.parse_regions = parse_regions_compat
```

Measure P(s) and its derivative
