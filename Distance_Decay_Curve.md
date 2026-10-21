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
view = bioframe.make_viewframe(df)
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

Define parse_regions because it doesn't exist

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

Makes cases match

```
chrom_lookup = {
    chrom.lower(): chrom
    for chrom in JU4118_clr.chromnames
}

view['chrom'] = view['chrom'].map(
    lambda chrom: chrom_lookup.get(chrom.lower(), chrom)
)

print(view['chrom'].head())
```

Measure P(s) and its derivative

```
X = scaling_agg(JU4118_clr,view.iloc[[0]].reset_index(drop=True))
V = scaling_agg(JU4118_clr,view.iloc[[1]].reset_index(drop=True))
IV = scaling_agg(JU4118_clr,view.iloc[[2]].reset_index(drop=True))
I = scaling_agg(JU4118_clr,view.iloc[[3]].reset_index(drop=True))
II = scaling_agg(JU4118_clr,view.iloc[[4]].reset_index(drop=True))
III = scaling_agg(JU4118_clr,view.iloc[[5]].reset_index(drop=True))
```

Determine local maxima of the derivative for each chromosome (used to annotate on the plot).

```
slope_I = I[1]
slope_I = slope_I.iloc[:25]
max_slope_I = slope_I[slope_I['slope'] == slope_I['slope'].max()]
loop_size_I = max_slope_I['s_bp'].values[0]
slope_II = II[1]
slope_II = slope_II.iloc[:23]
max_slope_II = slope_II[slope_II['slope'] == slope_II['slope'].max()]
loop_size_II = max_slope_II['s_bp'].values[0]
slope_III = III[1]
slope_III = slope_III.iloc[:25]
max_slope_III = slope_III[slope_III['slope'] == slope_III['slope'].max()]
loop_size_III = max_slope_III['s_bp'].values[0]
slope_IV = IV[1]
slope_IV = slope_IV.iloc[:25]
max_slope_IV = slope_IV[slope_IV['slope'] == slope_IV['slope'].max()]
loop_size_IV = max_slope_IV['s_bp'].values[0]
slope_V = V[1]
slope_V = slope_V.iloc[:24]
max_slope_V = slope_V[slope_V['slope'] == slope_V['slope'].max()]
loop_size_V = max_slope_V['s_bp'].values[0]
slope_X = X[1]
slope_X = slope_X.iloc[:24]
max_slope_X = slope_X[slope_X['slope'] == slope_X['slope'].max()]
loop_size_X = max_slope_X['s_bp'].values[0]

loop_size_I, loop_size_II, loop_size_III, loop_size_IV, loop_size_V, loop_size_X
```

Plot the P(s) and its derivative as a function of separation in bp
```
f, axs = plt.subplots(
    figsize=(8,10),
    ncols=1,
    nrows=2,
    sharex=True, sharey='row',
    gridspec_kw={'height_ratios':[6,3]}
)

lw=1.25
a=1
lw2=2.5
############################################ plot Ps ################################################
ax0 = axs[0]
ax0.loglog(
    I[0]['s_bp'],
    I[0]['balanced.avg'],
    'o-', markersize=0,
    label='I',
    alpha=a, color='#CCCCCC',linewidth=1
) #### plot P(s,I) on log-log space
ax0 = axs[0]
ax0.loglog(
    II[0]['s_bp'],
    II[0]['balanced.avg'],
    'o-', markersize=0,
    label='II',
    alpha=a, color='#AAAAAA',linewidth=lw
) #### plot P(s,III) on log-log space
ax0 = axs[0]
ax0.loglog(
    III[0]['s_bp'],
    III[0]['balanced.avg'],
    'o-', markersize=0,
    label='III',
    alpha=a, color='#888888',linewidth=lw
) #### plot P(s,V) on log-log space
ax0 = axs[0]
ax0.loglog(
    IV[0]['s_bp'],
    IV[0]['balanced.avg'],
    'o-', markersize=0,
    label='IV',
    alpha=a, color='#666666',linewidth=lw
) #### plot P(s,N) on log-log space
ax0 = axs[0]
ax0.loglog(
    V[0]['s_bp'],
    V[0]['balanced.avg'],
    'o-', markersize=0,
    label='V',
    alpha=a, color='#444444',linewidth=lw
) #### plot P(s,N) on log-log space
ax0.loglog(
    X[0]['s_bp'],
    X[0]['balanced.avg'],
    'o-', markersize=0,
    label='X',
    alpha=1.0, color='#D55E00',linewidth=lw
) #### plot P(s,X) on the same axis, ax0

############################################ plot Ps derivative #######################################
ax1 = axs[1]
ax1.semilogx(
    I[1]['s_bp'],
    I[1]['slope'],
    'o-', markersize=0,
    label='I',
    alpha=a,color='#CCCCCC',linewidth=lw
) #### plot log-derivative of P(s,I) on semi-log space
ax1.semilogx(
    II[1]['s_bp'],
    II[1]['slope'],
    'o-', markersize=0,
    label='II',
    alpha=a,color='#AAAAAA',linewidth=lw
) #### plot log-derivative of P(s,III) on semi-log space
ax1.semilogx(
    III[1]['s_bp'],
    III[1]['slope'],
    'o-', markersize=0,
    label='III',
    alpha=a,color='#888888',linewidth=lw
) #### plot log-derivative of P(s,V) on semi-log space
ax1.semilogx(
    IV[1]['s_bp'],
    IV[1]['slope'],
    'o-', markersize=0,
    label='IV',
    alpha=a,color='#666666',linewidth=lw
) #### plot log-derivative of P(s,N) on semi-log space
ax1.semilogx(
    V[1]['s_bp'],
    V[1]['slope'],
    'o-', markersize=0,
    label='V',
    alpha=a,color='#444444',linewidth=lw
) #### plot log-derivative of P(s,N) on semi-log space
ax1.semilogx(
    X[1]['s_bp'],
    X[1]['slope'],
    'o-', markersize=0,
    label='X',
    alpha=a,color='#D55E00',linewidth=lw
) #### plot log-derivative of P(s,X) on the same axis, ax1

######################################## do some labeling #############################################
# putting grids
ax0.grid(lw=0.5, axis = 'x')
ax1.grid(lw=0.5, axis = 'x')
#ax2.grid(lw=0.5)

# label y-axis
ax0.set_ylabel('P(s)', fontsize=20, labelpad=15)
ax1.set_ylabel('slopes', fontsize=20, labelpad=15)
ax0.tick_params(axis='y', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)
#ax2.set_ylabel('P(s,X)/P(s,A)')

# format x-axis, replacing exponentials with intuitive abbreviations
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
        
format_ticks(ax=ax1,x=True,y=False,rotate=False) # x-axis is shared, only one needs to be formatted
ax1.set_xlabel('separation, s', fontsize=20, labelpad=15)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)


# add legend
ax0.legend(fontsize=14)
plt.axvline(x = loop_size_I, color = '#CCCCCC', ymin = 0.7, ymax = 0.75, linewidth = lw2)
plt.axvline(x = loop_size_II, color = '#AAAAAA', ymin = 0.6, ymax = 0.65, linewidth = lw2)
plt.axvline(x = loop_size_III, color = '#888888', ymin = 0.3, ymax = 0.35, linewidth = lw2)
plt.axvline(x = loop_size_IV, color = '#666666', ymin = 0.2, ymax = 0.25, linewidth = lw2)
plt.axvline(x = loop_size_V, color = '#444444', ymin = 0.1, ymax = 0.15, linewidth = lw2)
plt.axvline(x = loop_size_X, color = '#D55E00', ymin = 0.3, ymax = 0.35, linewidth = lw2)
plt.tight_layout()
```

Done! But make sure you save your file! 

```
plt.savefig("U4118_distance_decay.pdf"]), format='pdf', dpi=300)
```
