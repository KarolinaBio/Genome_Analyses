# Jupyter Instructions

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

Next

```
resolution = "5000" 

# load cool file
JU4118_clr = cooler.Cooler(f"/home/data/jfierst/Karolina/JU4118/YaHS/JU4118_test_new_5000_raw.cool"

# Load chrom.sizes files

chr=bioframe.read_table('/home/data/jfierst/Karolina/JU4118/YaHS/yahs.out_scaffolds_final.chrom.sizes',schema='bed3')
chr=bioframe.make_viewframe(chr)

print(chr)
```
