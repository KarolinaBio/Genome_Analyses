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
file="JU4118"

bin = str(5000) # resolution of cool file to load in base pairs

# load cool file
ppaclr = cooler.Cooler("".join(["/home/aa7013/hic-to-cool/",file,"/",file,"_30_",bin,"_raw.cool"]))
Load chrom.sizes files

chr=bioframe.read_table('/home/aa7013/hic-to-cool/chrom-sizes_files/pristionchus_pacificus.PRJNA12644.WBPS16.genomic.RS.chrom.sizes',schema='bed3')
chr=bioframe.parse_regions(chr)

print(chr)
```
