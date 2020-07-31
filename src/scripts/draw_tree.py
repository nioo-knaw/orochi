# https://github.com/hardingnj/snakemake_pres/blob/master/deploy_cluster/scripts/draw_tree.py
# Fix couldn't connect to display "localhost:0.0" error
# http://python.omics.wiki/plot/matplotlib/tkinter-tclerror
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

# Fix for the error: qt.qpa.xcb: could not connect to display localhost:15.0
import os
os.environ['QT_QPA_PLATFORM']='offscreen'

dist = pd.read_table(snakemake.input[0], index_col=0, header=None, comment="#")
# refer to variables defined outside
SAMPLES = list(snakemake.params.samples.keys())
dist.columns = SAMPLES
dist.index = SAMPLES

d = squareform(np.array(dist))
z = linkage(d, method="complete")

r = dendrogram(
    z, labels=SAMPLES, count_sort=True,
    color_threshold=0,
    above_threshold_color='k')
plt.savefig(snakemake.output[0])
