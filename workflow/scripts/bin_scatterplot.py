import argparse
import pandas as pd
import matplotlib.pyplot as plt

# set up argument parser
parser = argparse.ArgumentParser(description="Read input file for bin scatterplot")
parser.add_argument(
    "-input", 
    type=str, 
    required=True, 
    help="Path to the input file (TSV) resulting from Orochi 'binning' module execution"
)
parser.add_argument(
    "-output",
    type=str,
    required=True,
    help="Path to save the output scatterplot"
)
args = parser.parse_args()

df = pd.read_csv(args.input, sep=",")
fig = plt.subplots(figsize=(8,6))

plt.figure(figsize=(8,6))
scatter = plt.scatter(
    df["completeness"],
    df["contamination"],
    c=pd.factorize(df["phylum"])[0],  # color
    cmap="tab20",
    alpha=0.8,
    edgecolors="k"
)

plt.xlabel("Completeness (%)")
plt.ylabel("Contamination (%)")
plt.title("Bin Quality (Completeness vs Contamination)")

# Legend: one color per taxonomy group
handles, labels = scatter.legend_elements(prop="colors")
plt.legend(handles, df["phylum"].unique(), title="Phylum", bbox_to_anchor=(1.05, 1))

plt.tight_layout()
plt.savefig(args.output, dpi=500)
