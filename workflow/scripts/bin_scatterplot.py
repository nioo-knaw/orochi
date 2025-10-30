import argparse
import pandas as pd
import plotly.express as px

parser = argparse.ArgumentParser(description="Read input file for bin scatterplot")
parser.add_argument("-input", type=str, required=True, help="Path to the input TSV file")
parser.add_argument("-output", type=str, required=True, help="Path to save the output HTML plot")
args = parser.parse_args()

# Read tab-separated file
df = pd.read_csv(args.input, sep=",")

# Create interactive scatter plot
fig = px.scatter(
    df,
    x="completeness",
    y="contamination",
    color="phylum",
    hover_data=df.columns,
    title="Bin Quality (Completeness vs Contamination)",
    labels={
        "completeness": "Completeness (%)",
        "contamination": "Contamination (%)"
    },
    color_discrete_sequence=px.colors.qualitative.Set2
)

# --- Add black borders around points ---
fig.update_traces(
    marker=dict(
        line=dict(width=0.8, color='black'),
        opacity=0.85,
        size=10
    )
)

# --- Clean scientific styling ---
fig.update_layout(
    template="simple_white",
    plot_bgcolor="white",
    paper_bgcolor="white",
    width=800,
    height=600,
    legend_title_text="Phylum",
    title_font=dict(size=20, family="Arial", color="black"),
    font=dict(size=14, family="Arial", color="black"),
    xaxis=dict(
        range=[0, 100],
        showline=True, linecolor="black",
        showgrid=True, gridcolor="lightgrey", zeroline=False
    ),
    yaxis=dict(
        range=[0, 100],
        showline=True, linecolor="black",
        showgrid=True, gridcolor="lightgrey", zeroline=False
    )
)

# Save as interactive HTML file
fig.write_html(args.output)
print(f"✅ Interactive plot saved to {args.output}")
