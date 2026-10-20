#!/bin/bash
# Requires: Python 3 (pandas, numpy, scipy, matplotlib)

python << 'PYCODE'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# =========================
# 1. Input data
# =========================
# Input: data.csv in the working directory, one row per sample/condition and one column
# per variable (see data_template.csv for the column layout used in the paper)
df = pd.read_csv("data.csv")
columns = list(df.columns)

# =========================
# 2. Correlation and p-values
# =========================
corr = df.corr(numeric_only=True)

pvals = pd.DataFrame(index=columns, columns=columns, dtype=float)
for i in columns:
    for j in columns:
        r, p = pearsonr(df[i], df[j])
        pvals.loc[i, j] = p

def stars(p):
    if p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return ""

# save tables
corr.to_csv("correlation_r.csv")
pvals.to_csv("correlation_p.csv")

display_tbl = pd.DataFrame(index=columns, columns=columns, dtype=object)
for i in columns:
    for j in columns:
        display_tbl.loc[i, j] = f"{corr.loc[i,j]:.2f}{stars(pvals.loc[i,j])}"
display_tbl.to_csv("correlation_matrix_with_stars.csv")

# =========================
# 3. Plot settings
# =========================
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["font.size"] = 9

n = len(columns)
fig, ax = plt.subplots(figsize=(14, 12), dpi=300)

# heatmap
im = ax.imshow(corr.values, cmap="coolwarm", vmin=-1, vmax=1)

# ticks
ax.set_xticks(np.arange(n))
ax.set_yticks(np.arange(n))
ax.set_xticklabels(columns, rotation=60, ha="right")
ax.set_yticklabels(columns)

# white background and clean frame
ax.set_facecolor("white")
for spine in ax.spines.values():
    spine.set_visible(False)

# minor gridlines
ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
ax.grid(which="minor", color="white", linestyle="-", linewidth=0.8)
ax.tick_params(which="minor", bottom=False, left=False)

# =========================
# 4. Annotate lower triangle only
# =========================
for i in range(n):
    for j in range(n):
        if i >= j:
            r = corr.iloc[i, j]
            p = pvals.iloc[i, j]
            label = f"{r:.2f}{stars(p)}"
            text_color = "black" if abs(r) < 0.65 else "white"
            ax.text(j, i, label, ha="center", va="center", color=text_color, fontsize=7)
        else:
            # blank upper triangle for cleaner publication style
            ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, color="white", ec="white"))

# title
ax.set_title("Correlation matrix of reactor performance and gene metrics", pad=18, fontsize=12)

# colorbar
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Pearson correlation coefficient (r)", rotation=90)

# footnote
fig.text(
    0.02, 0.015,
    "* p < 0.05, ** p < 0.01; correlation across reactors (n = 4), interpreted as exploratory.",
    ha="left", va="bottom", fontsize=8
)

plt.tight_layout(rect=[0, 0.03, 1, 1])

# save
plt.savefig("correlation_heatmap_WR_style.png", bbox_inches="tight", dpi=600)
plt.savefig("correlation_heatmap_WR_style.pdf", bbox_inches="tight")
plt.close()

print("Done.")
print("Saved files:")
print(" - correlation_r.csv")
print(" - correlation_p.csv")
print(" - correlation_matrix_with_stars.csv")
print(" - correlation_heatmap_WR_style.png")
print(" - correlation_heatmap_WR_style.pdf")
PYCODE
