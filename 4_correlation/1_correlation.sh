#!/bin/bash
# Requires: Python 3 (pandas, scipy)

python << 'PYCODE'
import pandas as pd
from scipy.stats import pearsonr

# Input: data.csv in the working directory, one row per sample/condition and one column
# per variable (see data_template.csv for the column layout used in the paper)
df = pd.read_csv("data.csv")
columns = list(df.columns)

corr = df.corr(numeric_only=True)
pvals = pd.DataFrame(index=columns, columns=columns, dtype=float)

for i in columns:
    for j in columns:
        r, p = pearsonr(df[i], df[j])
        pvals.loc[i, j] = p

def star(p):
    if p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return ""

result = pd.DataFrame(index=columns, columns=columns, dtype=object)

for i in columns:
    for j in columns:
        r = corr.loc[i, j]
        p = pvals.loc[i, j]
        result.loc[i, j] = f"{r:.2f}{star(p)}"

print("\n=== Correlation (r + significance) ===\n")
print(result)

corr.to_csv("correlation_r.csv")
pvals.to_csv("correlation_p.csv")
result.to_csv("correlation_matrix.csv")
PYCODE
