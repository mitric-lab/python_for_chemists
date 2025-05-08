#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv

languages = []
shares = []
with open('pypl_index.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader)  # skip header
    for row in reader:
        languages.append(row[1])
        shares.append(float(row[2].replace('%', '')) / 100)

bar_colors = sns.color_palette('husl', len(languages))[::-1]

fig, ax = plt.subplots(figsize=(8, 5.33))
ax.bar(languages, shares, color=bar_colors)
ax.set_xticks(
    range(len(languages)),
    labels=languages,
    rotation=90,
)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.set_ylabel('Share', fontsize=12)
ax.set_xlim(-0.7, len(languages) - 1.0 + 0.7)

fig.tight_layout()
fig.savefig('../../assets/figures/00-preface/popularity_pypl_202504.svg')

plt.show()

