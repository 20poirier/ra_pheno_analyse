import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 22})

import numpy as np


pd.set_option('expand_frame_repr', False)
df = pd.read_table('data/assoc1.assoc',delim_whitespace=True)
# print(df.head(10))

df['p_adj'] = -np.log10(df.P)
df.chr = df.CHR.astype('category')

# print(df.head(10))

df['ind'] = range(len(df))
df_grouped = df.groupby(('CHR'))

# df has as many lines as there are SNPs
num_snps = df.shape[0]
sig_thresh = 0.05/num_snps
print("Significance threshold = ", sig_thresh)

# Significant SNPs
sig_SNPs = df.SNP[df.P < sig_thresh]

# Number of significant SNPs
print(len(sig_SNPs))
print(sig_SNPs)
# print(df_grouped.head(10))

fig = plt.figure(figsize = (10, 6))
ax = fig.add_subplot(111)
colors = ['#E24E42','#008F95']
x_labels = []
x_labels_pos = []

for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='p_adj',color=colors[num % len(colors)], ax=ax, s=5)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))


# Plot Banefolini threshold
ax.plot([0, max(df.BP)], [-np.log10(sig_thresh), -np.log10(sig_thresh)], lw=1, color='black')
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df)])
ax.set_ylim([0,None])
ax.set_xlabel('Chromosome')
ax.set_ylabel(r'$-\log_{10}(P)$')
plt.xticks(fontsize = 22,rotation=60)
plt.yticks(fontsize = 22)
ax.set_title('Manhattan plot', fontsize = 28)



plt.show()
fig.savefig('data.png')