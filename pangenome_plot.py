import pandas as pd
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

# This function form table that count percent of samples (species) in orthogroup 
def core_varible(x, path): 
    max_len = len(x)
    full_ort = 0
    precent_90 = 0
    precent_70 = 0
    precent_50 = 0
    precent_30 = 0
    precent_10 = 0
    lower_10 = 0
    for w in x.columns:
        ort_len = len(x[x[w] != 0])
        percent_from_max = round((ort_len/max_len)*100)
        if ort_len == max_len:
            full_ort = full_ort + 1
        if percent_from_max in range(90, 99):
            precent_90 = precent_90 + 1
        if percent_from_max in range(70, 89):
            precent_70 = precent_70 + 1
        if percent_from_max in range(50, 69):
            precent_50 = precent_50 + 1
        if percent_from_max in range(30, 49):
            precent_30 = precent_30 + 1
        if percent_from_max in range(10, 29):
            precent_10 = precent_10 + 1
        if percent_from_max in range(0, 9):
            lower_10 = lower_10 + 1

    df = pd.DataFrame(
        {'Percent': ['100%', '90%', '70%' , '50%', '30%', '10%', '<10%'],
         'Number of orthogroups': [full_ort, precent_90, precent_70, precent_50, precent_30, precent_10, lower_10]
         })
    df.to_csv(str(path) + '/persent.tsv', sep='\t', index=False)

ortho = pd.read_csv(sys.argv[1], sep='\t')

path = sys.argv[1].rsplit('/', 1)

ortho_T = ortho.drop('Total', axis=1).T # Transpose dataframe
header=ortho_T.iloc[0]
number_of_samples = len(ortho_T.index) - 1 # Count number of samples 
ortho_T.columns = header
ortho_T = ortho_T[1:]
core_varible(ortho_T, path[0]) # Count percent of samples (species) in orthogroup
n = 0
iteration_core = {}
iteration_pan = {}

while n <= 20: # Number of bootstraps
    ortho_T = ortho_T.sample(frac=1, random_state=1) # Randomly shuffle dataframe rows
    first_row = ortho_T.iloc[:1].T # Transpose dataframe
    first_row = first_row.loc[(first_row!=0).any(axis=1)]
    first_row = first_row.reset_index() # Take first column or first row from dataframe
    # Use this first column as core and variable data for the first bootstraps
    core = first_row['Orthogroup'].to_list()
    pangenome = first_row['Orthogroup'].to_list()
    groups = ortho_T.groupby(ortho_T.index, sort=False)
    n_gr = groups.ngroups # Remove duplicated groups

    core_dict = defaultdict(list)
    pangenome_dict = defaultdict(list)
    number = 1
    # Searching for core and variable genes
    for k, v in groups:
        v = v.T
        v = v.loc[(v!=0).any(axis=1)]
        v = v.reset_index()
        for w in v['Orthogroup']:
            if w in core:
                core_dict[number].append(w)
            else:
                pangenome.append(w)
        pangenome_dict[number] = set(pangenome)
        core = core_dict[number] # Update of core set for next iteration
        number = number + 1


    core_number = {k:len(v) for k, v in core_dict.items()} # The number of core genes in every sample
    pangenome_number = {k:len(v) for k, v in pangenome_dict.items()} # The number of variable genes in every sample

    iteration_core[n] = core_number
    iteration_pan[n] = pangenome_number
    n = n + 1


df_core= pd.DataFrame.from_dict(iteration_core)
df_core = df_core.reset_index()
df_core = pd.melt(df_core, id_vars=['index'])
df_pangenome= pd.DataFrame.from_dict(iteration_pan)
df_pangenome = df_pangenome.reset_index()
df_pangenome = pd.melt(df_pangenome, id_vars=['index'])
df_pangenome['percent'] = abs((df_pangenome['value'] / max(df_pangenome['value']))* 100) # Percent of variable genes
df_core['percent'] = abs((df_core['value'] / max(df_core['value']))* 100) # Percent of core genes

df_pangenome.to_csv(str(path[0]) + '/variable.tsv', sep='\t', index=False) # File with variable genes at every interaction
df_core.to_csv(str(path[0]) + '/core.tsv', sep='\t', index=False) # File with core genes at every interaction

# Building a plot of core and variable genes 
plt.figure(figsize=(50,15))
sns.set(style="white", color_codes=True)
sns.pointplot(data=df_pangenome, x=df_pangenome["index"], y=df_pangenome['value'], linewidth=5, color="#720026", label='Variable')
sns.pointplot(data=df_core, x=df_core["index"], y=df_core['value'], linewidth=5, color='#00725F', label='Core')
sns.boxplot(data=df_pangenome, x=df_pangenome['index'], y=df_pangenome['value'], width=0.2, hue=df_pangenome['index'], legend=False, palette="rocket")
sns.stripplot(data=df_pangenome, x=df_pangenome['index'], y=df_pangenome['value'],
                   jitter=True,
                   marker='o',
                   alpha=0.5,
                   color='black')
sns.boxplot(data=df_core, x=df_core['index'], y=df_core['value'], width=0.2, hue=df_core['index'], legend=False, palette= "crest")
sns.stripplot(data=df_core, x=df_core['index'], y=df_core['value'],
                   jitter=True,
                   marker='o',
                   alpha=0.5,
                   color='black')

plt.legend(fontsize=30) 
# Choosing range of xticks
if number_of_samples in range(1, 50):
    step = 1
if number_of_samples in range(51, 100):
    step = 3
if number_of_samples in range(101, 500):
    step = 5
if number_of_samples in range(501, 1000):
    step = 7

plt.xticks(np.arange(0, number_of_samples, step=step), fontsize=30, rotation=90)
plt.yticks(fontsize=30)

plt.xlabel("Genome Number",fontsize=30)
plt.ylabel("Gene Cluster Number", fontsize=30)

plt.savefig(str(path[0]) + '/pangenome.png', dpi = 250)


