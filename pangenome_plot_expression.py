import pandas as pd
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
from itertools import chain
import warnings 
warnings.filterwarnings('ignore') 
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
    perc_100 = []
    perc_90 = []
    perc_70 = []
    perc_50 = []
    perc_30 = []
    perc_10 = []
    perc_less_10 = []
    core_id = []
    shell_id = []
    cloud_id = []
    x_dict = x.to_dict(orient='list')
    all_values = np.concatenate(list(x_dict.values()))
    for k,v in x_dict.items():
        ort_len = max_len - v.count(0)
        percent_from_max = round((ort_len/max_len)*100)
        if ort_len == max_len:
            core_id.append(k)
            perc_100.append(ort_len)
            full_ort = full_ort + 1
        elif percent_from_max in range(90, 101):
            shell_id.append(k)
            perc_90.append(ort_len)
            precent_90 = precent_90 + 1
        elif percent_from_max in range(70, 90):
            shell_id.append(k)
            perc_70.append(ort_len)
            precent_70 = precent_70 + 1
        elif percent_from_max in range(50, 70):
            shell_id.append(k)
            perc_50.append(ort_len)
            precent_50 = precent_50 + 1
        elif percent_from_max in range(30, 50):
            shell_id.append(k)
            perc_30.append(ort_len)
            precent_30 = precent_30 + 1
        elif percent_from_max in range(10, 30):
            shell_id.append(k)
            perc_10.append(ort_len)
            precent_10 = precent_10 + 1
        elif percent_from_max in range(0, 10):
            cloud_id.append(k)
            perc_less_10.append(ort_len)
            lower_10 = lower_10 + 1
    df = pd.DataFrame(
        {'Percent': ['100%', '90%', '70%' , '50%', '30%', '10%', '<10%'],
         'Number of genes': [full_ort, precent_90, precent_70, precent_50, precent_30, precent_10, lower_10],
         'From': [min(perc_100), min(perc_90), min(perc_70), min(perc_50), min(perc_30), min(perc_10), min(perc_less_10)],
         'To': [max(perc_100), max(perc_90), max(perc_70), max(perc_50), max(perc_30), max(perc_10), max(perc_less_10)]
         })
    df.to_csv(str(path) + '/persent_test.tsv', sep='\t', index=False)
    return df, core_id, shell_id, cloud_id

ortho = pd.read_csv(sys.argv[1], sep='\t')
path = sys.argv[1].rsplit('/', 1)

ortho_T = ortho.T # Transpose dataframe
header=ortho_T.iloc[0]
number_of_samples = len(ortho_T.index) - 1 # Count number of samples 
ortho_T.columns = header
ortho_T = ortho_T[1:]
# Create pan-genome or pan-transcriptome plot
print('Create pan-genome or pan-transcriptome plot')
n = 0
iteration_core = {}
iteration_pan = {}

while n <= 10: # Number of bootstraps
    ortho_T = ortho_T.sample(frac=1, random_state=1) # Randomly shuffle dataframe rows
    first_row = ortho_T.iloc[:1].T # Transpose dataframe
    first_row = first_row.loc[(first_row!=0).any(axis=1)]
    first_row = first_row.reset_index() # Take first column or first row from dataframe
    # Use this first column as core and variable data for the first bootstraps
    core = first_row['Genes'].to_list()
    pangenome = first_row['Genes'].to_list()
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
        if number == 1:
            core_dict[number] = set(core) & set(v['Genes'])
            pangenome = set(core) & set(v['Genes'])
        else:
            core_dict[number] = set(core) & set(v['Genes'])
            pangenome.update(set(core) ^ set(v['Genes']))
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
sns.pointplot(data=df_pangenome, x=df_pangenome["index"], y=df_pangenome['value'], linewidth= 10, color="#4177a1", label='Variable', markers='', errorbar=None)
sns.pointplot(data=df_core, x=df_core["index"], y=df_core['value'], linewidth= 10, color='#b97637', label='Core',  markers='', errorbar=None)
sns.boxplot(data=df_pangenome, x=df_pangenome['index'], y=df_pangenome['value'], width=0.2, hue=df_pangenome['index'], legend=False, palette='dark:#4177a1')
sns.stripplot(data=df_pangenome, x=df_pangenome['index'], y=df_pangenome['value'],
                   jitter=True,
                   marker='o',
                   alpha=0.5,
                   color='black')
sns.boxplot(data=df_core, x=df_core['index'], y=df_core['value'], width=0.2, hue=df_core['index'], legend=False, palette='dark:#b97637')
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

# Staticstic
print('Count statistic')
print('Count percent of samples (species)')
percent, c, s, cl = core_varible(ortho_T, path[0]) # Count percent of samples (species) in orthogroup

core = (percent['Number of genes'].iloc[0]/len(ortho))*100
shell = (percent['Number of genes'].iloc[1:6].sum()/len(ortho))*100
сloud = (percent['Number of genes'].iloc[6].sum()/len(ortho))*100
structure = ['Core', 'Shell', 'Cloud']
count = [round(core), round(shell), round(сloud)]

print('Create a pieplot')
plt.figure(figsize=(20,10))
plt.pie(count, labels=structure, autopct='%1.1f%%',  colors=['#b97637', '#97b7d8', '#fecd7d'], textprops={'fontsize':30})
plt.savefig(str(path[0]) + '/pieplot.png', dpi = 250)

frequencies = ['>10%', '10%', '30%', '50%', '70%', '90%', str(len(ortho_T))]
percent_reversed = percent.iloc[::-1].reset_index(drop=True)
values = percent_reversed['Number of genes'].to_list()
colors = ['#fecd7d', '#97b7d8', '#97b7d8', '#97b7d8', '#97b7d8', '#97b7d8', '#b97637']

data = pd.DataFrame({
    'Frequency': frequencies,
    'Values': values,
    'Colors': colors
})

# Create a bar chart
print('Create a bar plot')
plt.subplots()
plt.figure(figsize=(20,10))
plt.bar(data['Frequency'], data['Values'], color=data['Colors'])
plt.xlabel('Frequency', fontsize=30)
plt.ylabel('Number of genes', fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig(str(path[0]) + '/hist.png', dpi = 250)

ortho_T['Non-zero count'] = ortho_T.apply(lambda row: (row != 0).sum(), axis=1)
ortho_shell = ortho_T[s].apply(lambda row: (row != 0).sum(), axis=1)
ortho_cloud = ortho_T[cl].apply(lambda row: (row != 0).sum(), axis=1)
ortho_T['core_percent'] = (len(c)/ortho_T['Non-zero count'])*100
ortho_T['shell_percent'] = (ortho_shell/ortho_T['Non-zero count'])*100
ortho_T['cloud_percent'] = (ortho_cloud/ortho_T['Non-zero count'])*100
ortho_T.reset_index(inplace=True)
ortho_T_small = ortho_T[['index', 'core_percent', 'shell_percent', 'cloud_percent']]
ortho_T_small['sum'] = ortho_T_small['core_percent'] + ortho_T_small['shell_percent'] + ortho_T_small['cloud_percent'] 
ortho_T_small.to_csv(str(path[0]) +'/core_shell_cloud.tsv', sep='\t')

# Plotting a graph
print('Plotting a proportion graph')
width = 1  # Width of columns

plt.subplots()
plt.figure(figsize=(80,15))

# Bottom part - Core genes
plt.bar(ortho_T_small['index'], ortho_T_small['core_percent'], width, label='Core genes', color='#172540')

# Middle part - Shell genes
plt.bar(ortho_T_small['index'], ortho_T_small['shell_percent'], width, bottom=ortho_T_small['core_percent'], label='Shell genes', color='#97b7d8')

# Top part - Cloud genes
plt.bar(ortho_T_small['index'], ortho_T_small['cloud_percent'], width, bottom=ortho_T_small['core_percent'] + ortho_T_small['shell_percent'], label='Cloud genes', color='#fecd7d')

# Adding a legend
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=50)

# Setting up axis labels and title
plt.ylabel('Proportion (%)', fontsize=50)
plt.yticks(fontsize=50)
plt.xticks(ortho_T_small['index'], rotation=90)
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.3)
plt.savefig(str(path[0]) + '/proportion.png', dpi = 300)


