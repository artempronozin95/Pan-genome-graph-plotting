# Pan-genome graph plotting

## Introduction
To date, there are only few programs that can detect the core and variable parts of the pan-genome and build graphs based on them. Those that exist are more suitable for Windows than for Linux. Thus, we have written a script that allows convenient and fast calculation of the number of variable and core genes for pan-genome and pan-transcriptome. The script works based on the results of the Orthofinder program or any other orthologous group search programs. This script can be run on any operating system because it is written in Python. 

## Installation 
You need too install the next Python packages:
+ pandas
+ collections
+ seaborn
+ sys
+ numpy
+ matplotlib

Or you can use our Conda enviroment:
```
conda env create --file env/programs.yaml
conda activate Pangenome
```
## Input
The input data should have the following format (it is best to use Orthofinder results):

| Orthogroup | SRR765127 |	SRR765129 | SRR765130 |	SRR765150 |	SRR765151 |	Total |
|   :---:    | :---:     | :---:      |     :---: | :---:     | :---:     | :---: |
|  OG0000000 |	1	       | 1	        | 3	        | 1	        | 1	        | 6     |
|  OG0000001 |	2	       | 0          | 2	        | 1	        |  0	      | 5     |
|  OG0000002 | 0	       | 1	        | 1	        | 1	        | 2	        | 5     |

## Run
Example of run:
```
python pangenome_plot.py Orthogroups.tsv
```

The script performs 20 iterations of combining samples, to obtain the most accurate number of core and variable genes. 

## Results
As a result, you will get:

+ A table that allows you to estimate what percentage of the total number of samples (organisms) contain orthogroups: 100, 90, 70, 70, 50, 30, 10, and less than 10% of the organisms.

| Percent	| Number of orthogroups |
|   :---:    | :---:     |
| 100% | 1 |
| 90%	| 0 |
| 70% | 3 |
| 50%	| 22 |
| 30%	| 113 |
| 10%	| 8 |
| <10%	| 0 |

+ Files containing the core and variabel part for each iteration.

+ Visualization graph of the core and variable part of the pan-genome or pan-transcriptome. The x-axis shows the samples (organisms), the y-axis shows the number of orthogroups containing the core or variable part.
![image1](https://github.com/artempronozin95/Pan-genome-graph-plotting/blob/main/example/pangenome.png)
