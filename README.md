# Pan-genome graph plotting

## Introduction
To date, there are only few programs that can detect the core and variant parts of the pan-genome and build graphs based on them. Those that exist are more suitable for Windows than for Linux. Thus, we have written a script that allows convenient and fast calculation of the number of variable and core genes for pan-genome and pan-transcriptome. The script works based on the results of the Orthofinder program or any other orthologous group search programs. This script can be run on any operating system because it is written in Python. 

## Installation 

## Input
The input data should have the following format (it is best to use Orthofinder results):
```
| Orthogroup | SRR765127 |	SRR765129 | SRR765130 |	SRR765150 |	SRR765151 |	Total |
|   :---:    | :---:     | :---:      |     :---: | :---:     | :---:     | :---: |
|  OG0000000 |	1	       | 1	        | 3	        | 1	        | 1	        | 6     |
|  OG0000001 |	2	       | 0          | 2	        | 1	        |  0	      | 5     |
|  OG0000002 | 0	       | 1	        | 1	        | 1	        | 2	        | 5     |

```
