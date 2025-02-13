# pEA_Matrix_Generation
Generating pEA/sumEA matrix for any ML algorithm 



## Installation
1. git clone https://github.com/LichtargeLab/BigPipeline.git
2. conda env create -f ./BigPipeline/environment.yml
3. conda activate pyBigPipeline


## Usage
Required arguments:
| Argument                | Descripion |
| ---------------------- |--------------------- |
| --VCF                | Path to annotated VCF file |
| --samples            |Path to two-column CSV file with sample IDs in the first column and patient labels (cases=1, controls=0) in the second column. There should be no header row in the csv|
| --savepath           | Path for output files |
| --cores              | number of cpus to use |

Optional arguments:
| Argument                 | Descripion |
| ---------------------- |--------------------- |
| --maxaf  | sets maximum allele frequency threshold for ML pipelines (default: 0.01) |
| --minaf      |sets minimum allele frequency threshold for ML pipelines (default: 0)|
| --method           | pEA or sumEA matrix calculation. (options: pEA, sumEA)|
| --chrX       | Whether to include (1) or exclude (0) sex chromosomes in analysis (options: 1, 0 / default: 1 )|
