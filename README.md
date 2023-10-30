# Lecanosticta
Population Genetics of the Brown Spot Needle Blight Fungus

## Processing of nanopore reads
Use duplex mode in dorado to basecall using v10.4.1 superior accuracy model:
```bash
container=path/to/container; pod5=$1; singularity run --nv --app dorado034 $container dorado duplex --device 'cuda:all' dna_r10.4.1_e8.2_5khz_stereo@v1.1 --emit-fastq $pod5 > ${pod5/pod5/sup}_duplex.fastq
```
## Assemble genome using canu with various read length cutoffs:
```bash
#!/bin/bash

#SBATCH --time 96:00:00
#SBATCH --job-name=canu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=CAC48M192_L
#SBATCH --mem=180GB
#SBATCH --mail-type ALL
#SBATCH -A col_farman_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user xxx@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST

assembly=$1
nanoReads=$2
genomeSize=$3
minReadLength=$4

module load ccs/conda/canu-1.9

canu -d ${assembly}_canu_run -p $assembly genomeSize=$genomeSize useGrid=false gridOptionsOVS=" \
 --time 96:00:00 --partition=CAC48M192_L --ntasks=1 --cpus-per-task=16 " minReadLength=$minReadLength stopOnReadQuality=false -nanopore-raw $nanoReads
```
## Assemble genome using flye with corrected reads output from canu:
```bash
#!/bin/bash

#SBATCH --time 72:00:00
#SBATCH --job-name=flye
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --mem=80GB
#SBATCH --mail-type ALL
#SBATCH -A coa_farman_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user xxx@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST

correctedReads=$1
assembly=$2
genomeSize=$3
minOverlap=$4

singularity run --app flye291 /share/singularity/images/ccs/conda/amd-conda9-rocky8.sinf flye  --nano-raw $correctedReads --genome-size $genomeSize  --threads 8 --out-dir ${assembly}_flye --min-overlap $minOverlap
```
