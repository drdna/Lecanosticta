# Lecanosticta
Population Genetics of the Brown Spot Needle Blight Fungus

## Processing of nanopore reads
Use duplex mode in dorado to basecall using r10.4.1 superior accuracy model:
```bash
#!/bin/bash

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user xxx@uky.edu            #Where to send email
#SBATCH --account=gol_xxx_uksr              #Name of account to run under
#SBATCH --partition=V4V16_SKY32M192_L            #Partition

#SBATCH --job-name=dorado-dup_sup_gpu    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=28G                # total memory per node (4 GB per cpu-core is default)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --time=48:00:00          # total run time limit (HH:MM:SS)

pod5=$1
container=/share/singularity/images/ccs/conda/lcc-conda-8-rocky8.sinf

module load ccs/singularity-3.8.2

singularity run --nv --app dorado034 $container dorado duplex --device 'cuda:all' dna_r10.4.1_e8.2_400bps_sup@v4.2.0 --emit-fastq $pod5 > ${pod5/pod5/sup}_duplex.fastq
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
#SBATCH -A col_xxx_uksr
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
#SBATCH -A coa_xxx_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user xxx@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST

correctedReads=$1
assembly=$2
genomeSize=$3
minOverlap=$4

singularity run --app flye291 /share/singularity/images/ccs/conda/amd-conda9-rocky8.sinf flye  --nano-raw $correctedReads --genome-size $genomeSize  --threads 8 --out-dir ${assembly}_flye --min-overlap $minOverlap
```
