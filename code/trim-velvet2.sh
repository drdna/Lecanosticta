#!/bin/bash

#SBATCH --time 48:00:00
#SBATCH --job-name=trim-velvet
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=normal
#SBATCH --mem=500GB
#SBATCH --mail-type ALL
#SBATCH -A coa_farman_uksr		# change account ID as needed
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu	# change email to include those who should receive notifcations about submitted/started/completed/failed jobs (separate with commas, no spaces)

echo "SLURM_NODELIST: "$SLURM_NODELIST

## ARGUMENTS TO SCRIPT

# specify directory containing fastq data
dir=$1

# specify read IDs (only the prefixes that come at the start of the filenames: e.g. StrainX-1_1.fastq.gz would be StrainX-1
readID=$2

# specify starting k-mer value (must be odd number)
lowK=$3

# starting k-mer value (must be odd number)
highK=$4

# step size (must be even number)
step=$5


## SYSTEM COMMANDS

# Create directory for assembly
mkdir $f 

# copy forward reads into directory
cp $dir/$f*_1*f*q* $f/

# copy reverse reads into directory
cp $dir/$f*_2*f*q* $f/

# Change into assembly directory
cd $f

# Run trimmomatic (should always be set to "yes" unless input reads are already trimmed, or a prior assembly attempt failed after the trimming step)
if [ $trim == 'yes' ]
then
  singularity run --app trimmomatic039 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf trimmomatic PE \
  -threads 16 -phred33 -trimlog ${f}_errorlog.txt \
  $f*_1*.f*q* $f*_2*.f*q* \
  ${f}_R1_paired.fq ${f}_R1_unpaired.fq \
  ${f}_R2_paired.fq ${f}_R2_unpaired.fq \
  ILLUMINACLIP:/project/farman_uksr/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:20:20 MINLEN:90;
fi

# create an interleaved dataset from the paired reads resulting from trimmomatic

module load python-2.7.18-gcc-9.3.0-5efgwu4
python /project/farman_uksr/BASH_SCRIPTS/interleave-fastq.py *1_paired.fq *2_paired.fq > interleaved.fq
rm  *_*paired.fq
module unload python-2.7.18-gcc-9.3.0-5efgwu4


# run the assembly in singularity
singularity run --app perlvelvetoptimiser226 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf VelvetOptimiser.pl \
 -s $lowK -e $highK -x $step -d velvet_assembly -f ' -shortPaired -interleaved -fastq interleaved.fq'

# change the name of the assembly to the provided strain ID
mv velvet_assembly/contigs.fa velvet_assembly/${f}".fasta"

# change the name of the logfile to include the strain ID
prefix=`ls velvet_assembly/*logfile.txt'
mv $prefix /scratch/farman/ASSEMBLIES/${f}_${prefix/*\//}


# change back into the base directory 
cd ..

#rm -r $f



