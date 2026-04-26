#!/bin/bash
#SBATCH --mem 48G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

cd ~/sharedscratch/software/interproscan

module load python/3.9.12
module load openjdk/11.0.20.1_1

IPS_VERSION=5.77-108.0

# tar -pxvzf interproscan-${IPS_VERSION}-64-bit.tar.gz
	# uncomment if needed to download
cd interproscan-${IPS_VERSION}

## check
ls -lh interproscan.sh
test -x interproscan.sh && echo "interproscan.sh is executable"

# Initialise/index models
python3 setup.py -f interproscan.properties


# check interproscan
./interproscan.sh
./interproscan.sh -version

# run official test jobs
./interproscan.sh -i test_all_appl.fasta -f tsv -dp
./interproscan.sh -i test_all_appl.fasta -f tsv

# validate
ls -lh test_all_appl.fasta*.tsv
head test_all_appl.fasta.tsv
wc -l test_all_appl.fasta.tsv
cut -f1-8 test_all_appl.fasta.tsv | head

## FINAL CHECKLIST
uname -m
perl -version
python3 --version
java -version
echo $JAVA_HOME
test -x interproscan.sh && echo "Executable OK"
./interproscan.sh -version
ls -lh test_all_appl.fasta*.tsv

