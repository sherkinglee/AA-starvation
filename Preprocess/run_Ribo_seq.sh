#!/bin/bash
#SBATCH -J run_Ribo_seq.sh
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH --ntasks=1

module load bowtie/1.3.0
module load fastqc/0.11.9
module load samtools/1.16.1
module load sratoolkit/2.11.2
module load star/2.7.9a



snakemake -s Ribo-seq-snakemake.py  --cluster "sbatch -p amd-ep2 -N 2 -n 8 -e Ribo.err -o Ribo.out --mem=50G" --jobs 451

