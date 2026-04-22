#!/bin/bash

#SBATCH -J bam2bed.sh
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH -p intel-sc3
#SBATCH --ntasks=1
#SBATCH --mem=50G

module load fastqc/0.11.9
module load star/2.7.9a
module load sratoolkit/2.11.2
module load samtools/1.16.1
module load bcftools/1.14
module load bowtie/2.4.2
module load bedtools/2.30.0

for i in `ls|grep SRR`;do
       bedtools bamtobed -i $i/${i%_*}.Aligned.toTranscriptome.out.sorted.bam > $i/${i%_*}.Aligned.toTranscriptome.out.sorted.bed
done

for i in `ls|grep SRR`;do
        less -S $i/${i%_*}.Aligned.toTranscriptome.out.sorted.bed|awk -F "\t" 'BEGIN{OFS="\t"}{print $1,$2,$3}' > $i/${i%_*}.Aligned.toTranscriptome.out.sorted.final.bed
        bedtools coverage -a /home/dawolfLab/dongxingjian/lifj/Reference/human/ensemble109/longest_transcript_select_100_windows.bed -b $i/${i%_*}.Aligned.toTranscriptome.out.sorted.final.bed | cut -f1,2,3,4 > $i/${i%_*}.Aligned.toTranscriptome.out.sorted.bed_counts.txt
done
