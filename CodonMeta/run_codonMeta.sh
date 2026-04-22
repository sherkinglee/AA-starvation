#!/bin/bash
#SBATCH -J Meta.sh
#SBATCH -o Meta.out
#SBATCH -e Meta.err
#SBATCH -N 1

workdir=`pwd`
BamDir=$workdir/../07.STAR
Ref=/home/dawolfLab/dongxingjian/lifj/Reference/human/ensemble109
transcript=/home/dawolfLab/dongxingjian/lifj/Reference/human/ensemble109/longest_cds_sequences.fa
results=$workdir/MA
attribute=$workdir/configure3_A.txt
trans_info=$Ref/longest.transcripts.info.txt
groups='108Tcells-Trp-full-24h,108Tcells-Trp-depletion-24h,108Tcells-Trp-full-48h,108Tcells-Trp-depletion-48h,MD55A3cells-Tyr-full,MD55A3cells-Tyr-depletion-48h'
replicates='108Tcells-Trp-full-24h__108Tcells-Trp-depletion-24h-1,108Tcells-Trp-depletion-24h-2__108Tcells-Trp-full-48h__108Tcells-Trp-depletion-48h-1,108Tcells-Trp-depletion-48h-2__MD55A3cells-Tyr-full-1,MD55A3cells-Tyr-full-2__MD55A3cells-Tyr-depletion-48h-1,MD55A3cells-Tyr-depletion-48h-2'


mkdir -p codonMeta/codon codonMeta/AA


while read codon
do
        python CodonMetaAnalysis.py -f $attribute -c $trans_info -o codonMeta/codon/$codon -M RPKM -l 100 -n 1 --type1 $codon -F $transcript

done < codons.txt

while read AA
do
        python CodonMetaAnalysis.py -f $attribute -c $trans_info -o codonMeta/AA/$AA -M RPKM -l 100 -n 1 --type2 $AA -F $transcript

done < AAs.txt


python PlotCodonMeta.py -i codonMeta/codon/TGG_motifDensity_dataframe.txt -o codonMeta/codon/TGG_Trp_48h -g 108Tcells-Trp-full-48h,108Tcells-Trp-depletion-48h -r 108Tcells-Trp-full-48h__108Tcells-Trp-depletion-48h-1,108Tcells-Trp-depletion-48h-2 --mode mean --slide-window y
python PlotCodonMeta.py -i codonMeta/AA/W_motifDensity_dataframe.txt -o codonMeta/AA/W_Trp_48h -g 108Tcells-Trp-full-48h,108Tcells-Trp-depletion-48h -r 108Tcells-Trp-full-48h__108Tcells-Trp-depletion-48h-1,108Tcells-Trp-depletion-48h-2 --mode mean --slide-window y
python PlotCodonMeta.py -i codonMeta/codon/TGG_motifDensity_dataframe.txt -o codonMeta/codon/TGG_Trp_24h -g 108Tcells-Trp-full-24h,108Tcells-Trp-depletion-24h -r 108Tcells-Trp-full-24h__108Tcells-Trp-depletion-24h-1,108Tcells-Trp-depletion-24h-2 --mode mean --slide-window y
python PlotCodonMeta.py -i codonMeta/AA/W_motifDensity_dataframe.txt -o codonMeta/AA/W_Trp_24h -g 108Tcells-Trp-full-24h,108Tcells-Trp-depletion-24h -r 108Tcells-Trp-full-24h__108Tcells-Trp-depletion-24h-1,108Tcells-Trp-depletion-24h-2 --mode mean --slide-window y

python PlotCodonMeta.py -i codonMeta/codon/TAC_motifDensity_dataframe.txt -o codonMeta/codon/TAC_Tyr_48h -g MD55A3cells-Tyr-full,MD55A3cells-Tyr-depletion-48h -r MD55A3cells-Tyr-full-1,MD55A3cells-Tyr-full-2__MD55A3cells-Tyr-depletion-48h-1,MD55A3cells-Tyr-depletion-48h-2 --slide-window y --mode mean
python PlotCodonMeta.py -i codonMeta/codon/TAT_motifDensity_dataframe.txt -o codonMeta/codon/TAT_Tyr_48h -g MD55A3cells-Tyr-full,MD55A3cells-Tyr-depletion-48h -r MD55A3cells-Tyr-full-1,MD55A3cells-Tyr-full-2__MD55A3cells-Tyr-depletion-48h-1,MD55A3cells-Tyr-depletion-48h-2 --slide-window y --mode mean
python PlotCodonMeta.py -i codonMeta/AA/Y_motifDensity_dataframe.txt -o codonMeta/AA/Y_Tyr_48h -g MD55A3cells-Tyr-full,MD55A3cells-Tyr-depletion-48h -r MD55A3cells-Tyr-full-1,MD55A3cells-Tyr-full-2__MD55A3cells-Tyr-depletion-48h-1,MD55A3cells-Tyr-depletion-48h-2 --slide-window y --mode mean
