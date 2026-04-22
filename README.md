# Amino acid starvation
Translational Resistance to Amino Acid Starvation in Breast Cancer Cells Revealed by Ribosome Profiling

Fajin Li<sup>1,2†\*</sup>, Siqiong Zhang<sup>1,2†</sup>, Haoran Duan<sup>1,2</sup>,Dieter A. Wolf<sup>1,2</sup>,and Xuerui Yang<sup>3</sup>

<sup>1</sup> Westlake Laboratory of Life Sciences and Biomedicine, Hangzhou, Zhejiang, China. 
<sup>2</sup> School of Medicine, Westlake University, Hangzhou, Zhejiang, China.
<sup>3</sup> School of Life Sciences, Tsinghua University, Beijing, China.
<sup>†</sup> These authors contributed equally to this work
<sup>*</sup> Corresponding author. Email: lfj17[at]tsinghua[dot]org[dot]cn



---

## **Introduction**

This file is a description of how the results presented in the manuscript were generated, including codes of pre-processing and codes for generating the figures. All scripts used for processing BAM files should be excuted in Linux platform and other scripts can be both used in Linux and windows platform.

## **Dependencies**
The softwares and their versions are listed  in the **RiboMiner_requirements.txt** file.

## **Codes**

+ **Preprocess**: The codes for preprocessing of all Ribo-seq data.
```
$ bash run_Ribo_seq.sh
```
+ **Codes_for_regenerations.Rmd**: The codes for formal analysis, which were used for generating all figures in the manuscript.
+ **BumpFinder**: The codes for bump finder algorithm (for example: GSE142822). Original code for this part of analysis can be found at https://github.com/apataskar/bump_finder_example2
  ```
 	## 01. counting for Ribo-seq data
	$ bash bam2bed.sh

	## 02. merge the counts file
	$ R
	setwd("D:\\Projects\\02. AA starvation\\GSE142822\\13.metagenePlot\\BumpFinder")
	Trp_24h_1 <- read.table("./SRR10814089.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Trp_24h_2 <- read.table("./SRR10814090.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Trp_48h_1 <- read.table("./SRR10814093.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Trp_48h_2 <- read.table("./SRR10814094.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Trp_24h_full <- read.table("./SRR10814092.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Trp_48h_full <- read.table("./SRR10814096.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Tyr_48h_1 <- read.table("./SRR10814099.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Tyr_48h_2 <- read.table("./SRR10814100.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Tyr_48h_full_1 <- read.table("./SRR10814097.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	Tyr_48h_full_2 <- read.table("./SRR10814098.Aligned.toTranscriptome.out.sorted.bed_counts.txt",sep = "\t",stringsAsFactors = F)
	
	colnames(Trp_24h_1) <- c("id","start","end","Trp_24_1")
	colnames(Trp_24h_2) <- c("id","start","end","Trp_24_2")
	colnames(Trp_48h_1) <- c("id","start","end","Trp_48h_1")
	colnames(Trp_48h_2) <- c("id","start","end","Trp_48h_2")
	colnames(Trp_24h_full) <- c("id","start","end","Trp_24h_full")
	colnames(Trp_48h_full) <- c("id","start","end","Trp_48h_full")
	colnames(Tyr_48h_1) <- c("id","start","end","Tyr_48h_1")
	colnames(Tyr_48h_2) <- c("id","start","end","Tyr_48h_2")
	colnames(Tyr_48h_full_1) <- c("id","start","end","Tyr_48h_full_1")
	colnames(Tyr_48h_full_2) <- c("id","start","end","Tyr_48h_full_2")
	
	data_Tyr_48h <- Reduce(cbind,list(Tyr_48h_full_1,Tyr_48h_full_2[,4],Tyr_48h_1[,4],Tyr_48h_2[,4]))
	colnames(data_Tyr_48h) <- c(colnames(Tyr_48h_full_1),"Tyr_48h_full_2","Tyr_48h_1","Tyr_48h_2")
	
	data_Trp_48h <- Reduce(cbind,list(Trp_48h_full,Trp_48h_full[,4],Trp_48h_1[,4],Trp_48h_2[,4]))
	colnames(data_Trp_48h) <- c("id","start","end","Trp_48h_full_1","Trp_48h_full_2","Trp_48h_1","Trp_48h_2")
	
	data_Trp_24h <- Reduce(cbind,list(Trp_24h_full,Trp_24h_full[,4],Trp_24h_1[,4],Trp_24h_2[,4]))
	colnames(data_Trp_24h) <- c("id","start","end","Trp_24h_full_1","Trp_24h_full_2","Trp_48h_1","Trp_48h_2")
	
	write.table(data_Tyr_48h,"data_Tyr_48h.txt",quote = F,col.names = T,sep = "\t",row.names = F)
	write.table(data_Trp_48h,"data_Trp_48h.txt",quote = F,col.names = T,sep = "\t",row.names = F)
	write.table(data_Trp_24h,"data_Trp_24h.txt",quote = F,col.names = T,sep = "\t",row.names = F)

	## 03. get idx for each data set
	$ cut -f 1 data_Tyr_48h.txt > Tyr_48_idx.txt 
	$ cut -f 1 data_Trp_24h.txt > Trp_24_idx.txt 
	$ cut -f 1 data_Trp_48h.txt > Trp_48_idx.txt
	

	## 04. Profiler
	$ Rscript profiler.R
	setwd("D:\\Projects\\02. AA starvation\\GSE142822\\13.metagenePlot\\BumpFinder\\Tyr_48h/")
	list= read.table ("ids_Tyr_48h.txt", head=F)
	
	data_Trp_48h = read.table ("data_Tyr_48h.txt", head=T)
	
	for (gene2 in 1: nrow (list))
	{
	
	gene=  data.frame (id=list[gene2,"V1"])
	gene = merge (gene, data_Trp_48h, by="id")
	gene = gene[ order (gene$start),]
	
	
	
	gene$minus = rowMeans(gene[,c(4,5)])
	gene$plus = rowMeans(gene[,c(6,7)])
	
	gene$average = rowMeans (gene[,c("minus","plus")])
	
	if ( log2(sum(gene$average)) > 5)
	{
	gene$plus_log= ave(gene$plus, FUN=function(x) c(0, diff(x)))
	gene$minus_log= ave(gene$minus, FUN=function(x) c(0, diff(x)))
	
	
	gene$plus_log_norm= (gene$plus_log-min(gene$plus_log))/(max(gene$plus_log)-min(gene$plus_log))
	gene$minus_log_norm= (gene$minus_log-min(gene$minus_log))/(max(gene$minus_log)-min(gene$minus_log))
	
	gene$prof= gene$plus_log_norm -  gene$minus_log_norm
	gene$prof_norm= (gene$prof-min(gene$prof))/(max(gene$prof)-min(gene$prof))
	
	data=t(gene$prof)
	data= data.frame (data)
	data$name = unique (gene$id)
	data[,c(which(colnames(data)=="name"),which(colnames(data)!="name"))]
	
	data2=t(gene$prof_norm)
	data2= data.frame (data2)
	data2$name = unique (gene$id)
	data2[,c(which(colnames(data2)=="name"),which(colnames(data2)!="name"))]
	
	data3=t(gene$minus_log_norm)
	data3= data.frame (data3)
	data3$name = unique (gene$id)
	data3[,c(which(colnames(data3)=="name"),which(colnames(data3)!="name"))]
	
	data4=t(gene$plus_log_norm)
	data4= data.frame (data4)
	data4$name = unique (gene$id)
	data4[,c(which(colnames(data4)=="name"),which(colnames(data4)!="name"))]
	
	write.table (gene, "GENE_DETAILS_Tyr_48h.txt", sep="\t", quote=F, row.names=F, col.names=F, append=T)
	write.table (data, "PROF_Tyr_48h_RC_NONSCALED_DETAILS.txt",sep="\t", quote=F, row.names=F, col.names=F,append=T)
	write.table (data2, "PROF_Tyr_48h_RC_SCALED_DETAILS.txt",sep="\t", quote=F, row.names=F, col.names=F,append=T)
	write.table (data3, "PROF_Tyr_48h_MINUS_LOG_NORM.txt",sep="\t", quote=F, row.names=F, col.names=F,append=T)
	write.table (data4, "PROF_Tyr_48h_PLUS_LOG_NORM.txt",sep="\t", quote=F, row.names=F, col.names=F,append=T)
	

	## 05. script_1.R
	$ Rscript script_1.R

	library ("pracma")
	gene = read.table ("GENE_DETAILS_Tyr_48h.txt", head=F)
	colnames (gene)= c("id","start","end","m1","m2","p1","p2","minus","plus","average","plus_log","minus_log","plus_log_norm","minus_log_norm","prof","prof_norm")
	gene=gene[complete.cases (gene),]
	minus_log_peaks=findpeaks(as.numeric(gene$minus_log), nups=1, ndowns=1, minpeakheight=10 )
	plus_log_peaks=findpeaks(as.numeric(gene$plus_log), nups=1, ndowns=1, minpeakheight=10 )
	gene_plus_trial= gene[ plus_log_peaks[,2],]
	gene_minus_trial= gene[ minus_log_peaks[,2],]
	gene_minus_trial_inverse= gene[ -minus_log_peaks[,2],]
	gene_plus_trial_inverse= gene[ -plus_log_peaks[,2],]
	gene_minus_trial$class_minus = "MINUS_PEAK"
	gene_minus_trial_inverse$class_minus = "MINUS_NOPEAK"
	gene_minuss_details = rbind (gene_minus_trial, gene_minus_trial_inverse)
	gene_plus_trial_inverse$class_plus = "PLUS_NOPEAK"
	gene_plus_trial$class_plus = "PLUS_PEAK"
	gene_plus_details = rbind (gene_plus_trial, gene_plus_trial_inverse)
	write.table (gene_minus_trial[,c(1,2,3,8,9,11,12)], "MINUS_peaks.txt",sep="\t", quote=F, row.names=F, col.names=F)
	write.table (gene_plus_trial[,c(1,2,3,8,9,11,12)], "PLUS_peaks.txt",sep="\t", quote=F, row.names=F, col.names=F)
	gene$diff= gene$plus_log - gene$minus_log


	## 06. run bedtools
	$ bash bedtools.sh

	module load bedtools/2.30.0
	bedtools intersect -v -a PLUS_peaks.txt -b MINUS_peaks.txt > PLUS_only_peaks.txt
	bedtools closest -a PLUS_only_peaks.txt -b MINUS_peaks.txt -d > PLUS_ONLY_CLOSEST.txt
	bedtools intersect -v -a MINUS_peaks.txt -b PLUS_peaks.txt > MINUS_only_peaks.txt
	bedtools closest -a MINUS_only_peaks.txt -b PLUS_peaks.txt  -d > MINUS_ONLY_CLOSEST.txt

	## 07. run script_2.R
	$ Rscript script_2.R

	plus_only = read.table ("PLUS_ONLY_CLOSEST.txt", head=F)
	plus_only_away= plus_only [ plus_only$V15 < 0 | plus_only$V15> 9,]
	gene$id_chr = paste (gene$id,gene$start, gene$end, sep="_")
	plus_only_away$id_chr = paste (plus_only_away$V1, plus_only_away$V2, plus_only_away$V3, sep="_")
	plus_only_away_detailed=merge (plus_only_away, gene, by="id_chr")
	write.table (plus_only_away_detailed, "PLUS_ONLY_AWAY_DETAILED.xls",sep="\t", quote=F, row.names=F)
	pdf ("PLUS_MINUS_PROF_CUTOFF_REASON.pdf")
	par (mfrow=c(3,1))
	plot (density ( plus_only_away_detailed$minus_log_norm), type="l", col="red", main="MINUS", lwd=3)
	plot (density ( plus_only_away_detailed$plus_log_norm), type="l", col="green", main="PLUS", lwd=3)
	plot (density ( plus_only_away_detailed$plus_log_norm), type="l", col="green", lwd=3)
	lines (density ( plus_only_away_detailed$minus_log_norm), type="l", col="red", main="MINUS", lwd=3)
	dev.off()
	plus_only_detailed_prof10= plus_only_away_detailed [ plus_only_away_detailed$plus_log_norm - plus_only_away_detailed$minus_log_norm > 0.1, ]
	plus_only_detailed_prof10_diff = plus_only_detailed_prof10[ plus_only_detailed_prof10$diff > 20,]
	write.table (plus_only_detailed_prof10_diff, "PLUS_ONLY_DETAILED_PROF10_DIFF.xls",sep="\t", quote=F, row.names=F)

	#PART THREE: RUN bedtools.sh before this part. MINUS_SPECIFIC
	minus_only = read.table ("MINUS_ONLY_CLOSEST.txt", head=F)
	minus_only_away= minus_only [ minus_only$V15 < 0 | minus_only$V15> 9,]
	minus_only_away$id_chr = paste (minus_only_away$V1, minus_only_away$V2, minus_only_away$V3, sep="_")
	minus_only_away_detailed=merge (minus_only_away, gene, by="id_chr")
	write.table (minus_only_away_detailed, "MINUS_ONLY_AWAY_DETAILED.xls",sep="\t", quote=F, row.names=F)
	minus_only_detailed_prof10= minus_only_away_detailed [ minus_only_away_detailed$plus_log_norm - minus_only_away_detailed$minus_log_norm > 0.1, ]
	#NO diff cutoff- for stats)
	write.table (minus_only_detailed_prof10, "MINUS_ONLY_DETAILED_PROF10.xls",sep="\t", quote=F, row.names=F)
	minus_only_detailed_prof10_diff = minus_only_detailed_prof10[ minus_only_detailed_prof10$diff > 20,]
	write.table (minus_only_detailed_prof10_diff, "MINUS_ONLY_DETAILED_PROF10_DIFF.xls",sep="\t", quote=F, row.names=F)

	## 08. run pre.sh
	cut -f1 PLUS_ONLY_AWAY_DETAILED.xls | sed 's/_/\t/g' | grep "ENST" > PLUS_ONLY_AWAY_DETAILED.bed
	cut -f1 MINUS_ONLY_AWAY_DETAILED.xls | sed 's/_/\t/g' | grep "ENST" > PLUS_ONLY_AWAY_DETAILED.bed

	## 09. run script_3.R
	$ Rscript script_3.R

	bumps = read.table ("PLUS_ONLY_AWAY_DETAILED.bed", head=F)
	bumps$mid = round (bumps$V2 + ((bumps$V3-bumps$V2)/2) )
	cds = read.table ("CDS.bed", head=F)
	colnames (cds) = c("V1","CS","CE")
	bumps_cds=  merge (bumps, cds, by="V1")
	bumps_cds$aa_mid = round ((bumps_cds$mid - bumps_cds$CS)/3) +1
	bumps_cds_atleast30= bumps_cds[ bumps_cds$aa_mid >= 30,]
	bumps_cds$aa_length= round ((bumps_cds$CE - bumps_cds$CS)/3) 
	bumps_cds_atleast30= bumps_cds[ bumps_cds$aa_mid >= 30 & ( (bumps_cds$aa_length - bumps_cds$aa_mid) > 30),]
	write.table (bumps_cds_atleast30, "bumps_cds_atleast30.txt",sep="\t", quote=F, row.names=F)

	## 10. preprocess.sh
	$ bash preprocess.sh
	
	## 11. Plot
	$ Rscript plot.R

	## 12. get sequences around the bump
	$ python bumps_30_each_2fasta.py

	## 13. sequence logo with the seq2logo
  ```
+  **Metagene analysis and codon specific pausing**: the scripts for metagene analysis and codon-specific pausing analysis are available in RiboMiner (https://github.com/xryanglab/RiboMiner).
+  **CodonMeta**: Scripts for codon-level metagene analysis.
