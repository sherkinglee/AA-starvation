bedtools intersect -v -a PLUS_peaks.txt -b MINUS_peaks.txt > PLUS_only_peaks.txt
bedtools closest -a PLUS_only_peaks.txt -b MINUS_peaks.txt -d > PLUS_ONLY_CLOSEST.txt
bedtools intersect -v -a MINUS_peaks.txt -b PLUS_peaks.txt > MINUS_only_peaks.txt
bedtools closest -a MINUS_only_peaks.txt -b PLUS_peaks.txt  -d > MINUS_ONLY_CLOSEST.txt
