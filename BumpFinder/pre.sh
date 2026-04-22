cut -f1 PLUS_ONLY_AWAY_DETAILED.xls | sed 's/_/\t/g' | grep "ENST" > PLUS_ONLY_AWAY_DETAILED.bed
