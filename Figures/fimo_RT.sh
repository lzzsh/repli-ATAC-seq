# select symbol from txt
cat JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt| grep MOTIF | cut -d " " -f 2,3 > symbol.txt

# classification
awk 'NR>1 {print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10}' OFS='\t' fimo.txt > fimo_temp.txt
mv fimo_temp.txt fimo.txt
sed 's/\t*$//' fimo.txt > fimo_temp.txt
mv fimo_temp.txt fimo.txt
