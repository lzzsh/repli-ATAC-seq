awk 'NR>1 {print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10}' OFS='\t' fimo.txt > fimo_temp.txt
mv fimo_temp.txt fimo.txt
sed 's/\t*$//' fimo.txt > fimo_temp.txt
mv fimo_temp.txt fimo.txt
