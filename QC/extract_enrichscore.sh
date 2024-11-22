for file in *enrichscore; do
  prefix=$(echo "$file" | cut -d'_' -f1-3)
  echo "====== ${prefix} ====="
  cat "$file"
done