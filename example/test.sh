time (fast2bRAD-M pipeline \
  --mode db-only \
  --genome-list example_db_data/example_genomes_list.txt \
  --taxonomy example_db_data/example_genomes_taxonomy.txt \
  --site BcgI \
  --level species \
  --outdir example_db \
  --prefix example_db_data \
  --gscore 5 \
  --threads 4 \
  --pc 8 \
  --min-qual 15 \
  --resume no \
  --ko-mapping example_db_data/example_ko_mapping.txt)

time (fast2bRAD-M pipeline \
  --mode sample-only \
  --samples example_sample_data/example_samples_list.txt \
  --database example_db/02_db_qual \
  --site BcgI \
  --level species \
  --outdir example_results \
  --prefix example \
  --gscore 5 \
  --threads 4 \
  --pc 8 \
  --min-qual 15 \
  --resume no \
  --ko-mapping example_db_data/example_ko_mapping.txt)
