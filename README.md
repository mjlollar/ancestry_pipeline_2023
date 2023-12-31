# ancestry_pipeline_2023

Command to get Sync file from mpileup:

java -ea -Xmx7g -jar "your path"/mpileup2sync.jar --input "chrom".mpileup --output "chrom".sync --fastq-type sanger --min-qual 20 --threads 8


Use populate_panel_counts_recomb_12323.py to get ahmm_in.panel files

Example input:
python3 populate_panel_counts_recomb_123123.py --p X.panel --s F1_chromX.sync --c X

Expected output: X.ahmm_in.panel
