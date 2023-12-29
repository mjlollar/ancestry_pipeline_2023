# ancestry_pipeline_2023

Command to get Sync file from mpileup:

java -ea -Xmx7g -jar "your path"/mpileup2sync.jar --input "chrom".mpileup --output "chrom".sync --fastq-type sanger --min-qual 20 --threads 8
