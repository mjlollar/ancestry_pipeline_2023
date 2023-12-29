# ancestry_pipeline_2023

Command to get Sync file from mpileup:

java -ea -Xmx7g -jar <path>/mpileup2sync.jar --input <chr>.mpileup --output <chr>.sync --fastq-type sanger --min-qual 20 --threads 8
