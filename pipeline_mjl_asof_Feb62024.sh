#### MJL last update 01/29/23

##Preliminary alignment commands (run on chtc)
### Index reference genome
bwa index dmel-r5.57.fasta #use relevant flags if using ref genome >2GB
### Align sample
# -M and -R flags are needed for downstream picard/GATK compatibility
bwa mem -M -R '@RG\tID:sample1\tLB:Lane1\tPL:ILLUMINA\tSM:F1' DmelRef.fasta ${line}_R1.fastq ${line}_R2.fastq > ${line}_align.sam

#### conda activate condamatt for next set of commands ####
### Convert Sam to Bam and coordinate sort
java -Xmx5g -jar picard.jar SortSam -I ${line}_align.sam -O ${line}_align.bam -SORT_ORDER coordinate  
### Mark duplicates
java -Xmx5g -jar picard.jar MarkDuplicates -I ${line}_align.bam -O ${line}_align_md.bam -M duplicatemetrics.txt
###Build bai index for sorted bam
java -Xmx5g -jar picard.jar BuildBamIndex -I ${line}_align_md.bam

#### conda activate tiagojav8 for next set of commands ####
### Realign around Indels
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -I F1_align_md.bam -o F1.intervals -R DmelRef.fasta
java -jar GenomeAnalysisTK.jar -T IndelRealigner -I F1_align_md.bam -targetIntervals F1.intervals -R DmelRef.fasta -o F1_realign.bam 

### List the bam files
ls *.bam > bam_list.txt
sed -i 's/_realign.bam//' bam_list.txt # edit the file to use it as input to the for loop below

for BAM in $(cat bam_list.txt)
do	
	### Index the line (produces .bai file)
	samtools index ${BAM}_realign.bam

	### Each number in names corresponds to chromosome index, respectively: X, 2L, 2R, 3L, 3R
	names="4 3 7 5 8"
	### Generate mpileup files
	for name in $names
	do
		samtools mpileup -q 20 -Q 20 -r $name ${BAM}_realign.bam > ${name}.mpileup
	done

	### Generate sync files
	for name in $names
	do
		java -ea -Xmx7g -jar mpileup2sync.jar --input ${name}.mpileup --output ${name}.sync --fastq-type sanger --min-qual 20 --threads 8
		#mv ${name}.mpileup mpileups/${BAM}.${name}.mpileup
		#mv ${name}.sync syncs/${BAM}.${name}.sync
	done

	### Generate ahmm_input panels from population panel and sync files
  	### Use rec_map.py script to adjust recombination rates if desired
  	### conda activate sklearnpy (needs pandas)
	## X
	python3 populate_panel_counts_asof_020624.py --p X.panel --s ${BAM}.4.sync --c X --rc 0 --off -1
	#python3 rec_map.py X_recomb_0.1gen.csv X.ahmm_in.panel
	#mv X.ahmm_in.panel2 X.ahmm_in.panel
	## 2L
	python3 populate_panel_counts_asof_020624.py --p 2L.panel --s ${BAM}.3.sync --c 2L --rc 0 --off -1
	#python3 rec_map.py 2L_recomb_0.1gen.csv 2L.ahmm_in.panel
	#mv 2L.ahmm_in.panel2 2L.ahmm_in.panel
	## 2R
	python3 populate_panel_counts_asof_020624.py --p 2R.panel --s ${BAM}.7.sync --c 2R --rc 0 --off -1
	#python3 rec_map.py 2R_recomb_0.1gen.csv 2R.ahmm_in.panel
	#mv 2R.ahmm_in.panel2 2R.ahmm_in.panel
	## 3L
	python3 populate_panel_counts_asof_020624.py --p 3L.panel --s ${BAM}.5.sync --c 3L --rc 0 --off -1
	#python3 rec_map.py 3L_recomb_0.1gen.csv 3L.ahmm_in.panel
	#mv 3L.ahmm_in.panel2 3L.ahmm_in.panel
	## 3R
	python3 populate_panel_counts_asof_020624.py --p 3R.panel --s ${BAM}.8.sync --c 3R --rc 0 --off -1
	#python3 rec_map.py 3R_recomb_0.1gen.csv 3R.ahmm_in.panel
	#mv 3R.ahmm_in.panel2 3R.ahmm_in.panel

	mv X.ahmm_in.panel ${BAM}.X.ahmm_in.panel
	mv 2L.ahmm_in.panel ${BAM}.2L.ahmm_in.panel
	mv 2R.ahmm_in.panel ${BAM}.2R.ahmm_in.panel
	mv 3L.ahmm_in.panel ${BAM}.3L.ahmm_in.panel
	mv 3R.ahmm_in.panel ${BAM}.3R.ahmm_in.panel
	#rm *.sync #Remove intermediate sync files
	#rm *.mpileup #Remove intermediate mpileups

	### Create the ahmm sample files
	#X Chrom ploidy 1 for X (use if only males in pool/single male(hemizygous))
	ls ${BAM}_realign.bam | perl -p -e 's/\n/\t1\n/' > ahmm_X_in.samples
	#Ploidy 2 for remaining chroms
	ls ${BAM}_realign.bam | perl -p -e 's/\n/\t2\n/' > ahmm_in.samples

	#### Run Ancestry_hmm to generate posterior files
	#X Chrom
	ancestry_hmm -i ${BAM}.X.ahmm_in.panel -s ahmm_X_in.samples -a 2 0.5 0.5 -p 0 2 0.5 -p 1 2 0.5
	mv ${BAM}_realign.bam.posterior posteriors/X/${BAM}.posterior
	#Autosomes
	chroms="2L 2R 3L 3R"
	for name in $chroms
	do
		ancestry_hmm -i ${BAM}.${name}.ahmm_in.panel -s ahmm_in.samples -a 2 0.5 0.5 -p 0 2 0.5 -p 1 2 0.5 # May need to change Ancestry_hmm path
		## Move files into a folder "posteriors" and sub-directories (ex. "2L"). The posteriors must have the naming scheme "<line>.posterior
		## and be in a subdirectory with names X, 2L, 2R, 3L, and 3R for the script used in the next steps.
		mv ${BAM}_realign.bam.posterior posteriors/${name}/${BAM}.posterior
	done
	### Remove unecessary intermediate files
	rm ahmm_in.samples
	rm ahmm_X_in.samples
done

# From this output run "RIL_genotypes_windows.pl" script in posterior directory (with relevant window files) to call ancestry

