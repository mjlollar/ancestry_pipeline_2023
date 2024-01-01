#### MJL 12/31/23
#### ANCESTRY PIPELINE POST BWA,INDEX,REALIGN

# list the RIL bam files
ls *.bam > RIL_list.txt # list all the bam files (all the bam files should be from the RILs to be analized)
sed -i 's/_realign.bam//' RIL_list.txt # edit the file to use it as input to the for loop below

for RIL in $(cat RIL_list.txt)
do	
	# Index the line (produces .bam.bai file in the directory containing bams) This is used by mpileup call
	samtools index ${RIL}_realign.bam

	# Each number in names corresponds to chromosome index, respectively: X, 2L, 2R, 3L, 3R
	names="4 3 7 5 8"
	
	# Generate mpileup files
	for name in $names
	do
		samtools mpileup -q 20 -r $name ${RIL}_realign.bam > ${name}.mpileup 
	done

	# Generate sync files
	for name in $names
	do
		java -ea -Xmx7g -jar /raid10/Tiago/PigmCages/scripts/alignment_software/mpileup2sync.jar --input ${name}.mpileup --output ${name}.sync --fastq-type sanger --min-qual 20 --threads 8
		#rm ${name}.mpileup #remove mpileup intermediate
	done

	# Generate ahmm_input panels from population panel and sync files
	
	#X
	python3 populate_panel_counts_recomb_123123.py --p X.panel --s 4.sync --c X
	#python3 rec_map.py X_recomb_0.1gen.csv X.ahmm_in.panel
	#mv X.ahmm_in.panel2 X.ahmm_in.panel
	#2L
	python3 populate_panel_counts_recomb_123123.py --p 2L.panel --s 3.sync --c 2L
	#python3 rec_map.py 2L_recomb_0.1gen.csv 2L.ahmm_in.panel
	#mv 2L.ahmm_in.panel2 2L.ahmm_in.panel
	#2R
	python3 populate_panel_counts_recomb_123123.py --p 2R.panel --s 7.sync --c 2R
	#python3 rec_map.py 2R_recomb_0.1gen.csv 2R.ahmm_in.panel
	#mv 2R.ahmm_in.panel2 X.ahmm_in.panel
	#3L
	python3 populate_panel_counts_recomb_123123.py --p 3L.panel --s 5.sync --c 3L
	#python3 rec_map.py 3L_recomb_0.1gen.csv 3L.ahmm_in.panel
	#mv 3L.ahmm_in.panel2 3L.ahmm_in.panel
	#3R
	python3 populate_panel_counts_recomb_123123.py --p 3R.panel --s 8.sync --c 3R
	#python3 rec_map.py 3R_recomb_0.1gen.csv 3R.ahmm_in.panel
	#mv 3R.ahmm_in.panel2 3R.ahmm_in.panel
	
	#rm *.sync #Remove intermediate sync files

	#### Create the ahmm sample files
	#X Chrom ploidy 1 (use if only males in pool/single male(hemizygous))
	ls ${RIL}_realign.bam | perl -p -e 's/\n/\t1\n/' > ahmm_X_in.samples
	#Ploidy 2 for remaining chroms
	ls ${RIL}_realign.bam | perl -p -e 's/\n/\t2\n/' > ahmm_in.samples

	#### Run Ancestry_hmm to generate
	#X Chrom
	/raid10/Tiago/RILS/Ancestry_HMM/src/ancestry_hmm -i X.ahmm_in.panel -s ahmm_X_in.samples -a 2 0.5 0.5 -p 0 10000 0.5 -p 1 1 0.5
	mv ${RIL}_realign.bam.posterior posteriors/X/${RIL}.posterior
	#Autosomes
	chroms="2L 2R 3L 3R"
	for name in $chroms
	do
		#### email Matt if there are any questions about ancestry_hmm parameters
		/raid10/Tiago/RILS/Ancestry_HMM/src/ancestry_hmm -i ${name}.ahmm_in.panel -s ahmm_in.samples -a 2 0.5 0.5 -p 0 10000 0.5 -p 1 1 0.5 # May need to change Ancestry_hmm path
		### Move files into a folder "posteriors" and sub-directories (ex. "2L"). The posteriors must have the naming scheme "<line>.posterior
		### and be in a subdirectory with names X, 2L, 2R, 3L, and 3R for the script used in the next steps.
		mv ${RIL}_realign.bam.posterior posteriors/${name}/${RIL}.posterior
	done
	
	### Remove unecessary intermediate files
	#rm *.ahmm_in.panel
	#rm ahmm_in.samples
	#rm ahmm_X_in.samples
done
