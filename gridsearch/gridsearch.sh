#! /bin/bash

# execute this in /MMCI/MS/.../FAME/

#definition of program rebuild and run routine
build_progs() {
	q=$1
	# m=$2
	t=$2
	# echo $m
	echo $q
	echo $t
	# sed -i 's/^constexpr unsigned int WINLEN \=.*/constexpr unsigned int WINLEN \= '"$m"'\;/' CONST.h
	sed -i 's/^constexpr uint64_t KMERCUTOFF \=.*/constexpr uint64_t KMERCUTOFF \= '"$t"'\;/' CONST.h
	sed -i 's/^constexpr uint16_t QTHRESH \=.*/constexpr uint16_t QTHRESH \= '"$q"'\;/' CONST.h
	make clean
	make
	mv FAME "FAME_${q}_${t}"
}

run_par () {
	q=$1
	# m=$2
	t=$2
	/usr/bin/time ./FAME"_${q}_${t}" --load_index gridindex"_${t}" -r1 Synth/paired_grid_CHR22_p1.fastq -r2 Synth/paired_grid_CHR22_p2.fastq -o gridsearch2/FullOutput_q"${q}_t${t}" | tee gridsearch2/FullLog_q${q}_t${t}.txt
	Rscript results/plot_it_my.R Synth/paired_grid_CHR22_cpginfo_fwd.tsv Synth/paired_grid_CHR22_cpginfo_rev.tsv gridsearch2/FullOutput_q${q}_t${t}_cpg.tsv results2/FAME_Full_prediction_pairedend_q"${q}_t${t}".pdf "Methylation rates for mapping paired reads on CHR22" | tee gridsearch2/FullLog_plot_q${q}_t${t}
	# rm gridsearch2/FullOutput_q"${q}_t${t}"_cpg.tsv
	# rm FAME"_${q}_${t}"
}

build_ind () {
	# m=$1
	t=$1
	/usr/bin/time ./FAME"_4_${t}" --genome ../BS_hgref/hg19_ref/hg19.fa --store_index gridindex"_${t}" --human_opt | tee gridsearch2/FullLog_index_t${t}.txt
}

for q in "4" "5" "6"
do
	for t in "200" "500" "1500" "3500" "5000" "10000"
 	do
		# for m in "1024" "2048" "4096" "8192"
		# do 			
			build_progs $q $t
			sleep 0.1
		# done
	done
done 

for t in "200" "500" "1500" "3500" "5000" "10000"
do
	build_ind $t
	for q in "4" "5" "6"
 	do
		run_par $q $t
	done
	# rm gridindex"_${t}"
	# rm gridindex"_${t}"_strands
done
