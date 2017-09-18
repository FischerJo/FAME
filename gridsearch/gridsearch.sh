#! /bin/bash

# execute this in /MMCI/MS/.../Metal/

#definition of program rebuild and run routine
build_progs() {
	k=$1
	m=$2
	t=$3
	sed -i 's/^constexpr unsigned int WINLEN \=.*/constexpr unsigned int WINLEN \= '"$m"'\;/' CONST.h
	sed -i 's/^constexpr unsigned int KMERLEN \=.*/constexpr unsigned int KMERLEN \= '"$k"'\;/' CONST.h
	#sed -i 's/^constexpr uint64_t KMERCUTOFF \=.*/constexpr uint64_t KMERCUTOFF \= '"$t"'\;/' CONST.h
	make clean
	make
	mv Metal "Metal_${k}_${m}_${t}"
}

run_par () {
	k=$1
	m=$2
	t=$3
	#/usr/bin/time ./Metal"_${k}_${m}_${t}" --genome ../BS_hgref/GCA_000001405.25_GRCh38.p10_genomic.fna --store_index gridindex"_${k}_${m}_${t}"
	/usr/bin/time ./Metal"_${k}_${m}_${t}" --genome ../BS_hgref/GCA_000001405.25_GRCh38.p10_genomic.fna --store_index gridindex"_${k}_${m}_${t}" --no_loss
	/usr/bin/time ./Metal"_${k}_${m}_${t}" --load_index gridindex"_${k}_${m}_${t}" -r Synth/reads_2err_chr22_diffMeth_fwd.fastq -o gridsearch2/output_k"${k}_m${m}_t${t}" | tee gridsearch2/log_k${k}_m${m}_t${t}.txt
	rm gridindex"_${k}_${m}_${t}"
	rm gridindex"_${k}_${m}_${t}"_strands
	awk '$1 == 21' <gridsearch2/output_k"${k}_m${m}_t${t}"_cpg.tsv >output_k"${k}_m${m}_t${t}"_CHR22.tsv
	R3script results/plot_it_my.R Synth/reads_2err_chr22_diffMeth_cpginfo_fwd_SORTED.tsv output_k${k}_m${m}_t${t}_CHR22.tsv results2/Metal_prediction_diffMethData_k"${k}_m${m}_t${t}".pdf "Methylation rates for mapping on CHR22"
	rm output_k"${k}_m${m}_t${t}"_CHR22.tsv
	rm Metal"_${k}_${m}_${t}"
}

for k in "25" "28"
do
	#for t in "500" "2000" "5000"
	#do
		for m in "512" "1024" "2048"
		do
			#build_progs $k $m $t
			build_progs $k $m "no"
			#sleep 0.1
		done
	#done
done 

for k in "25" "28"
do
	#for t in "500" "2000" "5000"
	#do
		for m in "512" "1024" "2048"
		do
			#run_par $k $m $t
			run_par $k $m "no"
		done
	#done
done 
