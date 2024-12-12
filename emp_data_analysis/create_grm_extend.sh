#!/bin/bash

#Line 5-16 are commented out to run this script on any personal computer which has R and Rscript in computer's $PATH.

#BSUB -M 65536
#BSUB -W 48:00
#BSUB -n 8
#BSUB -R "span[ptile=2]"
#BSUB -e /data/predataSamit/empirical_data_analysis/%J.err
#BSUB -o /data/predataSamit/empirical_data_analysis/%J.out

#load libraries

#module load gcc/5.4.0
#module load glibc/2.14
#module load R/4.1.3

#add directory with binaries to $PATH

binpath=$(realpath ../bin)
export PATH="$PATH:$binpath"

#save current directory as variable

current_dir=$(pwd)

#read and go to the working directory

read -p "please provide the working directory with full path " wd
read -p "DO YOU WANT TO WORK WITH SINGLE CHROMOSOME, IF YES PLEASE PROVIDE CHR NUMBER " chr

cd $wd

#merge plink binary files from individual datasets (merge_list${file}.txt is a list of binary filesets (prefix) for other datasets with their path).
#This step is not necessary for individual datasets.

#for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
#	if [ -f "merge_list_${file}.txt" ]
#	then
#		fname=$(head -n 1 "merge_list_${file}.txt");
#		filepath=$(dirname "$fname")
#		plink --bfile ${filepath}/test1_extracted_${file} --merge-list merge_list_${file}.txt --indiv-sort 0 --keep-allele-order --make-bed --out test_extracted_${file} --threads 16;
#	else
#		echo -e "merge_list_${file}.txt must be present in $(pwd)";
#		exit 0
#	fi
#done

#create GRM and modify them
#Although GRMs corresponding to maternal, fetal genotypes and maternal transmitted, non-transmitted and paternal transmitted haplotypes are already created via create_grm.sh, these steps are repeated here to maintain consistency trhoughout the script. These steps can be commented out, if not needed.

if [ "$chr" ]
then
	for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
		plink --bfile test_extracted_${file} --chr "$chr" --make-grm-bin --out test_extracted_${file}_grm --threads 16;
	done
else
	for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
		plink --bfile test_extracted_${file} --make-grm-bin --out test_extracted_${file}_grm --threads 16;
	done

fi

for file in $(echo hapM1 hapM2 hapP1); do
	n_lines=$(< test_extracted_${file}.fam wc -l);
	divide_grm test_extracted_${file}_grm.grm.bin test_extracted_${file}_grm.grm.bin_tmp $n_lines 2;
	mv test_extracted_${file}_grm.grm.bin_tmp test_extracted_${file}_grm.grm.bin;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
        gcta64 --grm test_extracted_${file}_grm --grm-cutoff 0.05 --make-grm-bin --out tmp_${file}_grm0.05 --threads 16;
done

#extract mother-child pairs with set cut-off for relatedness

if [ -f "mduos_extract_list.txt" ]
then
	awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' tmp_genoF_grm0.05.grm.id tmp_hapM1_grm0.05.grm.id > t1.txt
	awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t1.txt tmp_hapM2_grm0.05.grm.id > t2.txt
	awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t2.txt tmp_hapP1_grm0.05.grm.id > t3.txt
	awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t3.txt mduos_extract_list.txt > t4.txt
	awk 'NR==FNR{x[$1]++; next} $3 in x{print $0}' tmp_genoM_grm0.05.grm.id t4.txt > mduos_grm0.05_extract_list.txt
	awk '{print $1, "\t", $1}' mduos_grm0.05_extract_list.txt > fetus_grm0.05_extract_list.txt
	awk '{print $3, "\t", $3}' mduos_grm0.05_extract_list.txt > mother_grm0.05_extract_list.txt
	cat mother_grm0.05_extract_list.txt fetus_grm0.05_extract_list.txt > mf_grm0.05_extract_list.txt
	awk '{print $1, "\t", $1}' mduos_grm0.05_extract_list.txt > fetus_grm0.05_extract_list.keep
	awk '{print $3, "\t", $3}' mduos_grm0.05_extract_list.txt > mother_grm0.05_extract_list.keep
	rm t1.txt t2.txt t3.txt t4.txt
else
	echo -e "pedigree file of duos/trios must be present in $(pwd)";
	exit 0
fi

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
        rm tmp_${file}_grm0.05.*;
done

gcta64 --grm test_extracted_genoM_grm --keep mother_grm0.05_extract_list.txt --make-grm-bin --out test_extracted_genoM_grm0.05 --threads 16;

for file in $(echo genoF hapM1 hapM2 hapP1); do
        gcta64 --grm test_extracted_${file}_grm --keep fetus_grm0.05_extract_list.txt --make-grm-bin --out test_extracted_${file}_grm0.05 --threads 16;
done

#create grm using M-GCTA approach

plink --bfile test_extracted_genoM --bmerge test_extracted_genoF --keep-allele-order --indiv-sort 0 --out test_extracted_mf

plink --bfile test_extracted_mf --make-grm-bin --out test_extracted_mf_grm --threads 16;

#create grm among relatives i.e., off-diagnonals that are < 0.05 are set to 0

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
        gcta64 --grm test_extracted_${file}_grm --make-bK 0.05 --out test_extracted_${file}_bK_grm --threads 16;
done

#thininng SNPs using LDAK-thin model

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak6.linux --thin test_extracted_${file}_thin --bfile test_extracted_${file} --window-prune 0.98 --window-kb 100 --max-threads 16;
done

awk 'NR==FNR{x[$1]++; next} !($1 in x){print $0}' test_extracted_genoM_thin.in test_extracted_genoF_thin.in > tmp.txt
cat test_extracted_genoM_thin.in tmp.txt | sort -h -k 1,1 > test_extracted_mf_thin.in

awk 'NR==FNR{x[$1]++; next} !($1 in x){print $0}' test_extracted_genoM_thin.out test_extracted_genoF_thin.out > tmp.txt
cat test_extracted_genoM_thin.out tmp.txt | sort -h -k 1,1 > test_extracted_mf_thin.out
rm tmp.txt

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       awk '{print $1, 1}' test_extracted_${file}_thin.in > weights_${file}.thin;
done

#calculate weights of SNPs (LDAK-weights mdodel)

if [ "$chr" ]
then
	#for single chromosome data

	for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
		mkdir -p ./test_extracted_${file}_sections/sections${chr}
		awk -v var=${chr} '{split($1, a, ":"); if(a[1] == var) print $1}' test_extracted_${file}_thin.in > ./test_extracted_${file}_sections/sections${chr}/thin.in;
		awk -v var=${chr} '{split($1, a, ":"); if(a[1] == var) print $1}' test_extracted_${file}_thin.out > ./test_extracted_${file}_sections/sections${chr}/thin.out;
	done

	for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
		ldak6.linux --cut-weights ./test_extracted_${file}_sections/sections${chr} --bfile test_extracted_${file} --chr ${chr} --no-thin DONE --max-threads 16;
		ldak6.linux --calc-weights-all ./test_extracted_${file}_sections/sections20 --bfile test_extracted_${file} --chr ${chr} --allow-many-samples YES --max-threads 16;
		cp ./test_extracted_${file}_sections/sections${chr}/weights.short ./test_extracted_${file}_sections/weights_${file}.short
	done
else
	#for whole genome

	for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
		for j in $(seq 1 22); do
			mkdir -p ./test_extracted_${file}_sections/sections$j
			awk -v var=$j '{split($1, a, ":"); if(a[1] == var) print $1}' test_extracted_${file}_thin.in > ./test_extracted_${file}_sections/sections$j/thin.in;
			awk -v var=$j '{split($1, a, ":"); if(a[1] == var) print $1}' test_extracted_${file}_thin.out > ./test_extracted_${file}_sections/sections$j/thin.out;
		done
	done

	for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
		for j in {1..22}; do
			ldak6.linux --cut-weights ./test_extracted_${file}_sections/sections$j --bfile test_extracted_${file} --chr $j --no-thin DONE --max-threads 16;
			ldak6.linux --calc-weights-all ./test_extracted_${file}_sections/sections$j --bfile test_extracted_${file} --chr $j --allow-many-samples YES --max-threads 16;
		done
		cat ./test_extracted_${file}_sections/sections{1..22}/weights.short > ./test_extracted_${file}_sections/weights_${file}.short
	done
fi

#calculate kinship matrix using LDAK-Thin model i.e. using SNPs after pruning

#using alpha = -1.0

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak6.linux --calc-kins-direct test_extracted_${file}_thin_alpha-1.0_grm --bfile test_extracted_${file} --weights weights_${file}.thin --power -1.0 --max-threads 16;
done

#using alpha = -0.25

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak6.linux --calc-kins-direct test_extracted_${file}_thin_alpha-0.25_grm --bfile test_extracted_${file} --weights weights_${file}.thin --power -0.25 --max-threads 16;
done

#calculate kinship matrix using LDAK-weights model (SNPs used for kinship calculation are generally less than the # of SNPs that were intially pruned-in i.e. SNPs in thin.in file, because many of them get zero weight and are ignored)

#using alpha = -1.0

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak6.linux --calc-kins-direct test_extracted_${file}_ldak_weights_alpha-1.0_grm --bfile test_extracted_${file} --weights ./test_extracted_${file}_sections/weights_${file}.short --power -1.0 --max-threads 16;
done

#using alpha = -0.25

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak6.linux --calc-kins-direct test_extracted_${file}_ldak_weights_alpha-0.25_grm --bfile test_extracted_${file} --weights ./test_extracted_${file}_sections/weights_${file}.short --power -0.25 --max-threads 16;
done

#calculate kinship matrix using GCTA model i.e. ignoring weights of SNPs

#using alpha = -1.0

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak6.linux --calc-kins-direct test_extracted_${file}_ldak_gcta_alpha-1.0_grm --bfile test_extracted_${file} --ignore-weights YES --power -1.0 --max-threads 16;
done

#using alpha = -0.25

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak6.linux --calc-kins-direct test_extracted_${file}_ldak_gcta_alpha-0.25_grm --bfile test_extracted_${file} --ignore-weights YES --power -0.25 --max-threads 16;
done

#filter individuals with set grm-cutoff (e.g. 0.05)

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak6.linux --filter test_extracted_${file}_thin_alpha-1.0_grm0.05 --grm test_extracted_${file}_thin_alpha-1.0_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak6.linux --filter test_extracted_${file}_thin_alpha-0.25_grm0.05 --grm test_extracted_${file}_thin_alpha-0.25_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak6.linux --filter test_extracted_${file}_ldak_weights_alpha-1.0_grm0.05 --grm test_extracted_${file}_ldak_weights_alpha-1.0_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak6.linux --filter test_extracted_${file}_ldak_weights_alpha-0.25_grm0.05 --grm test_extracted_${file}_ldak_weights_alpha-0.25_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak6.linux --filter test_extracted_${file}_ldak_gcta_alpha-1.0_grm0.05 --grm test_extracted_${file}_ldak_gcta_alpha-1.0_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak6.linux --filter test_extracted_${file}_ldak_gcta_alpha-0.25_grm0.05 --grm test_extracted_${file}_ldak_gcta_alpha-0.25_grm --max-rel 0.05;
done

#extract mother-child pairs with Kinship coefficient < 0.05. Respective GRMs using LDAK Thin/weights model are used for this purpose.

if [ -f "mduos_extract_list.txt" ]
then
    for i in $(echo thin ldak_weights ldak_gcta); do
        for j in $(echo -1.0 -0.25); do
            awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' test_extracted_genoF_${i}_alpha${j}_grm0.05.keep test_extracted_hapM1_${i}_alpha${j}_grm0.05.keep > t1.txt;
			awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t1.txt test_extracted_hapM2_${i}_alpha${j}_grm0.05.keep > t2.txt;
			awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t2.txt test_extracted_hapP1_${i}_alpha${j}_grm0.05.keep > t3.txt;
			awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t3.txt mduos_extract_list.txt > t4.txt;
			awk 'NR==FNR{x[$1]++; next} $3 in x{print $0}' test_extracted_genoM_${i}_alpha${j}_grm0.05.keep t4.txt > mduos_${i}_alpha${j}_grm0.05_extract_list.txt;
			awk '{print $1, "\t", $1}' mduos_${i}_alpha${j}_grm0.05_extract_list.txt > fetus_${i}_alpha${j}_grm0.05_extract_list.keep;
			awk '{print $3, "\t", $3}' mduos_${i}_alpha${j}_grm0.05_extract_list.txt > mother_${i}_alpha${j}_grm0.05_extract_list.keep;
		done
	done
	rm t1.txt t2.txt.t3.txt t4.txt;
else
	echo -e "pedigree file of duos/trios needed in the same directory";
	exit 0
fi

#modify joint maternal-fetal GRMs

Rscript ${current_dir}/mod_mf_grm.R $wd test_extracted_mf_grm
Rscript ${current_dir}/mod_mf_grm.R $wd test_extracted_mf_bK_grm

Rscript ${current_dir}/mod_mf_grm.R $wd test_extracted_mf_ldak_gcta_alpha-1.0_grm
Rscript ${current_dir}/mod_mf_grm.R $wd test_extracted_mf_ldak_gcta_alpha-0.25_grm

Rscript ${current_dir}/mod_mf_grm.R $wd test_extracted_mf_thin_alpha-1.0_grm
Rscript ${current_dir}/mod_mf_grm.R $wd test_extracted_mf_thin_alpha-0.25_grm

Rscript ${current_dir}/mod_mf_grm.R $wd test_extracted_mf_ldak_weights_alpha-1.0_grm
Rscript ${current_dir}/mod_mf_grm.R $wd test_extracted_mf_ldak_weights_alpha-0.25_grm

#copy .details files generated along with maternal-fetal grms in ldak model as *_mother_grm.details, *_fetus_grm.details, and *_mf_grm.details

for model in $(echo thin ldak_weights ldak_gcta); do
	for val in $(echo -1.0 -0.25); do
		for ext in $(echo adjust details); do
			cp test_extracted_mf_${model}_alpha${val}_grm.grm.${ext} test_extracted_mf_${model}_alpha${val}_mother_grm.grm.${ext};
			cp test_extracted_mf_${model}_alpha${val}_grm.grm.${ext} test_extracted_mf_${model}_alpha${val}_fetus_grm.grm.${ext};
			cp test_extracted_mf_${model}_alpha${val}_grm.grm.${ext} test_extracted_mf_${model}_alpha${val}_mf_grm.grm.${ext};
		done

	done
done

#calculate PCs using all samples

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
	plink --bfile test_extracted_${file} --indep-pairwise 200 50 0.2 --out test_extracted_${file}_ld --threads 16;
	plink --bfile test_extracted_${file} --extract test_extracted_${file}_ld.prune.in --pca 20 --out test_extracted_${file}_pc --threads 16;
done

#calculate pcs using unrelated samples (kinship coefficient > 0.05)

plink --bfile test_extracted_genoM --keep mother_grm0.05_extract_list.txt --indep-pairwise 200 50 0.2 --out test_extracted_genoM_kinship0.05_ld --threads 16;

plink --bfile test_extracted_genoM --keep mother_grm0.05_extract_list.txt --extract test_extracted_genoM_kinship0.05_ld.prune.in --pca 20 --out test_extracted_genoM_kinship0.05_pc --threads 16;

for file in $(echo genoF hapM1 hapM2 hapP1); do
	plink --bfile test_extracted_${file} --keep fetus_grm0.05_extract_list.txt --indep-pairwise 200 50 0.2 --out test_extracted_${file}_kinship0.05_ld --threads 16;
	plink --bfile test_extracted_${file} --keep fetus_grm0.05_extract_list.txt --extract test_extracted_${file}_kinship0.05_ld.prune.in --pca 20 --out test_extracted_${file}_kinship0.05_pc --threads 16;
done

plink --bfile test_extracted_mf --keep mf_grm0.05_extract_list.txt --indep-pairwise 200 50 0.2 --out test_extracted_mf_kinship0.05_ld --threads 16;

plink --bfile test_extracted_mf --keep mf_grm0.05_extract_list.txt --extract test_extracted_mf_kinship0.05_ld.prune.in --pca 20 --out test_extracted_mf_kinship0.05_pc --threads 16;


#create matrix of pcs for analysis using H-GCTA and M-GCTA approach

awk -F " " 'BEGIN{OFS=="\t"} NR==FNR{a[$1 OFS $2] = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12; next} {for(i in a){if(i == $1 OFS $2) print a[i] OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12} }' test_extracted_hapM1_pc.eigenvec test_extracted_hapM2_pc.eigenvec > tmp.txt

awk -F " " 'BEGIN{OFS=="\t"} NR==FNR{a[$1 OFS $2] = $0; next} {for(i in a){if(i == $1 OFS $2) print a[i] OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12} }' tmp.txt test_extracted_hapP1_pc.eigenvec > test_extracted_hap_pc.eigenvec

cut -d " " -f 1-12 test_extracted_genoF_pc.eigenvec | paste - <(cut -d " " -f 3-12 test_extracted_genoM_pc.eigenvec) > test_extracted_genoFM_pc.eigenvec

awk -F " " 'BEGIN{OFS=="\t"} NR==FNR{a[$1 OFS $2] = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12; next} {for(i in a){if(i == $1 OFS $2) print a[i] OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12} }' test_extracted_hapM1_kinship0.05_pc.eigenvec test_extracted_hapM2_kinship0.05_pc.eigenvec > tmp.txt

awk -F " " 'BEGIN{OFS=="\t"} NR==FNR{a[$1 OFS $2] = $0; next} {for(i in a){if(i == $1 OFS $2) print a[i] OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12} }' tmp.txt test_extracted_hapP1_kinship0.05_pc.eigenvec > test_extracted_hap_kinship0.05_pc.eigenvec

cut -d " " -f 1-12 test_extracted_genoF_kinship0.05_pc.eigenvec | paste - <(cut -d " " -f 3-12 test_extracted_genoM_kinship0.05_pc.eigenvec) > test_extracted_genoFM_kinship0.05_pc.eigenvec

rm tmp.txt
