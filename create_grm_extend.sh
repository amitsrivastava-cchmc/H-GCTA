#!/bin/bash

#Commands in line 5-15 are specific to HPC (IBM-LSF). To run this script on any machine please comment out/remove the commonds in line 5-15.

#BSUB -M 65536
#BSUB -W 48:00
#BSUB -n 8
#BSUB -R "span[ptile=2]"
#BSUB -e /data/predataSamit/empirical_data_analysis/pooled_data/%J.err
#BSUB -o /data/predataSamit/empirical_data_analysis/pooled_data/%J.out

#load libraries

module load gcc/5.4.0
module load glibc/2.14
module load R/3.6.1

read -p "please provide the working directory with full path" wd
cd $wd

#merge plink binary files from individual datasets (merge_list${file}.txt is a list of binary filesets (prefix) for other datasets with their path)

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
	if [ -f "merge_list_${file}.txt" ]
	then
		fname=$(head -n 1 "merge_list_${file}.txt");
		filepath=$(dirname "$fname")
		plink --bfile ${filepath}/ALSPAC_imputed_autosomes_extracted_${file}_comm --merge-list merge_list_${file}.txt --indiv-sort 0 --keep-allele-order --make-bed --out pooled_data_${file}_comm --threads 16;
	else
		echo -e "merge_list_${file}.txt must be present in $(pwd)";
		exit 0
	fi
done

#create GRM and modify them

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       plink --bfile pooled_data_${file}_comm --make-grm-bin --out pooled_data_${file}_comm_grm --threads 16;
done

for file in $(echo hapM1 hapM2 hapP1); do
	n_lines=$(< pooled_data_${file}_comm.fam wc -l);
	divide_grm pooled_data_${file}_comm_grm.grm.bin pooled_data_${file}_comm_grm.grm.bin_tmp $n_lines 2;
	mv pooled_data_${file}_comm_grm.grm.bin_tmp pooled_data_${file}_comm_grm.grm.bin;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
        gcta64 --grm pooled_data_${file}_comm_grm --grm-cutoff 0.05 --make-grm-bin --out tmp_${file}_grm0.05 --threads 16;
done

#extract mother-child pairs with set cut-off for relatedness

if [ -f "pooled_mduos_extract_list.txt" ]
then
	awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' tmp_genoF_grm0.05.grm.id tmp_hapM1_grm0.05.grm.id > t1.txt
	awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t1.txt tmp_hapM2_grm0.05.grm.id > t2.txt
	awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t2.txt tmp_hapP1_grm0.05.grm.id > t3.txt
	awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t3.txt pooled_mduos_extract_list.txt > t4.txt
	awk 'NR==FNR{x[$1]++; next} $3 in x{print $0}' tmp_genoM_grm0.05.grm.id t4.txt > pooled_mduos_comm_grm0.05_extract_list.txt
	awk '{print $1, "\t", $1}' pooled_mduos_comm_grm0.05_extract_list.txt > pooled_fetus_comm_grm0.05_extract_list.txt
	awk '{print $3, "\t", $3}' pooled_mduos_comm_grm0.05_extract_list.txt > pooled_mother_comm_grm0.05_extract_list.txt
	cat pooled_mother_comm_grm0.05_extract_list.txt pooled_fetus_comm_grm0.05_extract_list.txt > pooled_mf_comm_grm0.05_extract_list.txt
	awk '{print $1, "\t", $1}' pooled_mduos_comm_grm0.05_extract_list.txt > pooled_fetus_comm_grm0.05_extract_list.keep
	awk '{print $3, "\t", $3}' pooled_mduos_comm_grm0.05_extract_list.txt > pooled_mother_comm_grm0.05_extract_list.keep
	rm t1.txt t2.txt t3.txt t4.txt
else
	echo -e "pedigree file of duos/trios must be present in $(pwd)";
	exit 0
fi

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
        rm tmp_${file}_grm0.05.*;
done

gcta64 --grm pooled_data_genoM_comm_grm --keep pooled_mother_comm_grm0.05_extract_list.txt --make-grm-bin --out pooled_data_genoM_comm_grm0.05 --threads 16;

for file in $(echo genoF hapM1 hapM2 hapP1); do
        gcta64 --grm pooled_data_${file}_comm_grm --keep pooled_fetus_comm_grm0.05_extract_list.txt --make-grm-bin --out pooled_data_${file}_comm_grm0.05 --threads 16;
done

#create grm using M-GCTA approach

plink --bfile pooled_data_genoM_comm --bmerge pooled_data_genoF_comm --keep-allele-order --indiv-sort 0 --out pooled_data_mf_comm

plink --bfile pooled_data_mf_comm --make-grm-bin --out pooled_data_mf_comm_grm --threads 16;

#create grm among relatives i.e., off-diagnonals that are < 0.05 are set to 0

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
        gcta64 --grm pooled_data_${file}_comm_grm --make-bK 0.05 --out pooled_data_${file}_comm_bK_grm --threads 16;
done

#thininng SNPs using LDAK-thin model

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak5.1.linux --thin pooled_data_${file}_comm_thin --bfile pooled_data_${file}_comm --window-prune 0.98 --window-kb 100 --max-threads 16;
done

awk 'NR==FNR{x[$1]++; next} !($1 in x){print $0}' pooled_data_genoM_comm_thin.in pooled_data_genoF_comm_thin.in > tmp.txt
cat pooled_data_genoM_comm_thin.in tmp.txt | sort -h -k 1,1 > pooled_data_mf_comm_thin.in

awk 'NR==FNR{x[$1]++; next} !($1 in x){print $0}' pooled_data_genoM_comm_thin.out pooled_data_genoF_comm_thin.out > tmp.txt
cat pooled_data_genoM_comm_thin.out tmp.txt | sort -h -k 1,1 > pooled_data_mf_comm_thin.out
rm tmp.txt

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       awk '{print $1, 1}' pooled_data_${file}_comm_thin.in > weights_${file}.thin;
done

#calculate weights of SNPs (LDAK-weights mdodel)

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       for j in $(seq 1 22); do mkdir -p ./pooled_data_${file}_comm_sections/sections$j/; awk -v var=$j '{split($1, a, ":"); if(a[1] == var) print $1}' pooled_data_${file}_comm_thin.in > ./pooled_data_${file}_comm_sections/sections$j/thin.in;
       done
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       for j in $(seq 1 22); do awk -v var=$j '{split($1, a, ":"); if(a[1] == var) print $1}' pooled_data_${file}_comm_thin.out > ./pooled_data_${file}_comm_sections/sections$j/thin.out;
       done
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       for j in {1..22}; do
               ldak5.1.linux --cut-weights ./pooled_data_${file}_comm_sections/sections$j --bfile pooled_data_${file}_comm --chr $j --no-thin DONE --max-threads 16;
               ldak5.1.linux --calc-weights-all ./pooled_data_${file}_comm_sections/sections$j --bfile pooled_data_${file}_comm --chr $j --max-threads 16;
       done
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do 
       cat ./pooled_data_${file}_comm_sections/sections{1..22}/weights.short > ./pooled_data_${file}_comm_sections/weights_${file}.short
done

#calculate kinship matrix using LDAK-Thin model i.e. using SNPs after pruning

#using alpha = -1.0

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak5.1.linux --calc-kins-direct pooled_data_${file}_comm_thin_alpha-1.0_grm --bfile pooled_data_${file}_comm --weights weights_${file}.thin --power -1.0 --max-threads 16;
done

#using alpha = -0.25

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak5.1.linux --calc-kins-direct pooled_data_${file}_comm_thin_alpha-0.25_grm --bfile pooled_data_${file}_comm --weights weights_${file}.thin --power -0.25 --max-threads 16;
done

#calculate kinship matrix using LDAK-weights model (SNPs used for kinship calculation are generally less than the # of SNPs that were intially pruned-in i.e. SNPs in thin.in file, because many of them get zero weight and are ignored)

#using alpha = -1.0

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak5.1.linux --calc-kins-direct pooled_data_${file}_comm_ldak_weights_alpha-1.0_grm --bfile pooled_data_${file}_comm --weights ./pooled_data_${file}_comm_sections/weights_${file}.short --power -1.0 --max-threads 16;
done

#using alpha = -0.25

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak5.1.linux --calc-kins-direct pooled_data_${file}_comm_ldak_weights_alpha-0.25_grm --bfile pooled_data_${file}_comm --weights ./pooled_data_${file}_comm_sections/weights_${file}.short --power -0.25 --max-threads 16;
done

calculate kinship matrix using GCTA model i.e. ignoring weights of SNPs

#using alpha = -1.0

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak5.1.linux --calc-kins-direct pooled_data_${file}_comm_ldak_gcta_apha-1.0_grm --bfile pooled_data_${file}_comm --ignore-weights YES --power -1.0 --max-threads 16;
done

#using alpha = -0.25

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
       ldak5.1.linux --calc-kins-direct pooled_data_${file}_comm_ldak_gcta_alpha-0.25_grm --bfile pooled_data_${file}_comm --ignore-weights YES --power -0.25 --max-threads 16;
done

#filter individuals with set grm-cutoff (e.g. 0.05)

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak5.1.linux --filter pooled_data_${file}_comm_thin_alpha-1.0_grm0.05 --grm pooled_data_${file}_comm_thin_alpha-1.0_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak5.1.linux --filter pooled_data_${file}_comm_thin_alpha-0.25_grm0.05 --grm pooled_data_${file}_comm_thin_alpha-0.25_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak5.1.linux --filter pooled_data_${file}_comm_ldak_weights_alpha-1.0_grm0.05 --grm pooled_data_${file}_comm_ldak_weights_alpha-1.0_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak5.1.linux --filter pooled_data_${file}_comm_ldak_weights_alpha-0.25_grm0.05 --grm pooled_data_${file}_comm_ldak_weights_alpha-0.25_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak5.1.linux --filter pooled_data_${file}_comm_ldak_gcta_alpha-1.0_grm0.05 --grm pooled_data_${file}_comm_ldak_gcta_alpha-1.0_grm --max-rel 0.05;
done

for file in $(echo genoM genoF hapM1 hapM2 hapP1); do
       ldak5.1.linux --filter pooled_data_${file}_comm_ldak_gcta_alpha-0.25_grm0.05 --grm pooled_data_${file}_comm_ldak_gcta_alpha-0.25_grm --max-rel 0.05;
done

#extract mother-child pairs with Kinship coefficient < 0.05. Respective GRMs using LDAK Thin/weights model are used for this purpose.

if [ -f "7.txt" ]
then
    for i in $(echo thin ldak_weights ldak_gcta); do
        for j in $(echo -1.0 -0.25); do
            awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' pooled_data_genoF_comm_${i}_alpha${j}_grm0.05.keep pooled_data_hapM1_comm_${i}_alpha${j}_grm0.05.keep > t1.txt;
			awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t1.txt pooled_data_hapM2_comm_${i}_alpha${j}_grm0.05.keep > t2.txt;
			awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t2.txt pooled_data_hapP1_comm_${i}_alpha${j}_grm0.05.keep > t3.txt;
			awk 'NR==FNR{x[$1]++; next} $1 in x{print $0}' t3.txt pooled_mduos_extract_list.txt > t4.txt;
			awk 'NR==FNR{x[$1]++; next} $3 in x{print $0}' pooled_data_genoM_comm_${i}_alpha${j}_grm0.05.keep t4.txt > pooled_mduos_comm_${i}_alpha${j}_grm0.05_extract_list.txt;
			awk '{print $1, "\t", $1}' pooled_mduos_comm_${i}_alpha${j}_grm0.05_extract_list.txt > pooled_fetus_comm_${i}_alpha${j}_grm0.05_extract_list.keep;
			awk '{print $3, "\t", $3}' pooled_mduos_comm_${i}_alpha${j}_grm0.05_extract_list.txt > pooled_mother_comm_${i}_alpha${j}_grm0.05_extract_list.keep;
		done
	done
	rm t1.txt t2.txt.t3.txt t4.txt;
else
	echo -e "pedigree file of duos/trios needed in the same directory";
	exit 0
fi

#modify joint maternal-fetal GRMs

Rscript mod_mf_grm.R $wd pooled_data_mf_comm_grm
Rscript mod_mf_grm.R $wd pooled_data_mf_comm_bK_grm

Rscript mod_mf_grm.R $wd pooled_data_mf_comm_ldak_gcta_alpha-1.0_grm
Rscript mod_mf_grm.R $wd pooled_data_mf_comm_ldak_gcta_alpha-0.25_grm

Rscript mod_mf_grm.R $wd pooled_data_mf_comm_thin_alpha-1.0_grm
Rscript mod_mf_grm.R $wd pooled_data_mf_comm_thin_alpha-0.25_grm

Rscript mod_mf_grm.R $wd pooled_data_mf_comm_ldak_weights_alpha-1.0_grm
Rscript mod_mf_grm.R $wd pooled_data_mf_comm_ldak_weights_alpha-0.25_grm

#calculate PCs using all samples

for file in $(echo genoM genoF hapM1 hapM2 hapP1 mf); do
	plink --bfile pooled_data_${file}_comm --indep-pairwise 200 50 0.2 --out pooled_data_${file}_comm_ld --threads 16;
	plink --bfile pooled_data_${file}_comm --extract pooled_data_${file}_comm_ld.prune.in --pca 20 --out pooled_data_${file}_comm_pc --threads 16;
done

#calculate pcs using unrelated samples (kinship coefficient > 0.05)

plink --bfile pooled_data_genoM_comm --keep pooled_mother_comm_grm0.05_extract_list.txt --indep-pairwise 200 50 0.2 --out pooled_data_genoM_comm_kinship0.05_ld --threads 16;

plink --bfile pooled_data_genoM_comm --keep pooled_mother_comm_grm0.05_extract_list.txt --extract pooled_data_genoM_comm_kinship0.05_ld.prune.in --pca 20 --out pooled_data_genoM_comm_kinship0.05_pc --threads 16;

for file in $(echo genoF hapM1 hapM2 hapP1); do
	plink --bfile pooled_data_${file}_comm --keep pooled_fetus_comm_grm0.05_extract_list.txt --indep-pairwise 200 50 0.2 --out pooled_data_${file}_comm_kinship0.05_ld --threads 16;
	plink --bfile pooled_data_${file}_comm --keep pooled_fetus_comm_grm0.05_extract_list.txt --extract pooled_data_${file}_comm_kinship0.05_ld.prune.in --pca 20 --out pooled_data_${file}_comm_kinship0.05_pc --threads 16;
done

plink --bfile pooled_data_mf_comm --keep pooled_mf_comm_grm0.05_extract_list.txt --indep-pairwise 200 50 0.2 --out pooled_data_mf_comm_kinship0.05_ld --threads 16;

plink --bfile pooled_data_mf_comm --keep pooled_mf_comm_grm0.05_extract_list.txt --extract pooled_data_mf_comm_kinship0.05_ld.prune.in --pca 20 --out pooled_data_mf_comm_kinship0.05_pc --threads 16;


#create matrix of pcs for analysis using H-GCTA and M-GCTA approach

awk -F " " 'BEGIN{OFS=="\t"} NR==FNR{a[$1 OFS $2] = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12; next} {for(i in a){if(i == $1 OFS $2) print a[i] OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12} }' pooled_data_hapM1_comm_pc.eigenvec pooled_data_hapM2_comm_pc.eigenvec > tmp.txt

awk -F " " 'BEGIN{OFS=="\t"} NR==FNR{a[$1 OFS $2] = $0; next} {for(i in a){if(i == $1 OFS $2) print a[i] OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12} }' tmp.txt pooled_data_hapP1_comm_pc.eigenvec > pooled_data_hap_comm_pc.eigenvec

cut -d " " -f 1-12 pooled_data_genoF_comm_pc.eigenvec | paste - <(cut -d " " -f 3-12 pooled_data_genoM_comm_pc.eigenvec) > pooled_data_genoFM_comm_pc.eigenvec

awk -F " " 'BEGIN{OFS=="\t"} NR==FNR{a[$1 OFS $2] = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12; next} {for(i in a){if(i == $1 OFS $2) print a[i] OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12} }' pooled_data_hapM1_comm_kinship0.05_pc.eigenvec pooled_data_hapM2_comm_kinship0.05_pc.eigenvec > tmp.txt

awk -F " " 'BEGIN{OFS=="\t"} NR==FNR{a[$1 OFS $2] = $0; next} {for(i in a){if(i == $1 OFS $2) print a[i] OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12} }' tmp.txt pooled_data_hapP1_comm_kinship0.05_pc.eigenvec > pooled_data_hap_comm_kinship0.05_pc.eigenvec

cut -d " " -f 1-12 pooled_data_genoF_comm_kinship0.05_pc.eigenvec | paste - <(cut -d " " -f 3-12 pooled_data_genoM_comm_kinship0.05_pc.eigenvec) > pooled_data_genoFM_comm_kinship0.05_pc.eigenvec

rm tmp.txt
