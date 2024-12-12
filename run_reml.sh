#!/bin/bash

#${1} = directory where genotype data and corresponding grms are saved e.g. "/home/user/H-GCTA/data.test/test_extracted"
#${2} = directory where pheotypes are present e.g "/home/user/H-GCTA/simulation/simulation_mf_with_ppoe_1.0_m1_0.0_mf_cor_1.0"
#${3} = number of iterations e.g. 100

#outputs will be written in ${2}

echo $(pwd);
iterations=${3}

#execute program

#run REML using GRM created by LDAK-Thin model

#using alpha = -1.0

for i in $(seq 1 "$iterations"); do
	pheno_m=$(find ${2}/pheno${i} -type f -name "*_m.phen");
	ldak6.linux --reml ${2}/pheno${i}/sim_ldak_thin_alpha-1.0_genoM_pheno${i}_reml --grm ${1}/test_extracted_genoM_thin_alpha-1.0_grm --pheno ${pheno_m} --keep ${1}/mother_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_thin_alpha-1.0_genoF_pheno${i}_reml --grm ${1}/test_extracted_genoF_thin_alpha-1.0_grm --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_thin_alpha-1.0_hap_pheno${i}_reml --mgrm ${1}/hgcta_ldak_thin_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_thin_alpha-1.0_mf_pheno${i}_reml --mgrm ${1}/mgcta_ldak_thin_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

#using alpha = -0.25

for i in $(seq 1 "$iterations"); do
        pheno_m=$(find ${2}/pheno${i} -type f -name "*_m.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_thin_alpha-0.25_genoM_pheno${i}_reml --grm ${1}/test_extracted_genoM_thin_alpha-0.25_grm --pheno ${pheno_m} --keep ${1}/mother_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_thin_alpha-0.25_genoF_pheno${i}_reml --grm ${1}/test_extracted_genoF_thin_alpha-0.25_grm --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_thin_alpha-0.25_hap_pheno${i}_reml --mgrm ${1}/hgcta_ldak_thin_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_thin_alpha-0.25_mf_pheno${i}_reml --mgrm ${1}/mgcta_ldak_thin_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

#run REML using GRM created by LDAK-Weights model

#using alpha = -1.0

for i in $(seq 1 "$iterations"); do
	pheno_m=$(find ${2}/pheno${i} -type f -name "*_m.phen");
	ldak6.linux --reml ${2}/pheno${i}/sim_ldak_weights_alpha-1.0_genoM_pheno${i}_reml --grm ${1}/test_extracted_genoM_ldak_weights_alpha-1.0_grm --pheno ${pheno_m} --keep ${1}/mother_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_weights_alpha-1.0_genoF_pheno${i}_reml --grm ${1}/test_extracted_genoF_ldak_weights_alpha-1.0_grm --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
	ldak6.linux --reml ${2}/pheno${i}/sim_ldak_weights_alpha-1.0_hap_pheno${i}_reml --mgrm ${1}/hgcta_ldak_weights_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_weights_alpha-1.0_mf_pheno${i}_reml --mgrm ${1}/mgcta_ldak_weights_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

#using alpha = -.025

for i in $(seq 1 "$iterations"); do
       pheno_m=$(find ${2}/pheno${i} -type f -name "*_m.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_weights_alpha-0.25_genoM_pheno${i}_reml --grm ${1}/test_extracted_genoM_ldak_weights_alpha-0.25_grm --pheno ${pheno_m} --keep ${1}/mother_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_weights_alpha-0.25_genoF_pheno${i}_reml --grm ${1}/test_extracted_genoF_ldak_weights_alpha-0.25_grm --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_weights_alpha-0.25_hap_pheno${i}_reml --mgrm ${1}/hgcta_ldak_weights_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_weights_alpha-0.25_mf_pheno${i}_reml --mgrm ${1}/mgcta_ldak_weights_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

#run REML using GRM created by LDAK-GCTA model

#using alpha = -1.0

for i in $(seq 1 "$iterations"); do
	pheno_m=$(find ${2}/pheno${i} -type f -name "*_m.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_gcta_alpha-1.0_genoM_pheno${i}_reml --grm ${1}/test_extracted_genoM_ldak_gcta_alpha-1.0_grm --pheno ${pheno_m} --keep ${1}/mother_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_gcta_alpha-1.0_genoF_pheno${i}_reml --grm ${1}/test_extracted_genoF_ldak_gcta_alpha-1.0_grm --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_gcta_alpha-1.0_hap_pheno${i}_reml --mgrm ${1}/hgcta_ldak_gcta_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_gcta_alpha-1.0_mf_pheno${i}_reml --mgrm ${1}/mgcta_ldak_gcta_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

#using alpha = -0.25

for i in $(seq 1 "$iterations"); do
        pheno_m=$(find ${2}/pheno${i} -type f -name "*_m.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_gcta_alpha-0.25_genoM_pheno${i}_reml --grm ${1}/test_extracted_genoM_ldak_gcta_alpha-0.25_grm --pheno ${pheno_m} --keep ${1}/mother_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_gcta_alpha-0.25_genoF_pheno${i}_reml --grm ${1}/test_extracted_genoF_ldak_gcta_alpha-0.25_grm --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_gcta_alpha-0.25_hap_pheno${i}_reml --mgrm ${1}/hgcta_ldak_gcta_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 "$iterations"); do
        pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        ldak6.linux --reml ${2}/pheno${i}/sim_ldak_gcta_alpha-0.25_mf_pheno${i}_reml --mgrm ${1}/mgcta_ldak_gcta_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/fetus_grm0.05_extract_list.keep --max-threads 16;
done

#run REML using GRM created by GCTA model

for i in $(seq 1 "$iterations"); do
	pheno_m=$(find ${2}/pheno${i} -type f -name "*_m.phen");
       gcta64 --reml --out ${2}/pheno${i}/sim_genoM_pheno${i}_reml --grm ${1}/test_extracted_genoM_grm0.05 --pheno ${pheno_m} --reml-no-constrain --thread-num 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        gcta64 --reml --out ${2}/pheno${i}/sim_genoF_pheno${i}_reml --grm ${1}/test_extracted_genoF_grm0.05 --pheno ${pheno_f} --reml-no-constrain --thread-num 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        gcta64 --reml --out ${2}/pheno${i}/sim_hap_pheno${i}_reml --mgrm ${1}/hgcta_mduos_grm_list.txt --pheno ${pheno_f} --reml-lrt 1 2 3 --reml-no-constrain --thread-num 16;
done

for i in $(seq 1 "$iterations"); do
	pheno_f=$(find ${2}/pheno${i} -type f -name "*_f.phen");
        gcta64 --reml --out ${2}/pheno${i}/sim_mf_pheno${i}_reml --mgrm ${1}/mgcta_mduos_grm_list.txt --pheno ${pheno_f} --reml-lrt 1 2 3 --reml-no-constrain --thread-num 16;
done
