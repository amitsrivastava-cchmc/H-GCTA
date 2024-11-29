#!/bin/bash

#${1} = directory where grm is saved e.g. "/data/predataSamit/empirical_data_analysis/pooled_data_maf0.000001"
#${2} = directory where pheotypes are present e.g "/data/predataSamit/simulation_with_poe_and_mf_cor/simulation_mf_with_ppoe_1.0_m1_1.0_mf_cor_0"
#outputs will also be written in ${2}

wd=$(pwd);
cd ${2};
echo $(pwd);

#pheno_m=$(find ./pheno1 -type f -name "*_m.phen");
#pheno_f=$(find ./pheno1 -type f -name "*_f.phen");

#execute program

#run REML using GRM created by LDAK-Thin model

#using alpha = -1.0

for i in $(seq 1 100); do
	pheno_m=$(find ./pheno${i} -type f -name "*_m.phen");
	/data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-1.0_genoM_pheno${i}_reml --grm ${1}/pooled_data_genoM_comm_0.000001_thin_alpha-1.0_grm --pheno ${pheno_m} --keep ${1}/pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-1.0_genoF_pheno${i}_reml --grm ${1}/pooled_data_genoF_comm_0.000001_thin_alpha-1.0_grm --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-1.0_hap_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/hgcta_ldak_thin_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-1.0_mf_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/mgcta_ldak_thin_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

#using alpha = -0.25

for i in $(seq 1 100); do
        pheno_m=$(find ./pheno${i} -type f -name "*_m.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-0.25_genoM_pheno${i}_reml --grm ${1}/pooled_data_genoM_comm_0.000001_thin_alpha-0.25_grm --pheno ${pheno_m} --keep ${1}/pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-0.25_genoF_pheno${i}_reml --grm ${1}/pooled_data_genoF_comm_0.000001_thin_alpha-0.25_grm --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-0.25_hap_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/hgcta_ldak_thin_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-0.25_mf_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/mgcta_ldak_thin_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

#run REML using GRM created by LDAK-Weights model

#using alpha = -1.0

for i in $(seq 1 100); do
	pheno_m=$(find ./pheno${i} -type f -name "*_m.phen");
	/data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-1.0_genoM_pheno${i}_reml --grm ${1}/pooled_data_genoM_comm_0.000001_ldak_weights_alpha-1.0_grm --pheno ${pheno_m} --keep ${1}/pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-1.0_genoF_pheno${i}_reml --grm ${1}/pooled_data_genoF_comm_0.000001_ldak_weights_alpha-1.0_grm --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
	/data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-1.0_hap_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/hgcta_ldak_weights_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-1.0_mf_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/mgcta_ldak_weights_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

#using alpha = -.025

for i in $(seq 1 100); do
       pheno_m=$(find ./pheno${i} -type f -name "*_m.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-0.25_genoM_pheno${i}_reml --grm ${1}/pooled_data_genoM_comm_0.000001_ldak_weights_alpha-0.25_grm --pheno ${pheno_m} --keep ${1}/pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-0.25_genoF_pheno${i}_reml --grm ${1}/pooled_data_genoF_comm_0.000001_ldak_weights_alpha-0.25_grm --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-0.25_hap_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/hgcta_ldak_weights_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-0.25_mf_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/mgcta_ldak_weights_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

#run REML using GRM created by LDAK-GCTA model

#using alpha = -1.0

for i in $(seq 1 100); do
	pheno_m=$(find ./pheno${i} -type f -name "*_m.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-1.0_genoM_pheno${i}_reml --grm ${1}/pooled_data_genoM_comm_0.000001_ldak_gcta_alpha-1.0_grm --pheno ${pheno_m} --keep ${1}/pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-1.0_genoF_pheno${i}_reml --grm ${1}/pooled_data_genoF_comm_0.000001_ldak_gcta_alpha-1.0_grm --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-1.0_hap_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/hgcta_ldak_gcta_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-1.0_mf_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/mgcta_ldak_gcta_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

#using alpha = -0.25

for i in $(seq 1 100); do
        pheno_m=$(find ./pheno${i} -type f -name "*_m.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-0.25_genoM_pheno${i}_reml --grm ${1}/pooled_data_genoM_comm_0.000001_ldak_gcta_alpha-0.25_grm --pheno ${pheno_m} --keep ${1}/pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-0.25_genoF_pheno${i}_reml --grm ${1}/pooled_data_genoF_comm_0.000001_ldak_gcta_alpha-0.25_grm --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-0.25_hap_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/hgcta_ldak_gcta_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

for i in $(seq 1 100); do
        pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/ldak5.1.linux --reml ./pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-0.25_mf_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/mgcta_ldak_gcta_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_f} --keep ${1}/pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --max-threads 16;
done

#run REML using GRM created by GCTA model

for i in $(seq 1 100); do
	pheno_m=$(find ./pheno${i} -type f -name "*_m.phen");
       /data/predataSamit/simulation_with_poe_and_mf_cor/gcta64 --reml --out ./pheno${i}/sim_comm_grm0.05_genoM_pheno${i}_reml --grm ${1}/pooled_data_genoM_comm_0.000001_grm0.05 --pheno ${pheno_m} --reml-no-constrain --thread-num 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/gcta64 --reml --out ./pheno${i}/sim_comm_grm0.05_genoF_pheno${i}_reml --grm ${1}/pooled_data_genoF_comm_0.000001_grm0.05 --pheno ${pheno_f} --reml-no-constrain --thread-num 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/gcta64 --reml --out ./pheno${i}/sim_comm_grm0.05_hap_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/hgcta_mduos_grm_list.txt --pheno ${pheno_f} --reml-lrt 1 2 3 --reml-no-constrain --thread-num 16;
done

for i in $(seq 1 100); do
	pheno_f=$(find ./pheno${i} -type f -name "*_f.phen");
        /data/predataSamit/simulation_with_poe_and_mf_cor/gcta64 --reml --out ./pheno${i}/sim_comm_grm0.05_mf_pheno${i}_reml --mgrm ${wd}/pooled_data_maf0.000001_grm_lists/mgcta_mduos_grm_list.txt --pheno ${pheno_f} --reml-lrt 1 2 3 --reml-no-constrain --thread-num 16;
done

cd $wd;
