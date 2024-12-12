#!/bin/bash

#Lines 5-10 are commented out to run this script on any personal computer.

#BSUB -W 72:00
#BSUB -n 8
#BSUB -M 32768
#BSUB -R "span[ptile=2]"
#BSUB -e /data/predataSamit/empirical_data_analysis/%J.err
#BSUB -o /data/predataSamit/empirical_data_analysis/%J.out

read -p "please provide the working directory with full path " wd
cd $wd

#read phenotype name without extension [phenotype must have three columns in the order FID, IID, and Phenotype (without header)]

read -p "please provide mothers' phenotype file name without extension which must be present in the working directory " pheno_m
read -p "please provide children's phenotype file name without extension which must be present in the working directory " pheno_c

#create a directory within working directory for results

if [ -d "results" ]
then
	echo -e "results will saved in $wd/results.\n"
else
	mkdir -p ./results;
	echo -e "results will saved in $wd/results.\n"
fi

#run REML using GRM created by LDAK-Thin model
#alpha = -1.0

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_thin_alpha-1.0_genoM_20pcs_reml --grm test_extracted_genoM_thin_alpha-1.0_grm --pheno ${pheno_m}.phen --keep ./mother_grm0.05_extract_list.keep --covar test_extracted_genoM_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_thin_alpha-1.0_genoF_20pcs_reml --grm test_extracted_genoF_thin_alpha-1.0_grm --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoF_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_thin_alpha-1.0_hap_30pcs_reml --mgrm hgcta_ldak_thin_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_hap_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_thin_alpha-1.0_mf_20pcs_reml --mgrm mgcta_ldak_thin_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoFM_kinship0.05_pc.eigenvec --max-threads 16;

#alpha = -0.25

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_thin_alpha-0.25_genoM_20pcs_reml --grm test_extracted_genoM_thin_alpha-0.25_grm --pheno ${pheno_m}.phen --keep ./mother_grm0.05_extract_list.keep --covar test_extracted_genoM_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_thin_alpha-0.25_genoF_20pcs_reml --grm test_extracted_genoF_thin_alpha-0.25_grm --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep  --covar test_extracted_genoF_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_thin_alpha-0.25_hap_30pcs_reml --mgrm hgcta_ldak_thin_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_hap_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_thin_alpha-0.25_mf_20pcs_reml --mgrm mgcta_ldak_thin_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoFM_kinship0.05_pc.eigenvec --max-threads 16;

#run REML using GRM created by LDAK-Weights model
#alpha = -1.0

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_weights_alpha-1.0_genoM_20pcs_reml --grm test_extracted_genoM_ldak_weights_alpha-1.0_grm --pheno ${pheno_m}.phen --keep ./mother_grm0.05_extract_list.keep --covar test_extracted_genoM_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_weights_alpha-1.0_genoF_20pcs_reml --grm test_extracted_genoF_ldak_weights_alpha-1.0_grm --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoF_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_weights_alpha-1.0_hap_30pcs_reml --mgrm hgcta_ldak_weights_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_hap_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_weights_alpha-1.0_mf_20pcs_reml --mgrm mgcta_ldak_weights_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoFM_kinship0.05_pc.eigenvec --max-threads 16;

#alpha = -0.25

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_weights_alpha-0.25_genoM_20pcs_reml --grm test_extracted_genoM_ldak_weights_alpha-0.25_grm --pheno ${pheno_m}.phen --keep ./mother_grm0.05_extract_list.keep --covar test_extracted_genoM_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_weights_alpha-0.25_genoF_20pcs_reml --grm test_extracted_genoF_ldak_weights_alpha-0.25_grm --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoF_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_weights_alpha-0.25_hap_30pcs_reml --mgrm hgcta_ldak_weights_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_hap_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_weights_alpha-0.25_mf_20pcs_reml --mgrm mgcta_ldak_weights_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoFM_kinship0.05_pc.eigenvec --max-threads 16;

#run REML using GRM created by LDAK-GCTA model
#alpha = -1.0

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_gcta_alpha-1.0_genoM_20pcs_reml --grm test_extracted_genoM_ldak_gcta_alpha-1.0_grm --pheno ${pheno_m}.phen --keep ./mother_grm0.05_extract_list.keep --covar test_extracted_genoM_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_gcta_alpha-1.0_genoF_20pcs_reml --grm test_extracted_genoF_ldak_gcta_alpha-1.0_grm --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoF_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_gcta_alpha-1.0_hap_30pcs_reml --mgrm hgcta_ldak_gcta_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_hap_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_gcta_alpha-1.0_mf_20pcs_reml --mgrm mgcta_ldak_gcta_mduos_alpha-1.0_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoFM_kinship0.05_pc.eigenvec --max-threads 16;

#alpha = -0.25

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_gcta_alpha-0.25_genoM_20pcs_reml --grm test_extracted_genoM_ldak_gcta_alpha-0.25_grm --pheno ${pheno_m}.phen --keep ./mother_grm0.05_extract_list.keep --covar test_extracted_genoM_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_gcta_alpha-0.25_genoF_20pcs_reml --grm test_extracted_genoF_ldak_gcta_alpha-0.25_grm --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoF_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_gcta_alpha-0.25_hap_30pcs_reml --mgrm hgcta_ldak_gcta_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_hap_kinship0.05_pc.eigenvec --max-threads 16;

ldak6.linux --reml ./results/test_extracted_grm0.05_ldak_gcta_alpha-0.25_mf_20pcs_reml --mgrm mgcta_ldak_gcta_mduos_alpha-0.25_grm_list.txt --pheno ${pheno_c}.phen --keep ./fetus_grm0.05_extract_list.keep --covar test_extracted_genoFM_kinship0.05_pc.eigenvec --max-threads 16;

#run REML using GRM created by GCTA model

gcta64 --reml --out ./results/test_extracted_grm0.05_genoM_20pcs_reml --grm test_extracted_genoM_grm0.05 --pheno ${pheno_m}.phen --qcovar test_extracted_genoM_kinship0.05_pc.eigenvec --threads 16;

gcta64 --reml --out ./results/test_extracted_grm0.05_genoF_20pcs_reml --grm test_extracted_genoF_grm0.05 --pheno ${pheno_c}.phen --qcovar test_extracted_genoF_kinship0.05_pc.eigenvec --threads 16;

gcta64 --reml --out ./results/test_extracted_grm0.05_hap_30pcs_reml --mgrm hgcta_mduos_grm_list.txt --pheno ${pheno_c}.phen --qcovar test_extracted_hap_kinship0.05_pc.eigenvec --threads 16;

gcta64 --reml --out ./results/test_extracted_grm0.05_mf_20pcs_reml --mgrm mgcta_mduos_grm_list.txt --pheno ${pheno_c}.phen --keep fetus_grm0.05_extract_list.txt --qcovar test_extracted_genoFM_kinship0.05_pc.eigenvec --threads 16;


#run REML using Threshold-GRM method

gcta64 --reml --out ./results/test_extracted_grm0.05_threshold_genoM_20pcs_reml --mgrm conventional_threshold_genoM_grm_list.txt --pheno ${pheno_m}.phen --qcovar test_extracted_genoM_kinship0.05_pc.eigenvec --threads 16;

gcta64 --reml --out ./results/test_extracted_grm0.05_threshold_genoF_20pcs_reml --mgrm conventional_threshold_genoF_grm_list.txt --pheno ${pheno_c}.phen --qcovar test_extracted_genoF_kinship0.05_pc.eigenvec --threads 16;

gcta64 --reml --out ./results/test_extracted_grm0.05_threshold_hap_30pcs_reml --mgrm hgcta_threshold_mduos_grm_list.txt --pheno ${pheno_c}.phen --qcovar test_extracted_hap_kinship0.05_pc.eigenvec --threads 16;

gcta64 --reml --out ./results/test_extracted_grm0.05_threshold_mf_20pcs_reml --mgrm mgcta_threshold_mduos_grm_list.txt --pheno ${pheno_c}.phen --qcovar test_extracted_genoFM_kinship0.05_pc.eigenvec --threads 16;

#run HE regression using genetic covariance calculated from GREML i.e. all SNPs have same weights and alpha = -1.0

gcta64 --HEreg --out ./results/test_extracted_grm0.05_genoM_20pcs_HE --grm test_extracted_genoM_grm0.05 --pheno ${pheno_m}.phen --qcovar test_extracted_genoM_kinship0.05_pc.eigenvec --threads 16;

gcta64 --HEreg --out ./results/test_extracted_grm0.05_genoF_20pcs_HE --grm test_extracted_genoF_grm0.05 --pheno ${pheno_c}.phen --qcovar test_extracted_genoF_kinship0.05_pc.eigenvec --threads 16;

gcta64 --HEreg --out ./results/test_extracted_grm0.05_hap_HE --mgrm hgcta_mduos_grm_list.txt --pheno ${pheno_c}.phen --qcovar test_extracted_hap_pc.eigenvec --threads 16;

gcta64 --HEreg --out ./results/test_extracted_grm0.05_mf_HE --mgrm mgcta_mduos_grm_list.txt --pheno ${pheno_c}.phen --qcovar test_extracted_genoFM_kinship0.05_pc.eigenvec --threads 16;
