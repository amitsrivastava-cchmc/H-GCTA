#!/bin/bash

#Line 5-12 are commented out to run the script on any personal computer which has R and Rscript in computer's "$PATH".

#BSUB -W 24:00
#BSUB -n 8
#BSUB -M 32768
#BSUB -R "span[ptile=2]"
#BSUB -e /data/user/simulation_with_poe_and_mf_cor/%J.err
#BSUB -o /data/user/simulation_with_poe_and_mf_cor/%J.out

#module load R

#Lines 28-33 are the most crucial variables which decide what would be the number of iterations and relative maternal and fetal genetic contribution to dyadic traits.

#add directory with binaries to $PATH

binpath=$(realpath ../bin)
export PATH="$PATH:$binpath"

#read directories
#pheno_dir=$(echo "/home/user/H-GCTA/simulation/simulation_mf_with_ppoe_1.0_m1_0.0_mf_cor_1.0");
#geno_dir=$(echo "/home/user/H-GCTA/data.test");

read -p "provide the directory where simulated phenotypes and their results will be saved (full path) " pheno_dir
read -p "provide the directory where plink format files and corresponding GRMs are saved (full path) " geno_dir
current_dir=$(pwd);
iterations=5
nvar1=0;
nvar2=0;
nvar3=2000;
var_m=1;
var_f=1;

for i in $(seq 1 "$iterations"); do
	mkdir -p "$pheno_dir"/pheno${i};
done

# create plink format files with kinship coefficient above set cut-off used in create_grm_extend.sh

plink --bfile ${geno_dir}/test_extracted_genoM --keep ${geno_dir}/mother_grm0.05_extract_list.txt --make-bed --out ${geno_dir}/test_extracted_genoM_kinship0.05;
for gt in $(echo genoF hapM1 hapM2 hapP1); do
	plink --bfile ${geno_dir}/test_extracted_${gt} --keep ${geno_dir}/fetus_grm0.05_extract_list.txt --make-bed --out ${geno_dir}/test_extracted_${gt}_kinship0.05;
done

#created vector of parental, fetal and parent-of-origin effects

Rscript ${current_dir}/sim_effects_and_pheno_run.R ${pheno_dir} ${geno_dir} $iterations $nvar1 $nvar2 $nvar3 $var_m $var_f

#run reml analyses

${current_dir}/run_reml.sh ${geno_dir} ${pheno_dir} $iterations

#process results from GCTA

for i in $(seq 1 "$iterations"); do
	awk '$1 ~ /\/Vp/' ${pheno_dir}/pheno${i}/sim_genoM_pheno${i}_reml.hsq >> ${pheno_dir}/gcta_sim_results_genoM.txt;
done
for i in $(seq 1 "$iterations"); do
	awk '$1 ~ /\/Vp/' ${pheno_dir}/pheno${i}/sim_genoF_pheno${i}_reml.hsq >> ${pheno_dir}/gcta_sim_results_genoF.txt
done
for i in $(seq 1 "$iterations"); do
	awk '$1 ~ /\/Vp/' ${pheno_dir}/pheno${i}/sim_hap_pheno${i}_reml.hsq >> ${pheno_dir}/gcta_sim_results_hap.txt;
done
for i in $(seq 1 "$iterations"); do
	awk '$1 ~ /\/Vp/' ${pheno_dir}/pheno${i}/sim_mf_pheno${i}_reml.hsq >> ${pheno_dir}/gcta_sim_results_mf.txt;
done

awk '$1 ~ /V\(G1\)\/Vp/' ${pheno_dir}/gcta_sim_results_hap.txt > ${pheno_dir}/gcta_sim_results_M1.txt
awk '$1 ~ /V\(G2\)\/Vp/' ${pheno_dir}/gcta_sim_results_hap.txt > ${pheno_dir}/gcta_sim_results_M2.txt
awk '$1 ~ /V\(G3\)\/Vp/' ${pheno_dir}/gcta_sim_results_hap.txt > ${pheno_dir}/gcta_sim_results_P1.txt

awk '$1 ~ /V\(G1\)\/Vp/' ${pheno_dir}/gcta_sim_results_mf.txt > ${pheno_dir}/gcta_sim_results_M_var.txt
awk '$1 ~ /V\(G2\)\/Vp/' ${pheno_dir}/gcta_sim_results_mf.txt > ${pheno_dir}/gcta_sim_results_F_var.txt
awk '$1 ~ /V\(G3\)\/Vp/' ${pheno_dir}/gcta_sim_results_mf.txt > ${pheno_dir}/gcta_sim_results_MF_covar.txt

#process results from LADK-GCTA (alpha = -1.0)

for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_gcta_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_genoM.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_gcta_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_genoF.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_gcta_alpha-1.0_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_hap.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_gcta_alpha-1.0_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_MF_covar.txt

#process results from LADK-GCTA (alpha = -0.25)

for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_gcta_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_genoM.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_gcta_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_genoF.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_gcta_alpha-0.25_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_hap.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_gcta_alpha-0.25_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_MF_covar.txt

#process results from LADK-Thin (alpha = -1.0)

for i in $(seq 1 "$iterations"); do
	awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_thin_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_genoM.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_thin_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_genoF.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_thin_alpha-1.0_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_hap.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_thin_alpha-1.0_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_MF_covar.txt

#process results from LADK-Thin (alpha = -0.25)

for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_thin_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_genoM.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_thin_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_genoF.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_thin_alpha-0.25_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_hap.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_thin_alpha-0.25_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_MF_covar.txt

#process results from LADK-Weights (alpha = -1.0)

for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_weights_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_genoM.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_weights_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_genoF.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_weights_alpha-1.0_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_hap.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_weights_alpha-1.0_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_MF_covar.txt

#process results from LADK-Weights (alpha = -0.25)

for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_weights_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_genoM.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_ldak_weights_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_genoF.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_weights_alpha-0.25_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_hap.txt;
done
for i in $(seq 1 "$iterations"); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_ldak_weights_alpha-0.25_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_MF_covar.txt

#process results from simulated data to create tables with h2, SE, and p- values

Rscript ${current_dir}/process_sim_results.R $pheno_dir