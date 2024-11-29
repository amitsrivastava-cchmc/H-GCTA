#!/bin/bash

#commands in line 5-10 are specific to HPC (IBM-LSF). To run on any machine please comment out/remove commands in line 5-10.

#BSUB -W 24:00
#BSUB -n 8
#BSUB -M 32768
#BSUB -R "span[ptile=2]"
#BSUB -e /data/predataSamit/simulation_with_poe_and_mf_cor/%J.err
#BSUB -o /data/predataSamit/simulation_with_poe_and_mf_cor/%J.out

pheno_dir=$(echo "/data/predataSamit/simulation_with_poe_and_mf_cor/simulation_mat_pooled_data")
grm_dir=$(echo "/data/predataSamit/empirical_data_analysis/pooled_data_maf0.000001")
nvar1=10000
nvar2=0
nvar3=0
var_m=1
var_f=1

module load R

Rscript /data/predataSamit/simulation_with_poe_and_mf_cor/sim_effects_and_pheno_pooled_data_run.R $pheno_dir $nvar1 $nvar2 $nvar3 $var_m $var_f

/data/predataSamit/simulation_with_poe_and_mf_cor/run_reml_pooled_data.sh ${grm_dir} ${pheno_dir}

#process results from GCTA

for i in $(seq 1 100); do
	awk '$1 ~ /\/Vp/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_genoM_pheno${i}_reml.hsq >> ${pheno_dir}/gcta_sim_results_genoM.txt;
done
for i in $(seq 1 100); do
	awk '$1 ~ /\/Vp/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_genoF_pheno${i}_reml.hsq >> ${pheno_dir}/gcta_sim_results_genoF.txt
done
for i in $(seq 1 100); do
	awk '$1 ~ /\/Vp/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_hap_pheno${i}_reml.hsq >> ${pheno_dir}/gcta_sim_results_hap.txt;
done
for i in $(seq 1 100); do
	awk '$1 ~ /\/Vp/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_mf_pheno${i}_reml.hsq >> ${pheno_dir}/gcta_sim_results_mf.txt;
done

awk '$1 ~ /V\(G1\)\/Vp/' ${pheno_dir}/gcta_sim_results_hap.txt > ${pheno_dir}/gcta_sim_results_M1.txt
awk '$1 ~ /V\(G2\)\/Vp/' ${pheno_dir}/gcta_sim_results_hap.txt > ${pheno_dir}/gcta_sim_results_M2.txt
awk '$1 ~ /V\(G3\)\/Vp/' ${pheno_dir}/gcta_sim_results_hap.txt > ${pheno_dir}/gcta_sim_results_P1.txt

awk '$1 ~ /V\(G1\)\/Vp/' ${pheno_dir}/gcta_sim_results_mf.txt > ${pheno_dir}/gcta_sim_results_M_var.txt
awk '$1 ~ /V\(G2\)\/Vp/' ${pheno_dir}/gcta_sim_results_mf.txt > ${pheno_dir}/gcta_sim_results_F_var.txt
awk '$1 ~ /V\(G3\)\/Vp/' ${pheno_dir}/gcta_sim_results_mf.txt > ${pheno_dir}/gcta_sim_results_MF_covar.txt

#process results from LADK-GCTA (alpha = -1.0)

for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_genoM.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_genoF.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-1.0_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_hap.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-1.0_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-1.0_sim_results_MF_covar.txt

#process results from LADK-GCTA (alpha = -0.25)

for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_genoM.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_genoF.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-0.25_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_hap.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_gcta_alpha-0.25_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_gcta_alpha-0.25_sim_results_MF_covar.txt

#process results from LADK-Thin (alpha = -1.0)

for i in $(seq 1 100); do
	awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_genoM.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_genoF.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-1.0_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_hap.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-1.0_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-1.0_sim_results_MF_covar.txt

#process results from LADK-Thin (alpha = -0.25)

for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_genoM.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_genoF.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-0.25_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_hap.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_thin_alpha-0.25_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_thin_alpha-0.25_sim_results_MF_covar.txt

#process results from LADK-Weights (alpha = -1.0)

for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_genoM.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-1.0_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_genoF.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-1.0_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_hap.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-1.0_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-1.0_sim_results_MF_covar.txt

#process results from LADK-Weights (alpha = -0.25)

for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_genoM.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K1/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-0.25_genoM_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_genoF.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-0.25_hap_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_hap.txt;
done
for i in $(seq 1 100); do
        awk '$1 ~ /Her_K/' ${pheno_dir}/pheno${i}/sim_comm_grm0.05_ldak_weights_alpha-0.25_mf_pheno${i}_reml.reml >> ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_mf.txt;
done

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_M1.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_M2.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_hap.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_P1.txt

awk '$1 ~ /Her_K1/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_M_var.txt
awk '$1 ~ /Her_K2/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_F_var.txt
awk '$1 ~ /Her_K3/' ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_mf.txt > ${pheno_dir}/ldak_weights_alpha-0.25_sim_results_MF_covar.txt


Rscript /data/predataSamit/simulation_with_poe_and_mf_cor/process_sim_results.R $pheno_dir
