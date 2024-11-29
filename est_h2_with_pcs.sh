#!/bin/bash

#commands in line 5-10 are specific to HPC (IBM-LSF). To run it on any machine, please comment out/remove commands in line 5-10.

#BSUB -W 72:00
#BSUB -n 8
#BSUB -M 32768
#BSUB -R "span[ptile=2]"
#BSUB -e /data/predataSamit/empirical_data_analysis/pooled_data_maf0.000001/%J.err
#BSUB -o /data/predataSamit/empirical_data_analysis/pooled_data_maf0.000001/%J.out

cd /data/predataSamit/empirical_data_analysis/pooled_data_maf0.000001/

#run REML using GRM created by LDAK-Thin model

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_thin_alpha-1.0_genoM_${var}_20pcs_reml --grm pooled_data_genoM_comm_0.000001_thin_alpha-1.0_grm --pheno pooled_data_M_${var}.phen --keep ./pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_thin_alpha-1.0_genoF_${var}_20pcs_reml --grm pooled_data_genoF_comm_0.000001_thin_alpha-1.0_grm --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_thin_alpha-1.0_hap_${var}_30pcs_reml --mgrm hgcta_ldak_thin_mduos_alpha-1.0_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_hap_comm_0.000001_kinship0.05_pc-1.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_thin_alpha-1.0_mf_${var}_20pcs_reml --mgrm mgcta_ldak_thin_mduos_alpha-1.0_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_thin_alpha-0.25_genoM_${var}_20pcs_reml --grm pooled_data_genoM_comm_0.000001_thin_alpha-0.25_grm --pheno pooled_data_M_${var}.phen --keep ./pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_thin_alpha-0.25_genoF_${var}_20pcs_reml --grm pooled_data_genoF_comm_0.000001_thin_alpha-0.25_grm --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep  --covar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_thin_alpha-0.25_hap_${var}_30pcs_reml --mgrm hgcta_ldak_thin_mduos_alpha-0.25_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_hap_comm_0.000001_kinship0.05_pc-1.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_thin_alpha-0.25_mf_${var}_20pcs_reml --mgrm mgcta_ldak_thin_mduos_alpha-0.25_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

#run REML using GRM created by LDAK-Weights model

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_weights_alpha-1.0_genoM_${var}_20pcs_reml --grm pooled_data_genoM_comm_0.000001_ldak_weights_alpha-1.0_grm --pheno pooled_data_M_${var}.phen --keep ./pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_weights_alpha-1.0_genoF_${var}_20pcs_reml --grm pooled_data_genoF_comm_0.000001_ldak_weights_alpha-1.0_grm --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_weights_alpha-1.0_hap_${var}_30pcs_reml --mgrm hgcta_ldak_weights_mduos_alpha-1.0_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_hap_comm_0.000001_kinship0.05_pc-1.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_weights_alpha-1.0_mf_${var}_20pcs_reml --mgrm mgcta_ldak_weights_mduos_alpha-1.0_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_weights_alpha-0.25_genoM_${var}_20pcs_reml --grm pooled_data_genoM_comm_0.000001_ldak_weights_alpha-0.25_grm --pheno pooled_data_M_${var}.phen --keep ./pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_weights_alpha-0.25_genoF_${var}_20pcs_reml --grm pooled_data_genoF_comm_0.000001_ldak_weights_alpha-0.25_grm --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_weights_alpha-0.25_hap_${var}_30pcs_reml --mgrm hgcta_ldak_weights_mduos_alpha-0.25_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_hap_comm_0.000001_kinship0.05_pc-1.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_weights_alpha-0.25_mf_${var}_20pcs_reml --mgrm mgcta_ldak_weights_mduos_alpha-0.25_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

#run REML using GRM created by LDAK-GCTA model

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_gcta_alpha-1.0_genoM_${var}_20pcs_reml --grm pooled_data_genoM_comm_0.000001_ldak_gcta_alpha-1.0_grm --pheno pooled_data_M_${var}.phen --keep ./pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_gcta_alpha-1.0_genoF_${var}_20pcs_reml --grm pooled_data_genoF_comm_0.000001_ldak_gcta_alpha-1.0_grm --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_gcta_alpha-1.0_hap_${var}_30pcs_reml --mgrm hgcta_ldak_gcta_mduos_alpha-1.0_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_hap_comm_0.000001_kinship0.05_pc-1.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_gcta_alpha-1.0_mf_${var}_20pcs_reml --mgrm mgcta_ldak_gcta_mduos_alpha-1.0_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_gcta_alpha-0.25_genoM_${var}_20pcs_reml --grm pooled_data_genoM_comm_0.000001_ldak_gcta_alpha-0.25_grm --pheno pooled_data_M_${var}.phen --keep ./pooled_mother_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_gcta_alpha-0.25_genoF_${var}_20pcs_reml --grm pooled_data_genoF_comm_0.000001_ldak_gcta_alpha-0.25_grm --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_gcta_alpha-0.25_hap_${var}_30pcs_reml --mgrm hgcta_ldak_gcta_mduos_alpha-0.25_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_hap_comm_0.000001_kinship0.05_pc-1.eigenvec --constrain YES --max-threads 16;
done

for var in $(echo gday bw bl hc); do
	../ldak5.1.linux --reml ./pooled_data_comm_grm0.05_ldak_gcta_alpha-0.25_mf_${var}_20pcs_reml --mgrm mgcta_ldak_gcta_mduos_alpha-0.25_grm_list.txt --pheno pooled_data_F_${var}.phen --keep ./pooled_fetus_comm_maf0.000001_grm0.05_extract_list.keep --covar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --constrain YES --max-threads 16;
done

#run REML using GRM created by GCTA model

for var in $(echo gday bw bl hc); do
	./gcta64 --reml --out ./pooled_data_comm_grm0.05_genoM_${var}_20pcs_reml --grm pooled_data_genoM_comm_0.000001_grm0.05 --pheno pooled_data_M_${var}.phen --qcovar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
	./gcta64 --reml --out ./pooled_data_comm_grm0.05_genoF_${var}_20pcs_reml --grm pooled_data_genoF_comm_0.000001_grm0.05 --pheno pooled_data_F_${var}.phen --qcovar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
	./gcta64 --reml --out ./pooled_data_comm_grm0.05_hap_${var}_30pcs_reml --mgrm hgcta_mduos_grm_list.txt --pheno pooled_data_F_${var}.phen --qcovar pooled_data_hap_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
	./gcta64 --reml --out ./pooled_data_comm_grm0.05_mf_${var}_20pcs_reml --mgrm mgcta_mduos_grm_list.txt --pheno pooled_data_F_${var}.phen --qcovar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

#run REML using Threshold-GRM method

for var in $(echo gday bw bl hc); do
        ./gcta64 --reml --out ./pooled_data_comm_grm0.05_threshold_genoM_${var}_20pcs_reml --mgrm conventional_threshold_genoM_grm_list.txt --pheno pooled_data_M_${var}.phen --qcovar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
	./gcta64 --reml --out ./pooled_data_comm_grm0.05_threshold_genoF_${var}_20pcs_reml --mgrm conventional_threshold_genoF_grm_list.txt --pheno pooled_data_F_${var}.phen --qcovar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
	./gcta64 --reml --out ./pooled_data_comm_grm0.05_threshold_hap_${var}_30pcs_reml --mgrm hgcta_threshold_mduos_grm_list.txt --pheno pooled_data_F_${var}.phen --qcovar pooled_data_hap_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
	./gcta64 --reml --out ./pooled_data_comm_grm0.05_threshold_mf_${var}_20pcs_reml --mgrm mgcta_threshold_mduos_grm_list.txt --pheno pooled_data_F_${var}.phen --qcovar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

#run HE regression using genetic covariance calculated from GREML i.e. all SNPs have same weights and alpha = -1.0

for var in $(echo gday bw bl hc); do
       ./gcta64 --HEreg --out ./pooled_data_comm_grm0.05_genoM_${var}_20pcs_HE --grm pooled_data_genoM_comm_0.000001_grm0.05 --pheno pooled_data_M_${var}.phen --qcovar pooled_data_genoM_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
       ./gcta64 --HEreg --out ./pooled_data_comm_grm0.05_genoF_${var}_20pcs_HE --grm pooled_data_genoF_comm_0.000001_grm0.05 --pheno pooled_data_F_${var}.phen --qcovar pooled_data_genoF_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
       ./gcta64 --HEreg --out ./pooled_data_comm_grm0.05_hap_${var}_HE --mgrm hgcta_mduos_grm_list.txt --pheno pooled_data_F_${var}.phen --qcovar pooled_data_hap_comm_0.000001_pc.eigenvec --threads 16;
done

for var in $(echo gday bw bl hc); do
       ./gcta64 --HEreg --out ./pooled_data_comm_grm0.05_mf_${var}_HE --mgrm mgcta_mduos_grm_list.txt --pheno pooled_data_F_${var}.phen --qcovar pooled_data_genoFM_comm_0.000001_kinship0.05_pc.eigenvec --threads 16;
done

