#pooled_data_comm_maf0.000001_variants.txt is a list of randomly selected variants (one variant in line) from pooled data with all polymorphic SNPs.

args = commandArgs(trailingOnly = T);

dir_name = args[1]; 
nvar1 = as.numeric(args[2]);
nvar2 = as.numeric(args[3]);
nvar3 = as.numeric(args[4]);
var_m = as.numeric(args[5]);
var_f = as.numeric(args[6]);

cat(dir_name, nvar1, nvar2, nvar3, var_m, var_f, "\n");

source("/data/predataSamit/simulation_with_poe_and_mf_cor/sim_effects_and_pheno_pooled_data.R");

sim_eff_pheno(direc = dir_name, phen_rep = 100, var_ids = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_comm_maf0.000001_variants.txt", nvar = c(nvar1, nvar2, nvar3), var_m = var_m, var_f = var_f);

