#test_causal_variants.txt is a list of randomly selected variants (one variant in each line) from test data.

args = commandArgs(trailingOnly = T);

pheno_dir = args[1];
geno_dir = args[2];
phen_num = as.numeric(args[3]);
nvar1 = as.numeric(args[4]);
nvar2 = as.numeric(args[5]);
nvar3 = as.numeric(args[6]);
var_m = as.numeric(args[7]);
var_f = as.numeric(args[8]);

cat(pheno_dir, geno_dir, nvar1, nvar2, nvar3, var_m, var_f, "\n");

current_dir = getwd()
source(file.path(current_dir, "sim_effects_and_pheno.R"));

sim_eff_pheno(direc = pheno_dir, geno_dir = geno_dir, phen_rep = phen_num, var_ids = file.path(current_dir, "test_causal_variants.txt"), nvar = c(nvar1, nvar2, nvar3), var_m = var_m, var_f = var_f);
