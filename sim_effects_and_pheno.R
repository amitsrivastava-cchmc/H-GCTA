#Usage
#Rscript sim_effects_and_pheno_run.R ${pheno_dir} ${geno_dir} $nvar1 $nvar2 $nvar3 $var_m $var_f
#dir_name, nvar1, nvar2, nvar3, var_m and var_f are defined variables in a bash script (run_simulation.sh) in which above Rscript will be executed.

#sim_effects_and_pheno_run.R is as following...

#args = commandArgs(trailingOnly = T);

#pheno_dir = args[1];
#geno_dir = args[2];
#phen_rep = as.numeric(args[3]); # number of iterations
#nvar1 = as.numeric(args[4]); # number of variants with explicit maternal effects
#nvar2 = as.numeric(args[5]); # number of variants with explicit fetal effects
#nvar3 = as.numeric(args[6]); # number of variants with joint maternal-fetal effects
#var_m = as.numweic(args[7]) #variance of maternal effects
#var_f = as.numeric(args[8]) #variance of fetal effects

#source("./sim_effects_and_pheno.R");

#sim_eff_pheno(...);

## simulate effect sizes and phenotypes

sim_eff_pheno = function(direc, geno_dir, phen_rep, var_ids, nvar, var_m, var_f){
	current_dir = getwd();
	source(file.path(current_dir, "create_effects.R"));

	path = file.path(direc, "pheno1");
	base_name = basename(direc);
	zz = unlist(strsplit(base_name, "_"))
	if("mf" %in% zz){
		if("ppoe" %in% zz){
#			sub(".*?ppoe_(.*?)_.*", "\\1", base_name);
#			sub(".*?m1_(.*?)_.*", "\\1", base_name);
			p = which(zz %in% "ppoe");
			ppoe = zz[p +1];
			ppoe = as.numeric(ppoe);
			q = which(zz %in% "m1");
			rho = zz[q + 1];
			rho = as.numeric(rho);
		}

		if("cor" %in% zz){
#			sub(".*?cor_(.*?)_.*", "\\1", base_name);
			r = which(zz %in% "cor");
			correl = zz[r + 1];
			correl = as.numeric(correl);

			if("ppoe" %in% zz){

				create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  correl, POE = T, ppoe = ppoe, rho = rho);

				cat("Creating phenotypes using haplotypes in", path, "\n");

				system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_poes.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_fet_m1  --thread-num 16", geno_dir, direc, direc), intern = T);
			}else{
				cat("Proportion of causal variants with POEs in fetus is missing, cannot incorporate POEs\n")

				create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  correl, POE = F);

				cat("Creating phenotypes using haplotypes in", path, "\n");

				syste(command = sprintf("gcta64  --bfile %s/test_extracted_hapM1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_fet_m1  --thread-num 16", geno_dir, direc, direc), intern = T);
			}
		}else{
			if("ppoe" %in% zz){

				create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  NULL, POE = T, ppoe = ppoe, rho = rho);

				cat("Creating phenotypes using haplotypes in", path, "\n");

				system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_poes.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_fet_m1  --thread-num 16", geno_dir, direc, direc), intern = T);
			}else{
				cat("Proportion of causal variants with POEs in fetus is missing, cannot incorporate POEs\n")

				create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  NULL, POE = F);

				cat("Creating phenotypes using haplotypes in", path, "\n");

				system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_fet_m1  --thread-num 16", geno_dir, direc, direc), intern = T);
			}
		}

		system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_par_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_par_m1  --thread-num 16", geno_dir, direc, direc), intern = T);
		system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM2_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_par_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_par_m2  --thread-num 16", geno_dir, direc, direc), intern = T);
		system(command = sprintf("gcta64  --bfile %s/test_extracted_hapP1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_fet_p1  --thread-num 16", geno_dir, direc, direc), intern = T);

		sim_fet_m1 = read.table(file.path(direc, "pheno1/sim_fet_m1.phen"), head=F, as.is=T);
		sim_fet_p1 = read.table(file.path(direc, "pheno1/sim_fet_p1.phen"), head=F, as.is=T);
		sim_par_m1 = read.table(file.path(direc, "pheno1/sim_par_m1.phen"), head=F, as.is=T);
		sim_par_m2 = read.table(file.path(direc, "pheno1/sim_par_m2.phen"), head=F, as.is=T);
		iter = phen_rep;
		m = sapply(3:(iter+2), function(x) (sim_par_m1[, x] + sim_par_m2[, x] + sim_fet_m1[, x] + sim_fet_p1[, x]))
		mid = sapply(1:2, function(x) sub("A", "M", sim_fet_m1[, x]));
		mid = as.data.frame(mid);
		cid = sim_fet_m1[, 1:2];
		pheno_m = cbind(mid, m);
		pheno_f = cbind(cid, m);
		x = 1;
		while(x <= iter){
			pheno = paste0("pheno", x);
			if("cor" %in% zz){
				if("ppoe" %in% zz){
					pheno_out_m = paste("sim_mf_ppoe", ppoe, "m1", rho, "mf_cor", correl, "m.phen", sep="_");
					pheno_out_f = paste("sim_mf_ppoe", ppoe, "m1", rho, "mf_cor", correl, "f.phen", sep="_");
				}else{
					pheno_out_m = paste("sim_mf", "cor", correl, "m.phen", sep="_");
					pheno_out_f = paste("sim_mf", "cor", correl, "f.phen", sep="_");
				}
			}else{
				if("ppoe" %in% zz){
					pheno_out_m = paste("sim_mf_ppoe", ppoe, "m1", rho, "m.phen", sep="_");
					pheno_out_f = paste("sim_mf_ppoe", ppoe, "m1", rho, "f.phen", sep="_");
				}else{
					pheno_out_m = paste("sim_mf", "m.phen", sep="_");
					pheno_out_f = paste("sim_mf", "f.phen", sep="_");
				}
			}
			path_m = file.path(direc, paste0("pheno", x), pheno_out_m);
			path_f = file.path(direc, paste0("pheno", x), pheno_out_f);
			write.table(pheno_m[, c(1,2, (x+2))], path_m, sep="\t", row.names=F, col.names=F, quote=F);
			write.table(pheno_f[, c(1,2, (x+2))], path_f, sep="\t", row.names=F, col.names=F, quote=F);

			cat(pheno, ":", "Final phenotype is created in", file.path(direc, pheno), "\n");
			x = x+1;
		}
	}else if("fet" %in% zz){
		if("cor" %in% zz){
			message("Matrix of variance-covariance cannot be used in the absence of joint parental-fetal effects\n");
		}

		if("ppoe" %in% zz){
#			sub(".*?ppoe_(.*?)_.*", "\\1", base_name);
#			sub(".*?m1_(.*?)_.*", "\\1", base_name);
			p = which(zz %in% "ppoe");
			ppoe = zz[p +1];
			ppoe = as.numeric(ppoe);
			q = which(zz %in% "m1");
			rho = zz[q + 1];
			rho = as.numeric(rho);

			create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  NULL, POE = T, ppoe = ppoe, rho = rho);

			cat("Creating phenotypes using haplotypes in", path, "\n");

			system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_poes.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_fet_m1  --thread-num 16", geno_dir, direc, direc), intern = T);
		}else{

			create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  NULL, POE = F);

			system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_fet_m1  --thread-num 16", geno_dir, direc, direc), intern = T);
		}

		system(command = sprintf("gcta64  --bfile %s/test_extracted_hapP1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_fet_p1  --thread-num 16", geno_dir, direc, direc), intern = T);

		sim_fet_m1 = read.table(file.path(direc, "pheno1/sim_fet_m1.phen"), head=F, as.is=T);
		sim_fet_p1 = read.table(file.path(direc, "pheno1/sim_fet_p1.phen"), head=F, as.is=T);

		iter = phen_rep;
		m = sapply(3:(iter+2), function(x) (sim_fet_m1[, x] + sim_fet_p1[, x]))
		mid = sapply(1:2, function(x) sub("A", "M", sim_fet_m1[, x]));
		mid = as.data.frame(mid);
		cid = sim_fet_m1[, 1:2];
		pheno_m = cbind(mid, m);
		pheno_f = cbind(cid, m);
		x = 1;
		while(x <= iter){
			pheno = paste0("pheno", x);
			if("ppoe" %in% zz){
				pheno_out_m = paste("sim_fet_ppoe", ppoe, "m1", rho, "m.phen", sep="_");
				pheno_out_f = paste("sim_fet_ppoe", ppoe, "m1", rho, "f.phen", sep="_");
			}else{
				pheno_out_m = paste("sim_fet", "m.phen", sep="_");
				pheno_out_f = paste("sim_fet", "f.phen", sep="_");
			}
			path_m = file.path(direc, paste0("pheno", x), pheno_out_m);
			path_f = file.path(direc, paste0("pheno", x), pheno_out_f);
			write.table(pheno_m[, c(1,2, (x+2))], path_m, sep="\t", row.names=F, col.names=F, quote=F);
			write.table(pheno_f[, c(1,2, (x+2))], path_f, sep="\t", row.names=F, col.names=F, quote=F);

			cat(pheno, ":", "Final phenotype is created in", file.path(direc, pheno), "\n");
			x = x+1;
		}
	}else if("mat" %in% zz){

		create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  NULL, POE = F);

		system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM1_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_par_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_par_m1  --thread-num 16", geno_dir, direc, direc), intern = T);
		system(command = sprintf("gcta64  --bfile %s/test_extracted_hapM2_kinship0.05  --simu-qt  --simu-causal-loci %s/pheno1/sim_par_effects.txt --simu-hsq 0.5 --simu-rep 100  --out %s/pheno1/sim_par_m2  --thread-num 16", geno_dir, direc, direc), intern = T);

		sim_par_m1 = read.table(file.path(direc, "pheno1/sim_par_m1.phen"), head=F, as.is=T);
		sim_par_m2 = read.table(file.path(direc, "pheno1/sim_par_m2.phen"), head=F, as.is=T);

		iter = phen_rep;
		m = sapply(3:(iter+2), function(x) (sim_par_m1[, x] + sim_par_m2[, x]))
		mid = sapply(1:2, function(x) sub("A", "M", sim_par_m1[, x]));
		mid = as.data.frame(mid);
		cid = sim_par_m1[, 1:2];
		pheno_m = cbind(mid, m);
		pheno_f = cbind(cid, m);
		x = 1;
		while(x <= iter){
			pheno = paste0("pheno", x);
			pheno_out_m = paste("sim_mat", "m.phen", sep="_");
			pheno_out_f = paste("sim_mat", "f.phen", sep="_");
			path_m = file.path(direc, paste0("pheno", x), pheno_out_m);
			path_f = file.path(direc, paste0("pheno", x), pheno_out_f);
			write.table(pheno_m[, c(1,2, (x+2))], path_m, sep="\t", row.names=F, col.names=F, quote=F);
			write.table(pheno_f[, c(1,2, (x+2))], path_f, sep="\t", row.names=F, col.names=F, quote=F);

			cat(pheno, ":", "Final phenotype is created in", file.path(direc, pheno), "\n");
			x = x+1;
		}
	}
}
