#Use => module load R; Rscript sim_effects_and_pheno_run.R $dir_name $nvar1 $nvar2 $nvar3 var_m var_f
##dir_name, nvar1, nvar2, nvar3, var_m and var_f are defined variables in a bash script in which above Rscript will be run.

#sim_effects_and_pheno_run.R is as following...

#args = commandArgs(trailingOnly = T);

#dir_name = args[1];
#nvar1 = as.numeric(args[2]);
#nvar2 = as.numeric(args[3]);
#nvar3 = as.numeric(args[4]);
#var_m = variance of maternal effects
#var_f = variance of fetal effects
#source("/data/predataSamit/simulation_with_poe_and_mf_cor/sim_effects_and_pheno.R");
#sim_eff_pheno(...);

## simulate effect sizes and phenotypes

sim_eff_pheno = function(direc, phen_rep, var_ids, nvar, var_m, var_f){
	setwd(direc);
	source("/data/predataSamit/simulation_with_poe_and_mf_cor/create_effects.R");

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

				system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_poes.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_fet_m1  --thread-num 16");
			}else{
				cat("Proportion of causal variants with POEs in fetus is missing, cannot incorporate POEs\n")

				create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  correl, POE = F);

				cat("Creating phenotypes using haplotypes in", path, "\n");

				system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_fet_m1  --thread-num 16");
			}
		}else{
			if("ppoe" %in% zz){

				create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  NULL, POE = T, ppoe = ppoe, rho = rho);

				cat("Creating phenotypes using haplotypes in", path, "\n");

				system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_poes.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_fet_m1  --thread-num 16");
			}else{
				cat("Proportion of causal variants with POEs in fetus is missing, cannot incorporate POEs\n")

				create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  NULL, POE = F);

				cat("Creating phenotypes using haplotypes in", path, "\n");

				system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_fet_m1  --thread-num 16");
			}
		}

		system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_par_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_par_m1  --thread-num 16");
		system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM2_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_par_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_par_m2  --thread-num 16");
		system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapP1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_fet_p1  --thread-num 16");

		sim_fet_m1 = read.table("./pheno1/sim_fet_m1.phen", head=F, as.is=T);
		sim_fet_p1 = read.table("./pheno1/sim_fet_p1.phen", head=F, as.is=T);
		sim_par_m1 = read.table("./pheno1/sim_par_m1.phen", head=F, as.is=T);
		sim_par_m2 = read.table("./pheno1/sim_par_m2.phen", head=F, as.is=T);
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

			system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_poes.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_fet_m1  --thread-num 16");
		}else{

			create_effects(direc = direc, iter = 1, var_ids = var_ids, nvar = nvar, var_m = var_m, var_f = var_f, covar =  NULL, POE = F);

			system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_fet_m1  --thread-num 16");
		}

		system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapP1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_fet_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_fet_p1  --thread-num 16");

		sim_fet_m1 = read.table("./pheno1/sim_fet_m1.phen", head=F, as.is=T);
		sim_fet_p1 = read.table("./pheno1/sim_fet_p1.phen", head=F, as.is=T);

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

		system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM1_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_par_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_par_m1  --thread-num 16");
		system(command = "/data/predataSamit/simulation_with_poe_and_mf_cor/gcta64  --bfile /data/predataSamit/simulation_with_poe_and_mf_cor/gcta_simulation_pooled_data_maf0.000001/pooled_data_hapM2_comm_0.000001_kinship0.05  --simu-qt  --simu-causal-loci ./pheno1/sim_par_effects.txt --simu-hsq 0.5 --simu-rep 100  --out ./pheno1/sim_par_m2  --thread-num 16");

		sim_par_m1 = read.table("./pheno1/sim_par_m1.phen", head=F, as.is=T);
		sim_par_m2 = read.table("./pheno1/sim_par_m2.phen", head=F, as.is=T);

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

