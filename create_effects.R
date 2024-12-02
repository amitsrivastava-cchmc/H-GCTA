# Use => create_effects(direc = "dir-name-with-path" (pheno{1..iter} sub-directories must be present in the given directory),
#						iter = # of simulations,
#						var_ids = file with one column containing SNP IDs (second column of *.bim file which will be used to create phenotypes),
#								Use full path of "var_ids" file to avoid directory and sub-directory issue
#						nvar = number of variants with causal effects (vector of 3 elements - 
#																		First is # of variants with causal parental effect,
#																		second is # of variants with  causal fetal effect, and
#																		third is # of variants with causal joint effect),
#						var_m = variance of maternal effects
#						var_f = variance of fetal effects
#						covar = covariance of maternal and fetal effects
#				POE = TRUE/FALSE,
#				ppoe = proportion of causal variants with fetal effects showing poes,
#				rho = factor by which effect of m1 or p1 is reduced if eff_source is mf or pf respectively)

create_effects = function(direc, iter, var_ids, nvar, var_m, var_f, covar, POE, ppoe, rho){

	## reproducibility

	#set.seed(12345)   # Seeds random number generator (change integer to change sequence)

	## load packages

	require("MASS");

	## set working directory

	setwd(direc);

	## read vector of causal variants

	ids = read.table(var_ids, head=F, as.is=T);
	names(ids) = "snp_id";
	ids = ids$snp_id;
	nsnps = sum(nvar);
	if(nsnps == 0){
		stop("Select # of causal variants\n");
	}else if(nsnps > length(ids)){
		stop("Number of causal variants to be selected are more than the number of variants provided in the list\n");
	}

	# var_covar = a vector of variance-covariance for parental-fetal effects

	if(is.null(covar)){
		var_covar = NULL;
	}else{
		var_covar = c(var_m, covar, covar, var_f);
	}
	i = 1;
	while(i <= iter){

		## set directoy where effect sizes will be written

		pheno = paste0("pheno", i);
		ll = system(command = paste("find . -type d -name", pheno, sep = " "), intern = TRUE);
		if(is.null(ll)){
			cat("Skipping", pheno, ":", "Working directory must have sub-directory named as", pheno, "\n");
			next;
		}

		## select causal variants

		if(nvar[1]!= 0){
			cv_pt = sample(ids, nvar[1], replace = FALSE);
			lcv_pt = length(cv_pt);
			idx_pt = which(ids %in% cv_pt);
			if(nvar[2] != 0){
				cv_ft = sample(ids[-idx_pt], nvar[2], replace = FALSE);
				lcv_ft = length(cv_ft);
				idx_ft = which(ids %in% cv_ft);
			}else{
				cv_ft = character(0);
				lcv_ft = 0;
				idx_ft = integer(0);
			}
		}else{
			cv_pt = character(0);
			lcv_pt = 0;
			idx_pt = integer(0);
			if(nvar[2] != 0){
				cv_ft = sample(ids, nvar[2], replace = FALSE);
				lcv_ft = length(cv_ft);
				idx_ft = which(ids %in% cv_ft);
			}else{
				cv_ft = character(0);
				lcv_ft = 0;
				idx_ft = integer(0);
			}
		}
			
		if(nvar[3] != 0){
			if(nvar[1] == 0 & nvar[2] == 0){
				cv_joint = sample(ids, nvar[3], replace = FALSE);
			}else{
				cv_joint = sample(ids[-c(idx_pt, idx_ft)], nvar[3], replace = FALSE);
			}
			lcv_joint = length(cv_joint);
		}else{
			cv_joint = character(0);
			lcv_joint = 0;
		}

		## create effects based on eff_source (using standard normal distribution)

		u = matrix(0, nsnps, 2);
		row.names(u) = c(cv_pt, cv_ft, cv_joint);
		colnames(u) = c("par_eff", "fet_eff");

		if(lcv_pt != 0){
			u[cv_pt, 1] = rnorm(lcv_pt, 0, var_m);
		}
		if(lcv_ft != 0){
			u[cv_ft, 2] = rnorm(lcv_ft, 0, var_f);
		}
		if(lcv_joint != 0){
			if(is.null(var_covar)){
				stop("Provide a matrix of variance-covariance for parent-child joint effect\n");
			}
			u_joint = mvrnorm(lcv_joint, rep(0, 2), matrix(var_covar, 2, 2), empirical = TRUE);
			u[cv_joint, 1] = u_joint[, 1];
			u[cv_joint, 2] = u_joint[, 2];
		}
		if(POE){
			u_poe = matrix(u[, 2], nsnps, 1);
			row.names(u_poe) = row.names(u);
			lpoe = round(ppoe * sum(lcv_ft, lcv_joint), 4);
			cv_poe = sample(c(cv_ft, cv_joint), lpoe, replace = FALSE);
			u_poe[cv_poe, 1] = rho * u_poe[cv_poe, 1];
			write.table(u_poe, file.path(".", pheno, "sim_poes.txt"), sep = "\t", quote = FALSE, col.names = FALSE);
		}

		write.table(u[, 1], file.path(".", pheno, "sim_par_effects.txt"), sep = "\t", quote = FALSE, col.names = FALSE);
		write.table(u[, 2], file.path(".", pheno, "sim_fet_effects.txt"), sep = "\t", quote = FALSE, col.names = FALSE);

		cat(pheno, ":", "Causal effects created!\n");

		i = i + 1;
	}
}
