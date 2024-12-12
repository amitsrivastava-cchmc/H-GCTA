#following commands take output from GCTA or LDAK and create summary results based on mean heritability and error from all simulations.

args = commandArgs(trailingOnly = TRUE);
dir_name = args[1];

z.test = function(mean1, mean2, se1, se2){
		mean= abs(mean1 - mean2);
		se = sqrt((se1)^2 + (se2)^2);
		z = mean/se;
		p.onesided = 1-pnorm(z);
		p.twosided = 2 * p.onesided;
		cbind(p.onesided, p.twosided);
	}


#process results from LDAK-Thin (alpha = -1.0)

m = read.table(file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_genoM.txt"), head=F, as.is=T);
f = read.table(file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_genoF.txt"), head=F, as.is=T);
m_only = read.table(file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_M_var.txt"), head=F, as.is=T);
f_only = read.table(file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_F_var.txt"), head=F, as.is=T);
mf = read.table(file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_MF_covar.txt"), head=F, as.is=T);
m1 = read.table(file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_M1.txt"), head=F, as.is=T);
m2 = read.table(file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_M2.txt"), head=F, as.is=T);
p1 = read.table(file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_P1.txt"), head=F, as.is=T);

dd = data.frame("Var" = rep(0, 8), "S.E." = rep(0, 8), "P-Val (One-tailed)" = rep(0, 8), "P-Val (Two-tailed)" = rep(0, 8));

dd[1,c(1,2)] = c(mean(m[,2]), sd(m[, 2]));
dd[2,c(1,2)] = c(mean(f[,2]), sd(m[, 2]));
dd[3,c(1,2)] = c(mean(m_only[,2]), sd(m_only[, 2]));
dd[4,c(1,2)] = c(mean(f_only[,2]), sd(f_only[, 2]));
dd[5,c(1,2)] = c(mean(mf[,2]), sd(mf[, 2]));
dd[6,c(1,2)] = c(mean(m1[,2]), sd(m1[, 2]));
dd[7,c(1,2)] = c(mean(m2[,2]), sd(m2[, 2]));
dd[8,c(1,2)] = c(mean(p1[,2]), sd(p1[, 2]));

aa = sapply(1:nrow(dd), function(x) z.test(dd[x, 1], 0, dd[x, 2], 0));
dd[, c(3,4)] = t(aa);

write.table(dd, file.path(dir_name, "ldak_thin_alpha-1.0_sim_results_all.txt"), sep = "\t", quote = FALSE);

#process results from LDAK-Thin (alpha = -0.25)

m = read.table(file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_genoM.txt"), head=F, as.is=T);
f = read.table(file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_genoF.txt"), head=F, as.is=T);
m_only = read.table(file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_M_var.txt"), head=F, as.is=T);
f_only = read.table(file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_F_var.txt"), head=F, as.is=T);
mf = read.table(file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_MF_covar.txt"), head=F, as.is=T);
m1 = read.table(file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_M1.txt"), head=F, as.is=T);
m2 = read.table(file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_M2.txt"), head=F, as.is=T);
p1 = read.table(file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_P1.txt"), head=F, as.is=T);

dd = data.frame("Var" = rep(0, 8), "S.E." = rep(0, 8), "P-Val (One-tailed)" = rep(0, 8), "P-Val (Two-tailed)" = rep(0, 8));

dd[1,c(1,2)] = c(mean(m[,2]), sd(m[, 2]));
dd[2,c(1,2)] = c(mean(f[,2]), sd(m[, 2]));
dd[3,c(1,2)] = c(mean(m_only[,2]), sd(m_only[, 2]));
dd[4,c(1,2)] = c(mean(f_only[,2]), sd(f_only[, 2]));
dd[5,c(1,2)] = c(mean(mf[,2]), sd(mf[, 2]));
dd[6,c(1,2)] = c(mean(m1[,2]), sd(m1[, 2]));
dd[7,c(1,2)] = c(mean(m2[,2]), sd(m2[, 2]));
dd[8,c(1,2)] = c(mean(p1[,2]), sd(p1[, 2]));

aa = sapply(1:nrow(dd), function(x) z.test(dd[x, 1], 0, dd[x, 2], 0));
dd[, c(3,4)] = t(aa);

write.table(dd, file.path(dir_name, "ldak_thin_alpha-0.25_sim_results_all.txt"), sep = "\t", quote = FALSE);

#process results from LDAK-Weights (alpha = -1.0)

m = read.table(file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_genoM.txt"), head=F, as.is=T);
f = read.table(file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_genoF.txt"), head=F, as.is=T);
m_only = read.table(file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_M_var.txt"), head=F, as.is=T);
f_only = read.table(file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_F_var.txt"), head=F, as.is=T);
mf = read.table(file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_MF_covar.txt"), head=F, as.is=T);
m1 = read.table(file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_M1.txt"), head=F, as.is=T);
m2 = read.table(file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_M2.txt"), head=F, as.is=T);
p1 = read.table(file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_P1.txt"), head=F, as.is=T);

dd = data.frame("Var" = rep(0, 8), "S.E." = rep(0, 8), "P-Val (One-tailed)" = rep(0, 8), "P-Val (Two-tailed)" = rep(0, 8));

dd[1,c(1,2)] = c(mean(m[,2]), sd(m[, 2]));
dd[2,c(1,2)] = c(mean(f[,2]), sd(m[, 2]));
dd[3,c(1,2)] = c(mean(m_only[,2]), sd(m_only[, 2]));
dd[4,c(1,2)] = c(mean(f_only[,2]), sd(f_only[, 2]));
dd[5,c(1,2)] = c(mean(mf[,2]), sd(mf[, 2]));
dd[6,c(1,2)] = c(mean(m1[,2]), sd(m1[, 2]));
dd[7,c(1,2)] = c(mean(m2[,2]), sd(m2[, 2]));
dd[8,c(1,2)] = c(mean(p1[,2]), sd(p1[, 2]));

aa = sapply(1:nrow(dd), function(x) z.test(dd[x, 1], 0, dd[x, 2], 0));
dd[, c(3,4)] = t(aa);

write.table(dd, file.path(dir_name, "ldak_weights_alpha-1.0_sim_results_all.txt"), sep = "\t", quote = FALSE);

#process results from LDAK-Weights (alpha = -0.25)

m = read.table(file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_genoM.txt"), head=F, as.is=T);
f = read.table(file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_genoF.txt"), head=F, as.is=T);
m_only = read.table(file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_M_var.txt"), head=F, as.is=T);
f_only = read.table(file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_F_var.txt"), head=F, as.is=T);
mf = read.table(file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_MF_covar.txt"), head=F, as.is=T);
m1 = read.table(file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_M1.txt"), head=F, as.is=T);
m2 = read.table(file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_M2.txt"), head=F, as.is=T);
p1 = read.table(file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_P1.txt"), head=F, as.is=T);

dd = data.frame("Var" = rep(0, 8), "S.E." = rep(0, 8), "P-Val (One-tailed)" = rep(0, 8), "P-Val (Two-tailed)" = rep(0, 8));

dd[1,c(1,2)] = c(mean(m[,2]), sd(m[, 2]));
dd[2,c(1,2)] = c(mean(f[,2]), sd(m[, 2]));
dd[3,c(1,2)] = c(mean(m_only[,2]), sd(m_only[, 2]));
dd[4,c(1,2)] = c(mean(f_only[,2]), sd(f_only[, 2]));
dd[5,c(1,2)] = c(mean(mf[,2]), sd(mf[, 2]));
dd[6,c(1,2)] = c(mean(m1[,2]), sd(m1[, 2]));
dd[7,c(1,2)] = c(mean(m2[,2]), sd(m2[, 2]));
dd[8,c(1,2)] = c(mean(p1[,2]), sd(p1[, 2]));

aa = sapply(1:nrow(dd), function(x) z.test(dd[x, 1], 0, dd[x, 2], 0));
dd[, c(3,4)] = t(aa);

write.table(dd, file.path(dir_name, "ldak_weights_alpha-0.25_sim_results_all.txt"), sep = "\t", quote = FALSE);

#process results from LDAK-GCTA (alpha = -1.0)

m = read.table(file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_genoM.txt"), head=F, as.is=T);
f = read.table(file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_genoF.txt"), head=F, as.is=T);
m_only = read.table(file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_M_var.txt"), head=F, as.is=T);
f_only = read.table(file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_F_var.txt"), head=F, as.is=T);
mf = read.table(file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_MF_covar.txt"), head=F, as.is=T);
m1 = read.table(file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_M1.txt"), head=F, as.is=T);
m2 = read.table(file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_M2.txt"), head=F, as.is=T);
p1 = read.table(file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_P1.txt"), head=F, as.is=T);

dd = data.frame("Var" = rep(0, 8), "S.E." = rep(0, 8), "P-Val (One-tailed)" = rep(0, 8), "P-Val (Two-tailed)" = rep(0, 8));

dd[1,c(1,2)] = c(mean(m[,2]), sd(m[, 2]));
dd[2,c(1,2)] = c(mean(f[,2]), sd(m[, 2]));
dd[3,c(1,2)] = c(mean(m_only[,2]), sd(m_only[, 2]));
dd[4,c(1,2)] = c(mean(f_only[,2]), sd(f_only[, 2]));
dd[5,c(1,2)] = c(mean(mf[,2]), sd(mf[, 2]));
dd[6,c(1,2)] = c(mean(m1[,2]), sd(m1[, 2]));
dd[7,c(1,2)] = c(mean(m2[,2]), sd(m2[, 2]));
dd[8,c(1,2)] = c(mean(p1[,2]), sd(p1[, 2]));

aa = sapply(1:nrow(dd), function(x) z.test(dd[x, 1], 0, dd[x, 2], 0));
dd[, c(3,4)] = t(aa);

write.table(dd, file.path(dir_name, "ldak_gcta_alpha-1.0_sim_results_all.txt"), sep = "\t", quote = FALSE);

#process results from LDAK-GCTA (alpha = -0.25)

m = read.table(file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_genoM.txt"), head=F, as.is=T);
f = read.table(file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_genoF.txt"), head=F, as.is=T);
m_only = read.table(file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_M_var.txt"), head=F, as.is=T);
f_only = read.table(file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_F_var.txt"), head=F, as.is=T);
mf = read.table(file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_MF_covar.txt"), head=F, as.is=T);
m1 = read.table(file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_M1.txt"), head=F, as.is=T);
m2 = read.table(file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_M2.txt"), head=F, as.is=T);
p1 = read.table(file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_P1.txt"), head=F, as.is=T);

dd = data.frame("Var" = rep(0, 8), "S.E." = rep(0, 8), "P-Val (One-tailed)" = rep(0, 8), "P-Val (Two-tailed)" = rep(0, 8));

dd[1,c(1,2)] = c(mean(m[,2]), sd(m[, 2]));
dd[2,c(1,2)] = c(mean(f[,2]), sd(m[, 2]));
dd[3,c(1,2)] = c(mean(m_only[,2]), sd(m_only[, 2]));
dd[4,c(1,2)] = c(mean(f_only[,2]), sd(f_only[, 2]));
dd[5,c(1,2)] = c(mean(mf[,2]), sd(mf[, 2]));
dd[6,c(1,2)] = c(mean(m1[,2]), sd(m1[, 2]));
dd[7,c(1,2)] = c(mean(m2[,2]), sd(m2[, 2]));
dd[8,c(1,2)] = c(mean(p1[,2]), sd(p1[, 2]));

aa = sapply(1:nrow(dd), function(x) z.test(dd[x, 1], 0, dd[x, 2], 0));
dd[, c(3,4)] = t(aa);

write.table(dd, file.path(dir_name, "ldak_gcta_alpha-0.25_sim_results_all.txt"), sep = "\t", quote = FALSE);

#process results from GCTA (alpha = -1.0)

m = read.table(file.path(dir_name, "gcta_sim_results_genoM.txt"), head=F, as.is=T);
f = read.table(file.path(dir_name, "gcta_sim_results_genoF.txt"), head=F, as.is=T);
m_only = read.table(file.path(dir_name, "gcta_sim_results_M_var.txt"), head=F, as.is=T);
f_only = read.table(file.path(dir_name, "gcta_sim_results_F_var.txt"), head=F, as.is=T);
mf = read.table(file.path(dir_name, "gcta_sim_results_MF_covar.txt"), head=F, as.is=T);
m1 = read.table(file.path(dir_name, "gcta_sim_results_M1.txt"), head=F, as.is=T);
m2 = read.table(file.path(dir_name, "gcta_sim_results_M2.txt"), head=F, as.is=T);
p1 = read.table(file.path(dir_name, "gcta_sim_results_P1.txt"), head=F, as.is=T);

dd = data.frame("Var" = rep(0, 8), "S.E." = rep(0, 8), "P-Val (One-tailed)" = rep(0, 8), "P-Val (Two-tailed)" = rep(0, 8));

dd[1,c(1,2)] = c(mean(m[,2]), sd(m[, 2]));
dd[2,c(1,2)] = c(mean(f[,2]), sd(m[, 2]));
dd[3,c(1,2)] = c(mean(m_only[,2]), sd(m_only[, 2]));
dd[4,c(1,2)] = c(mean(f_only[,2]), sd(f_only[, 2]));
dd[5,c(1,2)] = c(mean(mf[,2]), sd(mf[, 2]));
dd[6,c(1,2)] = c(mean(m1[,2]), sd(m1[, 2]));
dd[7,c(1,2)] = c(mean(m2[,2]), sd(m2[, 2]));
dd[8,c(1,2)] = c(mean(p1[,2]), sd(p1[, 2]));

aa = sapply(1:nrow(dd), function(x) z.test(dd[x, 1], 0, dd[x, 2], 0));
dd[, c(3,4)] = t(aa);

write.table(dd, file.path(dir_name, "gcta_sim_results_all.txt"), sep = "\t", quote = FALSE);
