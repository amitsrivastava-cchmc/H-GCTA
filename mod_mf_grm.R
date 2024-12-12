#!/usr/bin/env Rscript

#following are commands to extract GRMs bsed on maternal genotypes, fetal genotypes and maternal-fetal genetic correlation from the mother-child joint GRM.
#function "ReadGRMBin" is adoopted from the GCTA.

args = commandArgs(trailingOnly = T);
dir = args[1];
input = args[2];
input = paste(dir, input, sep = "/");
output = sub("_grm.*", "", input);
m.bin = paste(output, "mother_grm.grm.bin", sep =  "_");
f.bin = paste(output, "fetus_grm.grm.bin", sep = "_");
mf.bin = paste(output, "mf_grm.grm.bin", sep = "_");
m.N.bin = paste(output, "mother_grm.grm.N.bin", sep =  "_");
f.N.bin = paste(output, "fetus_grm.grm.N.bin", sep = "_");
mf.N.bin = paste(output, "mf_grm.grm.N.bin", sep = "_");
m.id = paste(output, "mother_grm.grm.id", sep =  "_");
f.id = paste(output, "fetus_grm.grm.id", sep = "_");
mf.id = paste(output, "mf_grm.grm.id", sep = "_");

setwd(dir);

ReadGRMBin=function(prefix, AllN=F, size=4){
		sum_i=function(i){
			return(sum(1:i));
		}
		BinFileName=paste(prefix,".grm.bin",sep="");
		NFileName=paste(prefix,".grm.N.bin",sep="");
		IDFileName=paste(prefix,".grm.id",sep="");
		id = read.table(IDFileName);
		n=dim(id)[1];
		BinFile=file(BinFileName, "rb");
		grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size);
                i=sapply(1:n, sum_i);
                close(BinFile);
		if(file.exists(NFileName)){
			NFile=file(NFileName, "rb");
			if(AllN==T){
				M=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size);
			}else{
				M=readBin(NFile, n=1, what=numeric(0), size=size);
			}
			close(NFile);
			return(list(diag=grm[i], off=grm[-i], id=id, M=M));
		}else{
			return(list(diag=grm[i], off=grm[-i], id=id));
		}
}

grm = ReadGRMBin(input);
ind = length(grm$id$V2);
#mum = grm$id[1:(ind/2), ];
child = grm$id[(ind/2 + 1):ind, ];
#head(mum, 10); head(child, 10);
mat = matrix(0, ind, ind);
mat[upper.tri(mat, diag=F)] = grm$off;
mat = t(mat) + mat;
diag(mat) = grm$diag;
#mat[1:5,1:5];
m = mat[1:(ind/2), 1:(ind/2)];
f = mat[(ind/2 + 1):ind, (ind/2 + 1):ind];
mf = mat[(ind/2 + 1):ind, 1:(ind/2)] + t(mat[(ind/2 + 1):ind, 1:(ind/2)]);
#m[1:5,1:5]; f[1:5,1:5]; mf[1:5,1:5];
m_1 = m[upper.tri(m, diag=T)];
f_1 = f[upper.tri(f, diag=T)];
mf_1 = mf[upper.tri(mf, diag=T)];
writeBin(m_1, m.bin, size = 4);
writeBin(f_1, f.bin, size = 4);
writeBin(mf_1, mf.bin, size = 4);
write.table(child, m.id, quote=F, sep="\t", row.names=F, col.names=F);
write.table(child, f.id, quote=F, sep="\t", row.names=F, col.names=F);
write.table(child, mf.id, quote=F, sep="\t", row.names=F, col.names=F);

if(!is.null(grm$M)){
	writeBin(grm$M, m.N.bin, size = 4);
	writeBin(grm$M, f.N.bin, size = 4);
	writeBin(grm$M, mf.N.bin, size = 4);
}
