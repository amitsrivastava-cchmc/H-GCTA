# H-GCTA
TO start the analyses - simulation or empricial data analyses, first step is to split the haplotypes from imputed vcf files using create_grm.sh which utilizes split_haplotypes and divide_grm binaries.
split_haplotypes and divide_grm are written in c++ and source code can be obtained at https://github.com/PerinatalLab/HAPLOTYPES.
Once haplotypes and corrsponding GRMs are created, further analyses can be performed by using any software like GCTA or LDAK.
Following is the brief instruction for the use of each script. Scripts are written for Linux64-bit operating system on the x86_64 CPU platform and use binaries like Plink1.9, gcta1.26 or above.

# create_grm.sh

This is a shell script which takes five arguments -
#arg[1] - working directory where input-output files are present
#arg[2] - prefix of vcf file name
#arg[3] - list of individuals/duos/trios
#arg]4] - "trios" or "pat-duos" or "mat-duos"
#arg[5] - maf cutoff (say 0.01)

The script utilizes bcftools, plink, split_haplotypes and divide_grm and expects them to be present in either working directory provided as arg[1] or default path like /usr/bin or usr/local/bin.
Binaries are kept in the working directory to facilitate them running from the default path also, without repeating the codes for calling them from some directory (say ./bin) or without adding that directory to default path.
Pedigree file provided as arg[3] must have 3 columns (without header) with first column as individual id, second column as paternal id and third column as maternal id.
Data type is provided as arg[4]. The script also allows to set a MAF cut-off (arg[5]) for creating plink-format files and creating GRMs.

The first step extracts trios/duos present in imputed vcf file (arg[2]) using pedigree file provided as arg[3] and saves them another pedigree file named as trios_extract_list.txt in case of trios, mduos_extract_list.txt/pduos_extract_list.txt in case of mat-duos/pat-duos.

In the next step, it extract data for those trios, duos from the imputed vcf file using bcftools.
If bcftools is not present in either working directory or default path, an alternative awk script is used to extract those samples from the vcf file.
This step is computationally expensive, particularly for large data. However it is optional, one can omit this step and directly split haplotypes using raw imputed data.

Once extracted vcf file is created, script asks if the user want to create plink-format haplotype files.
This step only runs if user responds [y/Y] otherwise it only writes the extractd vcf file along with trios, duos pedigree file.

Further, script asks if the user wants to create GRM from extracted vcf file and plink-format haplotype files.
Like above step, this step also runs only if user responds {y/Y] otherwise extracted vcf, pedigree files and plink compatible haplotype files are saved in working directory.

# create_grm_extend.sh

This is another shell script which uses plink-format genotype and haplotype files from individual dataset each with common set of SNPs across all datasets.
The script was originally written for HPC (LSF platform). To run it on personal computer with Linux64 operating system on x86_64 CPU platform, please comment out/remove lines 5-15 from the script.
The script uses an example file of pooled dataset with all polymorphic SNPs which is created by merging individual datasets provided in a text file merge_list_{file}.txt (e.g: merge_list_genoM.txt).
It does several things such as creating and modifying GRMs using different models like LDAK-Thin, LDAK (we refer it as LDAK-Weights as it uses SNP specific weights) and Threshold GRM.
It extracts the set of duos using kinship-coefficient cut-off (0.05). For the purpose, script searches for "pooled_mduos_extract_list.txt" in the working directory which is a concatenated file created from "mduos_extract_list.txt" created in individual datasets by the script "create_grm.sh".
It also modifies the GRMs which are used in M-GCTA approach by extracting the GRMs based on mother's genotypes, fetuses' genotypes and maternal-fetal genetic correlation. A R script (mod_mf_grm.R) is called for the purpose which reads the GRM created from maternal-fetal joint data and writes plink/gcta compatible GRMs based on only mothers, fetuses and maternal-fetal genetic correlation.
Eventually, the script also creates matrices of PCs to be used as quantitative covariates in GCTA, H-GCTA and M-GCTA analyses.

# est_h2_with_pcs.sh

Like previous script (create_grm_extend.sh), this script is also specific to our manuscript (https://www.biorxiv.org/content/10.1101/2020.05.12.079863v2.full), which utilizes the GRMs and PCs created from create_grm_extend.sh and runs REML analyses using different approach and models for four traits (gday, bw, blen, hc).
This script was also originally created for HPC (LSF platform) and to run it on personal computer with Linux64 operating system on x86_64 CPU platform, please comment out/remove lines 5-10.

