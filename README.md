# H-GCTA
TO start the analyses - simulation or empricial data analyses, first step is to split the haplotypes from imputed vcf files using create_grm.sh which utilizes split_haplotypes and divide_grm binaries.
split_haplotypes and divide_grm are written in c++ and source code can be obtained at https://github.com/PerinatalLab/HAPLOTYPES.
Once haplotypes and corrsponding GRMs are created, further analyses can be performed by using any software like GCTA or LDAK.
Following is the brief instruction for the use of each script. Scripts are written for Linux64-bit operating system on the x86_64 CPU platform.

# create_grm.sh

This is a shell script which takes five arguments -
#arg[1] - working directory where input-output files are present
#arg[2] - prefix of vcf file name
#arg[3] - list of duos/trios
#arg]4] - "trios" or "pat-duos" or "mat-duos"
#arg[5] - maf cutoff (say 0.01)

The script utilizes bcftools, plink, split_haplotypes and divide_grm and expects them to be present in either working directory provided as arg[1] or default path like /usr/bin or usr/local/bin.
Binaries are kept in the working directory to facilitate them running from the default path also, without repeating the codes for calling them from some directory (say ./bin) or without adding that directory to default path.
Pedigree file provided as arg[3] must have 3 columns (without header) with first column as individual id, second column as paternal id and third column as maternal id.
Data type is provided as arg[4]. The script also allows to set a MAF cut-off (arg[5]) for creating plink-format files and creating GRMs.

The first step extracts trios/duos present in imputed vcf file (arg[2]) using pedigree file provided as arg[3] and saves them to another pedigree file named as trios_extract_list.txt in case of trios, mduos_extract_list.txt/pduos_extract_list.txt in case of mat-duos/pat-duos.

In the next step, it extract data for those trios, duos from the imputed vcf file using bcftools.
If bcftools is not present in either working directory or default path, an alternative awk script is used to extract those samples from the vcf file.
This step is computationally expensive, particularly for large data. However it is optional, one can omit this step and directly split haplotypes using raw imputed data.

Once extracted vcf file is created, script asks if the user want to create plink-format haplotype files.
This step only runs if user responds [y/Y] otherwise it only writes the extractd vcf file along with trios, duos pedigree file.

Further, script asks if the user wants to create GRM from extracted vcf file and plink-format haplotype files.
Like above step, this step also runs only if user responds {y/Y] otherwise extracted vcf, pedigree files and plink compatible haplotype files are saved in working directory.

# create_grm_extend.sh

This is another shell script which uses plink-format genotype and haplotype files from individual datasets each with common set of SNPs across all datasets.
The script was originally written for HPC (LSF platform). To run it on personal computer with Linux64 operating system on x86_64 CPU platform, please comment out/remove lines 5-15 in the script.
The script uses an example file of pooled dataset with all polymorphic SNPs (pooled_data_${file}_comm, e.g.: pooled_data_genoM_comm) which is created by merging individual datasets provided in a text file in the working directory - merge_list_${file}.txt (e.g: merge_list_genoM.txt).
**Please replace the input file name (optional for output) in line 28 before running the script!**
**If you choose to rename the output also in line 28, please make sure to find and relace it throughout the script**
**We strongly recommend to retain other output names as they can be directly used in est_h2_with_pcs.sh**

The script creates and modifies GRMs using different models like LDAK-Thin, LDAK (we refer it as LDAK-Weights as it uses SNP specific weights) and Threshold GRM, extracts the set of duos using kinship-coefficient cut-off (0.05). 
For the purpose, script searches for "pooled_mduos_extract_list.txt" in the working directory which is a concatenated file created from "mduos_extract_list.txt" created in individual datasets by the script "create_grm.sh".
It also modifies the GRMs which are used in M-GCTA approach by extracting the GRMs based on mother's genotypes, fetuses' genotypes and maternal-fetal genetic correlation. A R script (mod_mf_grm.R) is called for the purpose which reads the GRM created from maternal-fetal joint data and writes plink/gcta compatible GRMs based on only mothers, fetuses and maternal-fetal genetic correlation.
Eventually, the script also creates matrices of PCs to be used as quantitative covariates in further analyses via GCTA, H-GCTA and M-GCTA analyses.

# est_h2_with_pcs.sh

Like create_grm_extend.sh, this script was also originally created for HPC (LSF platform) and to run it on personal computer with Linux64 operating system on x86_64 CPU platform, please comment out/remove lines 5-10.
The script utilizes the GRMs and PCs created from create_grm_extend.sh and runs REML as well as HE regression analyses using different approach and models for four traits (gday, bw, blen, hc).
H-GCTA and M-GCTA accept list of GRMs without extension saved as hgcta_*_list.txt and mgcta_*_list.txt, respectively.
* represent grm component (e.g.: genoM, genoF, hapM1, hapM2, hapP1, mf), model (e.g.: thin, ldak_weights, ldak_gcta) and alpha value (e.g.: alpha-1.0 or alpha-0.25). Please refer to the script for detail on *.
If user retained the output file names in create_grm_extend.sh, this script can be used just by replacing the trait names (gday, bw, blen, hc).
**please, don't forget to replace the trait names before running the script!**

# simulation

simulation utilizes few R and shell scripts to create simulated data with varrying contribution and correlation of maternal and fetal genetic contribution to phenotypic variance with different levels of parent-of-origin effects (POEs).
POEs were created by simulating matrnal imprinting where maternally transmitted alleles have less effect as compared to paternal transmitted alleles. Various combinations were run; details can found (https://www.biorxiv.org/content/10.1101/2020.05.12.079863v2.full).
This is just an example and user can modify it as per their need.

Simulation utilizes following scripts - 1) create_effects.R, 2) sim_effects_and_pheno_pooled_data.R, 3) sim_effects_and_pheno_pooled_data_run.R, 4) run_reml_pooled_data.sh, 5) process_sim_results.R, and 6) run_simulation.sh
Simulation also uses LDAK5.1.linux or above and gcta64 (version 1.26 or above). Due to large size, we couldn't upload them here. We recommend the user to keep them whrea above scripts are kept.

The example is used for simulating data using pooled dataset and maternal trait but the name of the dataset can be changed and parameters can be modified in run_simulation.sh.
Script automatically captures correlation coefficient of maternal and fetal genetic effects, proportion of causal variants with POEs and imprinting factor (value used to reduce the effect of maternal transmitted alleles as compared to paternal transmitted alleles) from the variable "pheno_dir" (directory name where simulated phenotypes and corresponding results are saved) in run_simulation.sh

1) create_effects.R is a function which creates parental, fetal and parent-of-origin effects for the provided list of causal variants based on the parameters like mean correlation among maternal and fetal genetic effects, relative contribution of maternal and fetal genotypes, proportion of causal variants with POEs and level of maternal imprinting.
2) sim_effects_and_pheno_pooled_data.R is another function which calls create_effects.R to create effects and then simulate phenotypes using empirical data with the help of gcta. This is the function which picks mean correlation among maternal-fetal genetic effects, proportion of causal variants with POEs and level of maternal imprinting from the variable "pheno_dir" (created in run_simulation.sh).
3) sim_effects_and_pheno_pooled_data_run.R calls for previous scripts.
4) run_reml.sh is compilation of commands for running REML using LDAK and GCTA models.
5) process_sim_results.R is used to create table with mean estimates, standard error and corrsponding p values using one sided or two sided z test.
6) run_simulaiton.sh is the script which calls for above scripts and eventually creates tabular results with estimates, SE, and p values.
