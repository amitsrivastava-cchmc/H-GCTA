# H-GCTA
To start the analyses - simulation or empricial data analyses, first step is to split the haplotypes from imputed vcf files using create_grm.sh which utilizes split_haplotypes and divide_grm binaries.
split_haplotypes and divide_grm are written in c++ and source code can be obtained at https://github.com/PerinatalLab/HAPLOTYPES.
Once haplotypes and corrsponding GRMs are created, further analyses can be performed by using any software like GCTA or LDAK. We have provided scripts like create_grm_extend.sh and est_h2_with_pcs.sh (which also calls for mod_mf.R) for that.
Like empirical data analyses, simulation involves several scripts - create_effects.R, sim_effects_and_pheno.R sim_effects_and_pheno_run.R, run_reml.sh, process_sim_results.R, and run_simulation.sh. The last scripts (run_simulation.sh) is the main executable which calls for rest scripts.
# General Recommendations:
Most of the scripts except create_grm.sh and create_effects.R, use "test_extracted" in their input file names, therefore, please use name "test.vcf.gz" for your initial imputed data which automatically generates test_extracted.vcf.gz and future files with same names.
Please keep indexed test.vcf.gz and a list of trio/duos (file must have 3 columsn without header, same as column 2-4 of .bim file) in ~/H-GCTA/data.test.
List of GRMs must have GRM names without extension and with full path (examples of GRM lists are included in data.test, please replace "user" with your actual "user name").
Phenotype file must have 3 columns - FID, IID and Phenotype (without header).
Scripts are written for Linux64-bit operating system on the x86_64 CPU platform.

Following is the brief instruction for the use of each script.  

# emp_data_analysis
  # create_grm.sh

  This is a shell script which takes five arguments -
  #arg[1] - working directory where input-output files are present
  #arg[2] - prefix of vcf file name
  #arg[3] - list of duos/trios
  #arg]4] - "trios" or "pat-duos" or "mat-duos"
  #arg[5] - maf cutoff (say 0.01)

  The script utilizes bcftools, plink, split_haplotypes and divide_grm and expects them to be present in either ~/H-GCTA/bin or default path like /usr/bin or usr/local/bin.
  Pedigree file provided as arg[3] must have 3 columns (without header) with first column as individual id, second column as paternal id and third column as maternal id.
  Data type is provided as arg[4]. The script also allows to set a MAF cut-off (arg[5]) for creating plink-format files and creating GRMs.

  The first step extracts trios/duos present in imputed vcf file (arg[2]) using pedigree file provided as arg[3] and saves them to another pedigree file named as trios_extract_list.txt in case of trios, mduos_extract_list.txt/pduos_extract_list.txt in case of mat-duos/pat-duos.

  In the next step, it extracts data for those trios/duos from the imputed vcf file using bcftools.
  If bcftools is not present in either ~/H-GCTA/bin or default path, an alternative awk script invokes to extract those samples from the vcf file.
  This step is computationally expensive, particularly for large data. However it is optional, one can omit this step and directly split haplotypes using raw imputed data (In case, you don't want to create another vcf file with only trios/duos, please save the original indexed vcf file as "test_extracted.vcf.gz").

  Once extracted vcf file is created, script asks if the user wants to create plink-format haplotype files.
  This step only runs if user responds [y/Y] otherwise the script only writes the extractd vcf file along with trios, duos pedigree file.

  Further, script asks if the user wants to create GRM from extracted vcf file and plink-format haplotype files.
  Like above step, this step also runs only if user responds {y/Y] otherwise extracted vcf, pedigree files and plink compatible haplotype files are saved in working directory. One can say no at this step as this step will be repeated in create_grm_extend.sh.

  # create_grm_extend.sh

  This is another shell script which uses plink-format genotype and haplotype files from mother-child duos.
  The script was originally written for HPC (LSF platform). Lines 5-16 are commented out to run it on personal computer with Linux64 operating system on x86_64 CPU platform which has R and Rscript in computer's $PATH.

  **This script uses "test_extracted" in all input and output names, therefore we strongly recommend to store the original data as indexed test.vcf.gz.**

  The script creates and modifies GRMs using different models like GREML, LDAK-Thin, LDAK (we refer it to as LDAK-Weights as it uses SNP specific weights), extracts the set of duos using kinship-coefficient cut-off (0.05). 
  For the latter one, script searches for "mduos_extract_list.txt" which was created via "create_grm.sh".
  It also modifies the GRMs which are used in M-GCTA approach by extracting the GRMs based on mother's genotypes, fetuses' genotypes and maternal-fetal genetic correlation. An R script (mod_mf_grm.R) is called for the purpose which reads the GRM created from maternal-  fetal joint data and writes plink/gcta compatible GRMs based on only mothers, fetuses and maternal-fetal genetic correlation.
  Eventually, the script also creates matrices of PCs to be used as quantitative covariates in further analyses via GCTA, H-GCTA and M-GCTA approach.

  # est_h2_with_pcs.sh

  Like create_grm_extend.sh, this script was also originally created for HPC (LSF platform) and lines 5-16 are commented out to run it on personal computer with Linux64 operating system on x86_64 CPU platform.
  The script utilizes the GRMs and PCs created from create_grm_extend.sh and runs REML as well as HE regression analyses using different approach and models for input phenotype (provided via user input).
  H-GCTA and M-GCTA accept list of GRMs without extension saved as hgcta_[model]_mduos_[alpha-1.0/alpha-0.25]_grm_list.txt and mgcta_[model]_mduos_[alpha-1.0/alpha-0.25]_grm_list.txt, respectively.
  Model (e.g.: ladka_thin, ldak_weights, ldak_gcta) and alpha value (e.g.: alpha-1.0 or alpha-0.25) depend on the analysis. Please refer to the example files provided in ~/H-GCTA/data.test.

  **This script also uses "test_extracted" in all input and output names, therefore we strongly recommend to store the original data as indexed test.vcf.gz.**

# simulation

simulation utilizes few R and shell scripts to create simulated data with varrying contribution and correlation of maternal and fetal genetic contribution to phenotypic variance along with different levels of parent-of-origin effects (POEs).
POEs were created by simulating matrnal imprinting where maternally transmitted alleles have less effect as compared to paternal transmitted alleles. Analyses were performed using various simulation conditions; details can be found (https://www.biorxiv.org/content/10.1101/2020.05.12.079863v2.full).
Simulation utilizes following scripts - 1) create_effects.R, 2) sim_effects_and_pheno.R, 3) sim_effects_and_pheno_run.R, 4) run_reml.sh, 5) process_sim_results.R, and 6) run_simulation.sh

**Most of the simulaiton scripts (except create_effects.R) use "test_extracted" in their input names, therefore please save your original idexed vcf file as "test.vcf.gz".**
**User needs to execute just one script run_simulation.sh which calls for rest scripts.**

Analyses performed via run_simulation.sh use LDAK6.linux and gcta64 (version 1.94.1-linux-kernel-3-x86_64).
Like previous scripts, lines 5-12 are commented out to run the script on any personal computer which has R and Rscript in computer's "$PATH".
Lines 28-33 are the most crucial variables which decide what would be the number of iterations and relative maternal and fetal genetic contribution to dyadic traits.
When prompted to provide phenotype directory, please provide the directory name where you want the simulated phenotypes and corresponding results to be saved (e.g. /home/user/H-GCTA/simulation/simulation_mf_with_ppoe_1.0_m1_0.0_mf_cor_1.0)
Phenotype directory name has a unique style to represent the type of trait like mat/fet/mf, ppoe (0-1.0), effect of m1 as compared to p1, and correlation (mf_cor) for example, [(-1.0) - 1.0].
Script automatically captures correlation coefficient of maternal and fetal genetic effects, proportion of causal variants with POEs (ppoe) and imprinting factor (value used to reduce the effect of maternal transmitted alleles as compared to paternal transmitted alleles) from the variable "pheno_dir" (directory name where simulated phenotypes and corresponding results are saved) in run_simulation.sh.
Script also asks the user to provide the directory where plink format files and corresponding GRMs are saved (/home/user/H-GCTA/data.test).

  # Brief description of scripts used in simulation

  1) create_effects.R is a function which creates parental, fetal and parent-of-origin effects using the provided list of causal variants based on the parameters like mean correlation among maternal and fetal genetic effects (mf_cor), relative contribution of maternal and fetal genotypes, proportion of causal variants with POEs (ppoe) and level of maternal imprinting (effect of m1 as compared to p1).
  2) sim_effects_and_pheno.R is another function which calls create_effects.R to create effects and then simulate phenotypes using empirical data with the help of gcta. This is the function which picks mean correlation among maternal-fetal genetic effects (mf_cor), proportion of causal variants with POEs (ppoe) and level of maternal imprinting from the variable "pheno_dir" (created in run_simulation.sh).
  3) sim_effects_and_pheno_run.R calls for the above script ("sim_effects_and_pheno.R").
  4) run_reml.sh is compilation of commands for running REML using LDAK and GCTA models. This script expects that GRMs should be saved in directory provided as arg[1] and phenotypes in a directory provided as arg[2]. arg[3] is the number of iterations.
  5) process_sim_results.R is used to create table with mean estimates, standard error (SE) and corrsponding p values using one sided and two sided z test.
  6) run_simulaiton.sh is the main script which calls for above scripts and eventually creates tabular results with estimates, SE, and p values.
