# H-GCTA
TO start the analyses - simulation or empricial data analyses, first step is to split the haplotypes from imputed vcf files using create_grm.sh which utilizes split_haplotypes and divide_grm binaries.
split_haplotypes and divide_grm are written in c++ and source code can be obtained at https://github.com/PerinatalLab/HAPLOTYPES.
Once haplotypes and corrsponding GRMs are created, further analyses can be performed.
Following are the brief instruction for the use of each script.

create_grm.sh

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
