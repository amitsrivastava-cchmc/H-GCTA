#!/bin/bash

#if running on cluster first load bcftools, gcc and glibc as below
#module load bcftools/1.8
#module load gcc/5.4.0
#module load glibc/2.14

#arg[1] - working directory where input-output files are present
#arg[2] - prefix of vcf file name
#arg[3] - list of individuals/duos/trios
#arg]4] - "trios" or "pat-duos" or "mat-duos" "individuals"
#arg[5] - maf cutoff (say 0.000001)

#read working directory from input

cd ${1}
echo -e "Outputs will be generated in $(pwd)\n"

#list of samples present in imputed vcf file

gunzip -c ${2}.vcf.gz | head -n 200 | awk '$1 ~ /^#CHROM/{
								for(i=10; i<=NF; i++){
									print $i
								}
							}' | sort -h -k 1,1 > ${2}_ids.txt

#extract samples common between input samples list (arg[3]) and vcf_ids.txt
#sample list must have at least 3 columns: col I - iid, col II - pid, col III - mid

col=$(echo $(awk 'NR==1{print NF}' ${3}))

if [[ $col -lt 3 ]]
then
	echo -e "\nsample list must have at least 3 columns, check sample file\n"
	exit 0
fi

if [ ${4} = "trios" ]
then
	awk 'BEGIN{OFS="\t"}
		NR==FNR{
			x[$1] = $1
			next
		}
		{
			if($1 in x && $2 in x && $3 in x){
				print $0
			}
		}' ${2}_ids.txt ${3} | sort -h -k 1,1 > trios_extract_list.txt
	awk '{print $1 "\n" $2 "\n" $3}' trios_extract_list.txt | sort -h -k 1,1 > ind_to_extract.txt
elif [	${4} = "pat-duos" ]
then
	awk 'BEGIN{OFS="\t"}
		NR==FNR{
			x[$1] = $1
			next
		}
		{
			if($1 in x && $2 in x){
				print $0
			}
		}' ${2}_ids.txt ${3} | sort -h -k 1,1 > pduos_extract_list.txt
	awk '{print $1 "\n" $2}' pduos_extract_list.txt | sort -h -k 1,1 > ind_to_extract.txt
elif [  ${4} = "mat-duos" ]
then
	awk 'BEGIN{OFS="\t"}
		NR==FNR{
			x[$1] = $1
			next
		}
		{
			if($1 in x && $3 in x){
				print $0
			}
		}' ${2}_ids.txt ${3} | sort -h -k 1,1 > mduos_extract_list.txt
	awk '{print $1 "\n" $3}' mduos_extract_list.txt | sort -h -k 1,1 > ind_to_extract.txt
#else
#	awk 'BEGIN{OFS="\t"}
#		NR==FNR{
#			x[$1] = $1
#			next
#		}
#		{
#			if($1 in x && !($2 in x || $3 in x)){
#				print $1 OFS 0 OFS 0
#			}
#		}' ${2}_ids.txt ${3} | sort -h -k 1,1 > unrelatd_extract_list.txt
#	awk '{print $1}' | sort -h -k 1,1 > ind_to_extract.txt
fi

#extract list of samples from vcf file using either bcftools or awk script. 
#awk script is 2-3 times faster than bcftools.
#however, bcftools can also extract [variants maf[0] > 0.000001] simultaneously.

if [ "$(command -v ./bcftools)" ]
then
	echo -e "bcftools is executed from ${1}\n"
	../bin/bcftools view -S ind_to_extract.txt ${2}.vcf.gz --threads 32 -Oz > ${2}_extracted.vcf.gz
elif [ "$(command -v bcftools)" ]
then
	echo -e "bcftools is executed from " $(command -v bcftools)"\n"
	bcftools view -S ind_to_extract.txt ${2}.vcf.gz --threads 32 -Oz > ${2}_extracted.vcf.gz
else
	echo -e "bcftools isn't present in ${1} or any default path\n"
	echo -e "awk script will be used\n"
	awk 'BEGIN{OFS = "\t"}
		NR==FNR{
			x[$1]++;
			next
		}
		$1 ~ /^##/{
			print;
			next
		}
		$1 ~ /#CHROM/{
			for(i=10; i<=NF; i++){
				if($i in x){
					y[i]++
				}
			}
		}
		{
			for(i=1; i<=9; i++){
				printf "%s\t", $i
			}
			for (j in y){
				printf "%s\t", $j
			}
			printf "\n"
		}'  ind_to_extract.txt <(zcat ${2}.vcf.gz) | gzip > ${2}_extracted.vcf.gz
fi

#split extracted vcf file into plink compatible haplotype files.

read -p "Do you want to split extracted vcf file into plink compatible haplotype files [y/n]: "
if [[ $REPLY == "Y" || $REPLY == "y" ]]
then
	if [ "$(command -v ./split_haplotypes)" ] || [ "$(command -v split_haplotypes)" ]
	then
		if [ ${4} = "trios" ]
		then
			echo -e "\nConverting VCFs (trios) to haplotype BEDs..."
			zcat ${2}_extracted.vcf.gz | split_haplotypes trios_extract_list.txt ${2}_extracted_hap.log ${2}_extracted_hap
			echo -e "Haplotype BEDs produced\n"
		elif [ ${4} = "pat-duos" ]
		then
			echo -e "\nConverting VCFs (pat-duos) to haplotype BEDs..."
			zcat ${2}_extracted.vcf.gz | split_haplotypes pduos_extract_list.txt ${2}_extracted_hap.log ${2}_extracted_hap
			echo -e "Haplotype BEDs produced\n"
		elif [ ${4} = "mat-duos" ]
		then
			echo -e "\nConverting VCFs (mat-duos) to haplotype BEDs..."
			zcat ${2}_extracted.vcf.gz | split_haplotypes mduos_extract_list.txt ${2}_extracted_hap.log ${2}_extracted_hap
			echo -e "Haplotype BEDs produced\n"
		fi
	else
		echo -e "'split_haplotype' isn't present in ${1} or any default path\n"
		echo -e "extracted vcf file is saved in ${1}\n"
		exit 0
	fi
else
	echo -e "extracted vcf file is saved in ${1}\n"
	exit 0
fi

#modify bim files (use it only for genotype bim files not for haplotypes)

function modbim {
	#Enter bim file name (without extension)

	awk 'BEGIN{OFS="\t"}
		{
			if($2 == "."){
				$2 = $1 ":" $4 ":" $5 ":" $6
			}else{
				$2 = $1 ":" $2 ":" $4 ":" $5 ":" $6
			}
			print
		}' ${1}.bim > ${1}_tmp.txt
	mv ${1}_tmp.txt ${1}.bim
}
#create and adjust GRM

#Since 0/2 encoding of haplotype files overestimates the kinship coefficient by a factor of 2, therefore GRM created from splitted haplotypes need to be divided by a factor of 2

#function for adjusting GRMs by 2
#arg[1] - input GRM
#arg[2] - number of samples in fam file
#arg[3] - division constant

function adjustGRM {
	echo -e "\nAdjusting GRM ${1}.grm.bin by a factor of ${3}...\n"
	if [ "$(command -v ./divide_grm)" ]
	then
		./divide_grm ${1}.grm.bin ${1}.grm.bin_temp ${2} ${3} 
		mv ${1}.grm.bin_temp ${1}.grm.bin
	elif [ "$(command -v divide_grm)" ]
	then
		divide_grm ${1}.grm.bin ${1}.grm.bin_temp ${2} ${3} 
		mv ${1}.grm.bin_temp ${1}.grm.bin
	else
		"'divide_grm' isn't present in "$(pwd)" or any default path\n"
		exit 0
	fi
}

read -p "Do you want to create GRM from extracted vcf file and plink compatible haplotype files [y/n]: "
if [[ $REPLY == "Y" || $REPLY == "y" ]]
then
	if [ "$(command -v ./plink)" ] || [ "$(command -v plink)" ]
	then
        if [ ${4} = "trios" ]
        then
			echo -e "\nCalculating GRMs from BEDs...\n"

			#produce standard GRM for genotypes
			#(use --keep-allele-order if want to keep ref/alt assignments)

			echo -e "\nProducing maternal-genotype BED above set MAF cutoff...\n"

			plink --vcf ${2}_extracted.vcf.gz --keep ${2}_extracted_hapM.tokeep \
				--maf ${5} --keep-allele-order --indiv-sort f ${2}_extracted_hapM.tokeep --make-bed --out ${2}_extracted_genoM

			awk '{print $1 ":" $4 ":" $6 ":" $5}' ${2}_extracted_genoM.bim > ${2}_snps_gt_${5}.txt

			modbim ${2}_extracted_genoM
			awk '{print $2}' ${2}_extracted_genoM > ${2}_extracted_geno_snps_list.txt

			echo -e "\nProducing maternal-genotype GRM...\n"

			plink --bfile ${2}_extracted_genoM --make-grm-bin --out ${2}_extracted_genoM_grm

			echo -e "\nProducing paternal-genotype BED above set MAF cutoff...\n"

			plink --vcf ${2}_extracted.vcf.gz --keep ${2}_extracted_hapP.tokeep \
				--keep-allele-order --make-bed --out ${2}_extracted_genoP_tmp

			modbim ${2}_extracted_genoP_tmp

			echo -e "\n"

			plink --bfile ${2}_extracted_genoP_tmp --extract ${2}_extracted_geno_snps_list.txt \
				--keep-allele-order --indiv-sort f ${2}_extracted_hapP.tokeep --make-bed --out ${2}_extracted_genoP

			echo -e "\nProducing paternal-genotype GRM...\n"

			plink --bfile ${2}_extracted_genoP --make-grm-bin --out ${2}_extracted_genoP_grm

			rm ${2}_extracted_genoP_tmp.*

			echo -e "\nProducing fetal-genotype BED above set MAF cutoff...\n"

			plink --vcf ${2}_extracted.vcf.gz --keep ${2}_extracted_hapF.fam \
				--maf ${5} --keep-allele-order --make-bed --out ${2}_extracted_genoF

			modbim ${2}_extracted_genoF_tmp

			echo -e "\n"

			plink --bfile ${2}_extracted_genoF_tmp --extract ${2}_extracted_geno_snps_list.txt \
				--keep-allele-order --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_genoF

			rm ${2}_extracted_genoF_tmp.*

			echo -e "\nProducing fetal-genotype GRM...\n"

			plink --bfile ${2}_extracted_genoF --make-grm-bin --out ${2}_extracted_genoF_grm

			#produce GRM matrices for each haplotype

			flines="$(< ${2}_extracted_hapF.fam wc -l)"
			mlines="$(< ${2}_extracted_hapM.fam wc -l)"
			plines="$(< ${2}_extracted_hapP.fam wc -l)"

			echo -e "\nProducing maternal-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapM1 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapF.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_hapM1_tmp

			mv ${2}_extracted_hapM1_tmp.bed ${2}_extracted_hapM1.bed
			mv ${2}_extracted_hapM1_tmp.bim ${2}_extracted_hapM1.bim
			mv ${2}_extracted_hapM1_tmp.fam ${2}_extracted_hapM1.fam
			mv ${2}_extracted_hapM1_tmp.nosex ${2}_extracted_hapM1.nosex
			mv ${2}_extracted_hapM1_tmp.log ${2}_extracted_hapM1.log
			
			echo -e "\nProducing maternal-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapM1 --make-grm-bin --out ${2}_extracted_hapM1_grm
			adjustGRM ${2}_extracted_hapM1_grm $flines 2

			echo -e "\nProducing paternal-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapP1 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapF.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_hapP1_tmp

			mv ${2}_extracted_hapP1_tmp.bed ${2}_extracted_hapP1.bed
			mv ${2}_extracted_hapP1_tmp.bim ${2}_extracted_hapP1.bim
			mv ${2}_extracted_hapP1_tmp.fam ${2}_extracted_hapP1.fam
			mv ${2}_extracted_hapP1_tmp.nosex ${2}_extracted_hapP1.nosex
			mv ${2}_extracted_hapP1_tmp.log ${2}_extracted_hapP1.log

			echo -e "\nProducing paternal-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapP1 --make-grm-bin --out ${2}_extracted_hapP1_grm
			adjustGRM ${2}_extracted_hapP1_grm $flines 2

			echo -e "\nProducing maternal non-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapM2 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapM.fam \
				--extract ${2}_snps_gt_${5}.txt  --indiv-sort f ${2}_extracted_hapM.fam --make-bed --out ${2}_extracted_hapM2_tmp

			mv ${2}_extracted_hapM2_tmp.bed ${2}_extracted_hapM2.bed
			mv ${2}_extracted_hapM2_tmp.bim ${2}_extracted_hapM2.bim
			mv ${2}_extracted_hapM2_tmp.fam ${2}_extracted_hapM2.fam
			mv ${2}_extracted_hapM2_tmp.nosex ${2}_extracted_hapM2.nosex
			mv ${2}_extracted_hapM2_tmp.log ${2}_extracted_hapM2.log

			echo -e "\nProducing maternal non-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapM2 --make-grm-bin --out ${2}_extracted_hapM2_grm
			adjustGRM ${2}_extracted_hapM2_grm $mlines 2

			echo -e "\nProducing paternal non-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapP2 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapP.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapP.fam --make-bed --out ${2}_extracted_hapP2_tmp

			mv ${2}_extracted_hapP2_tmp.bed ${2}_extracted_hapP2.bed
			mv ${2}_extracted_hapP2_tmp.bim ${2}_extracted_hapP2.bim
			mv ${2}_extracted_hapP2_tmp.fam ${2}_extracted_hapP2.fam
			mv ${2}_extracted_hapP2_tmp.nosex ${2}_extracted_hapP2.nosex
			mv ${2}_extracted_hapP2_tmp.log ${2}_extracted_hapP2.log

			echo -e "\nProducing paternal non-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapP2 --make-grm-bin --out ${2}_extracted_hapP2_grm
			adjustGRM ${2}_extracted_hapP2_grm $plines 2

			rm ${2}_extracted_hap.bim ${2}_extracted_hapF.fam ${2}_extracted_hapM.fam ${2}_extracted_hapP.fam
		elif [ ${4} = "pat-duos" ]
		then
			echo -e "\nCalculating GRMs from BEDs...\n"

			#produce standard GRM for genotypes
			#(use --keep-allele-order if want to keep ref/alt assignments)

			echo -e "\nProducing paternal-genotype BED above set MAF cutoff...\n"

			plink --vcf ${2}_extracted.vcf.gz --keep ${2}_extracted_hapP.tokeep \
				--maf ${5} --keep-allele-order --indiv-sort f ${2}_extracted_hapP.tokeep --make-bed --out ${2}_extracted_genoP

			awk '{print $1 ":" $4 ":" $6 ":" $5}' ${2}_extracted_genoP.bim > ${2}_snps_gt_${5}.txt

			modbim ${2}_extracted_genoM
			awk '{print $2}' ${2}_extracted_genoP > ${2}_extracted_geno_snps_list.txt

			echo -e "\nProducing paternal-genotype GRM...\n"

			plink --bfile ${2}_extracted_genoP --make-grm-bin --out ${2}_extracted_genoP_grm

			echo -e "\nProducing fetal-genotype BED above set MAF cutoff...\n"

			plink --vcf ${2}_extracted.vcf.gz --keep ${2}_extracted_hapF.fam \
				--keep-allele-order --make-bed --out ${2}_extracted_genoF_tmp

			modbim ${2}_extracted_genoF_tmp

			echo -e "\n"

			plink --bfile ${2}_extracted_genoF_tmp --extract ${2}_extracted_geno_snps_list.txt \
				--keep-allele-order --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_genoF

			echo -e "\nProducing fetal-genotype GRM...\n"

			plink --bfile ${2}_extracted_genoF --make-grm-bin --out ${2}_extracted_genoF_grm

			rm ${2}_extracted_genoF_tmp.*

			#produce GRM matrices for each haplotype

			flines="$(< ${2}_extracted_hapF.fam wc -l)"
			plines="$(< ${2}_extracted_hapP.fam wc -l)"

			echo -e "\nProducing maternal-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapM1 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapF.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_hapM1_tmp

			mv ${2}_extracted_hapM1_tmp.bed ${2}_extracted_hapM1.bed
			mv ${2}_extracted_hapM1_tmp.bim ${2}_extracted_hapM1.bim
			mv ${2}_extracted_hapM1_tmp.fam ${2}_extracted_hapM1.fam
			mv ${2}_extracted_hapM1_tmp.nosex ${2}_extracted_hapM1.nosex
			mv ${2}_extracted_hapM1_tmp.log ${2}_extracted_hapM1.log
			
			echo -e "\nProducing maternal-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapM1 --make-grm-bin --out ${2}_extracted_hapM1_grm
			adjustGRM ${2}_extracted_hapM1_grm $flines 2

			echo -e "\nProducing paternal-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapP1 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapF.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_hapP1_tmp

			mv ${2}_extracted_hapP1_tmp.bed ${2}_extracted_hapP1.bed
			mv ${2}_extracted_hapP1_tmp.bim ${2}_extracted_hapP1.bim
			mv ${2}_extracted_hapP1_tmp.fam ${2}_extracted_hapP1.fam
			mv ${2}_extracted_hapP1_tmp.nosex ${2}_extracted_hapP1.nosex
			mv ${2}_extracted_hapP1_tmp.log ${2}_extracted_hapP1.log

			echo -e "\nProducing paternal-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapP1 --make-grm-bin --out ${2}_extracted_hapP1_grm
			adjustGRM ${2}_extracted_hapP1_grm $flines 2

			echo -e "\nProducing paternal non-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapP2 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapP.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapP.fam --make-bed --out ${2}_extracted_hapP2_tmp

			mv ${2}_extracted_hapP2_tmp.bed ${2}_extracted_hapP2.bed
			mv ${2}_extracted_hapP2_tmp.bim ${2}_extracted_hapP2.bim
			mv ${2}_extracted_hapP2_tmp.fam ${2}_extracted_hapP2.fam
			mv ${2}_extracted_hapP2_tmp.nosex ${2}_extracted_hapP2.nosex
			mv ${2}_extracted_hapP2_tmp.log ${2}_extracted_hapP2.log

			echo -e "\nProducing paternal non-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapP2 --make-grm-bin --out ${2}_extracted_hapP2_grm
			adjustGRM ${2}_extracted_hapP2_grm $plines 2

			rm ${2}_extracted_hap.bim ${2}_extracted_hapF.fam ${2}_extracted_hapM.fam ${2}_extracted_hapP.fam
		elif [ ${4} = "mat-duos" ]
		then
			echo -e "\nCalculating GRMs from BEDs...\n"

			#produce standard GRM for genotypes
			#(use --keep-allele-order if want to keep ref/alt assignments)

			echo -e "\nProducing maternal-genotype BED above set MAF cutoff...\n"

			plink --vcf ${2}_extracted.vcf.gz --keep ${2}_extracted_hapM.tokeep \
				--maf ${5} --keep-allele-order --indiv-sort f ${2}_extracted_hapM.tokeep --make-bed --out ${2}_extracted_genoM

			awk '{print $1 ":" $4 ":" $6 ":" $5}' ${2}_extracted_genoM.bim > ${2}_snps_gt_${5}.txt

			modbim ${2}_extracted_genoM
			awk '{print $2}' ${2}_extracted_genoM > ${2}_extracted_geno_snps_list.txt

			echo -e "\nProducing maternal-genotype GRM...\n"

			plink --bfile ${2}_extracted_genoM --make-grm-bin --out ${2}_extracted_genoM_grm

			echo -e "\nProducing fetal-genotype BED above set MAF cutoff...\n"

			plink --vcf ${2}_extracted.vcf.gz --keep ${2}_extracted_hapF.fam \
				--keep-allele-order --make-bed --out ${2}_extracted_genoF_tmp

			modbim ${2}_extracted_genoF_tmp

			echo -e "\n"

			plink --bfile ${2}_extracted_genoF_tmp --extract ${2}_extracted_geno_snps_list.txt \
				--keep-allele-order --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_genoF

			echo -e "\nProducing fetal-genotype GRM...\n"

			plink --bfile ${2}_extracted_genoF --make-grm-bin --out ${2}_extracted_genoF_grm

			rm ${2}_extracted_genoF_tmp.*

			#produce GRM matrices for each haplotype

			flines="$(< ${2}_extracted_hapF.fam wc -l)"
			mlines="$(< ${2}_extracted_hapM.fam wc -l)"

			echo -e "\nProducing maternal-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapM1 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapF.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_hapM1_tmp

			mv ${2}_extracted_hapM1_tmp.bed ${2}_extracted_hapM1.bed
			mv ${2}_extracted_hapM1_tmp.bim ${2}_extracted_hapM1.bim
			mv ${2}_extracted_hapM1_tmp.fam ${2}_extracted_hapM1.fam
			mv ${2}_extracted_hapM1_tmp.nosex ${2}_extracted_hapM1.nosex
			mv ${2}_extracted_hapM1_tmp.log ${2}_extracted_hapM1.log
			
			echo -e "\nProducing maternal-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapM1 --make-grm-bin --out ${2}_extracted_hapM1_grm
			adjustGRM ${2}_extracted_hapM1_grm $flines 2

			echo -e "\nProducing paternal-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapP1 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapF.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapF.fam --make-bed --out ${2}_extracted_hapP1_tmp

			mv ${2}_extracted_hapP1_tmp.bed ${2}_extracted_hapP1.bed
			mv ${2}_extracted_hapP1_tmp.bim ${2}_extracted_hapP1.bim
			mv ${2}_extracted_hapP1_tmp.fam ${2}_extracted_hapP1.fam
			mv ${2}_extracted_hapP1_tmp.nosex ${2}_extracted_hapP1.nosex
			mv ${2}_extracted_hapP1_tmp.log ${2}_extracted_hapP1.log

			echo -e "\nProducing paternal-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapP1 --make-grm-bin --out ${2}_extracted_hapP1_grm
			adjustGRM ${2}_extracted_hapP1_grm $flines 2

			echo -e "Producing maternal non-transmitted BED above set MAF cutoff...\n"

			plink --bfile ${2}_extracted_hapM2 --bim ${2}_extracted_hap.bim --fam ${2}_extracted_hapM.fam \
				--extract ${2}_snps_gt_${5}.txt --indiv-sort f ${2}_extracted_hapM.fam --make-bed --out ${2}_extracted_hapM2_tmp

			mv ${2}_extracted_hapM2_tmp.bed ${2}_extracted_hapM2.bed
			mv ${2}_extracted_hapM2_tmp.bim ${2}_extracted_hapM2.bim
			mv ${2}_extracted_hapM2_tmp.fam ${2}_extracted_hapM2.fam
			mv ${2}_extracted_hapM2_tmp.nosex ${2}_extracted_hapM2.nosex
			mv ${2}_extracted_hapM2_tmp.log ${2}_extracted_hapM2.log

			echo -e "\nProducing maternal non-transmitted GRM...\n"

			plink --bfile ${2}_extracted_hapM2 --make-grm-bin --out ${2}_extracted_hapM2_grm
			adjustGRM ${2}_extracted_hapM2_grm $mlines 2

			rm ${2}_extracted_hap.bim ${2}_extracted_hapF.fam ${2}_extracted_hapM.fam ${2}_extracted_hapP.fam
		fi
	else
		echo -e "plink isn't present in ${1} or any default path\n"
		echo -e "splitted haplotype BED files and corresponding FAM & BIM files are saved in ${1}\n"
		exit 0
	fi
else
	echo -e "splitted haplotype BED files and corresponding FAM & BIM files are saved in ${1}\n"
	exit 0
fi
