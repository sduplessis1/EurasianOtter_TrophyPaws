# Bioinformatic analyses for Eurasian otter trophy paws

Scripts used to conduct bioinformatic analyses for manuscript (in review): 'Historical genomic variation of Eurasian otters (*Lutra lutra*) in Britain, from hunting trophies'.

Briefly, we used a dataset of modern, high coverage variation (n=45) to interrogate two low coverage, historical samples to contextualise historical variation of Eurasian otters in Britain. 


## Summary of 01_paleomix/paleomix.sh

Aim: Script to run Paleomix, with all the details and settinsg in AncientTemplate_illumina.yaml. 

``` 
#!/bin/bash

# activate python virtual environment
. /path/paleomix/py3-venv-paleomix/bin/activate

# load up required software
module load bwa/v0.7.17 
module load samtools/1.15 
module load R/4.0.0 
module load picard/2.22.2 
export PATH="/path/bin:$PATH"

# run command
paleomix bam run --jar-root /path/bin \
	--jar-root /trinity/shared/modulefiles/site-local \
	--jar-root /path \
	--log-level info --max-threads 128 --adapterremoval-max-threads 64 \
	--bwa-max-threads 128 --log-file log_ancient AncientTemplate_illumina.yaml
 
deactivate
 ```

## Summary of 04_dip2hap/dip2hap.sh

Aim: Takes the vcf of high coverage, modern data, and finds the major consensus for each SNP, converting the diploid data to haploid.

``` 
#!/bin/bash
module load bcftools/1.14

WDIR=/path/04_dip2hap
cd $WDIR

# Modern samples, diploid vcf
VCF_1=/path/mLutLut_renamed_autosomes_bisnps.vcf.gz

# Find major consensus for each SNP
bcftools +setGT $VCF_1 -- -t q -i 'GT="0/1" & (FMT/AD[:0] >= FMT/AD[:1])' -n 0 | \
	bcftools +setGT - -- -t q -i 'GT="0/1" & (FMT/AD[:0] < FMT/AD[:1])' -n c:'1/1' | \
	bcftools view -o mLutLut_renamed_autosomes_bisnps_major.vcf.gz -O z 

 ```

## Summary of 05_formatting/file_formatting.sh

Aim: Combing modern and historical samples. Specifically, random reads from each of the historical bam files are sampled at each variant site in the modern data.


``` 
#!/bin/bash

module load angsd/0.933
PLINK=/path/plink/plink
WDIR=/path/05_formatting
cd $WDIR

## Modern, haploid vcf ##
VCF_MOD=/path/04_dip2hap/mLutLut_renamed_autosomes_bisnps_major.vcf.gz
PLINK_OUT=modern_major_chromcodes

# File prep, convert vcf to plink format
$PLINK --vcf $VCF_MOD --make-bed --vcf-half-call missing --out $PLINK_OUT --allow-extra-chr
# re-write ped file and make bed with populations and sex of ids
$PLINK --bfile $PLINK_OUT --update-ids pop_ids --make-bed --out $PLINK_OUT --allow-extra-chr
$PLINK --bfile $PLINK_OUT --update-sex sex --make-bed --out $PLINK_OUT --allow-extra-chr

## Trophy paws, random reads ##
BAM1=/path/01_paleomix/LUT01.mLutLut1.2.bam
BAM2=/path/01_paleomix/LUT02.mLutLut1.2.bam
VCF_3=mLutLut_trophypaws

# Sample random read for every modern SNP
# list of chromosomes (based on ref panel)
cut -f1 $PLINK_OUT".bim" | sort | uniq > chrs.txt 
# list of positions and possible alleles
cut -f 1,4,5,6 $PLINK_OUT".bim" > sites_ref.txt
# index list of positions using angsd
angsd sites index sites_ref.txt
# list bam files we're adding
echo -e "${BAM1}\n${BAM2}" > bam_list
# randomly sample a read per site from the BAM file of our historical sample
angsd -dohaplocall 1 -bam bam_list -out $VCF_3 \
	-doCounts 1 -nThreads 4 -minMapQ 30 -minQ 20 \
	-doMajorMinor 3 -sites sites_ref.txt

# convert angsd output to plink
/path/angsd/misc/haploToPlink $VCF_3".haplo.gz" $VCF_3"_haploid"
# convert from tfam... to bed, bim, fam
$PLINK --tfile $VCF_3"_haploid" --allow-extra-chr --make-bed --out $VCF_3"_haploid"
# update fam, to add the historical samples
echo "LUT01 LUT01 0 0 0 1" > $VCF_3"_haploid.fam"
echo "LUT02 LUT02 0 0 0 1" >> $VCF_3"_haploid.fam"

# rename snps so that they match between PLINK/angsd files
awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' ${PLINK_OUT}.bim > ${PLINK_OUT}.bim_new
mv ${PLINK_OUT}.bim_new ${PLINK_OUT}.bim 

# rename chromosomes for admixture
sed -i 's/\bLR738403.1\b/1/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738404.1\b/2/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738405.1\b/3/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738406.1\b/4/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738407.1\b/5/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738408.1\b/6/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738409.1\b/7/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738410.1\b/8/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738411.1\b/9/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738412.1\b/10/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738413.1\b/11/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738414.1\b/12/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738415.1\b/13/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738416.1\b/14/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738417.1\b/15/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738418.1\b/16/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738419.1\b/17/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738420.1\b/18/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738421.1\b/X/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738422.1\b/Y/g' ${PLINK_OUT}.bim
sed -i 's/\bLR738403.1\b/1/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738404.1\b/2/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738405.1\b/3/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738406.1\b/4/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738407.1\b/5/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738408.1\b/6/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738409.1\b/7/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738410.1\b/8/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738411.1\b/9/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738412.1\b/10/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738413.1\b/11/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738414.1\b/12/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738415.1\b/13/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738416.1\b/14/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738417.1\b/15/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738418.1\b/16/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738419.1\b/17/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738420.1\b/18/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738421.1\b/X/g' ${VCF_3}_haploid.bim
sed -i 's/\bLR738422.1\b/Y/g' ${VCF_3}_haploid.bim

# only use SNP genotypes in trophy paws i.e. snp list from ${VCF_3}_haploid.bim
cut -f 2 ${VCF_3}_haploid.bim > snps_to_include
# subset these from both files
$PLINK --bfile $PLINK_OUT --extract snps_to_include --make-bed --out ${PLINK_OUT}_subset --allow-extra-chr
# want to subset out only biallelic snps (these were identified when I tried to merge the 2 plink files, for the PCA, but we're gonna use it here first and still merge them for pca later)
# how the multiallelic snps were identified:
TEST=merge_test
$PLINK --bfile ${PLINK_OUT}_subset --make-bed -bmerge $VCF_3"_haploid" \
        --allow-extra-chr --allow-no-sex --out $TEST
# says we have multiallelic sites which can't be merged
# multiallelic snps are given in this file:
# merge_test-merge.missnp
# so we can remove them from both inputs
$PLINK --bfile ${PLINK_OUT}_subset --exclude ${TEST}-merge.missnp \
	--make-bed --out ${PLINK_OUT}_subset_bi
$PLINK --bfile ${VCF_3}_haploid --exclude ${TEST}-merge.missnp \
	--make-bed --out ${VCF_3}_haploid_bi

# merge files to single plink (for pca)
MERGED=mLutLut_merged_n47
$PLINK --bfile ${PLINK_OUT}_subset_bi --bmerge ${VCF_3}_haploid_bi --make-bed \
        --out $MERGED --allow-no-sex


## Converting to vcf, and filtering 
module load bcftools/1.14
# can convert merged file into vcf if needed
$PLINK --bfile mLutLut_merged_n47 --recode vcf -out mLutLut_merged_n47
# see what variants are there (specifically transitions and transversions)
bcftools stats mLutLut_merged_n47.vcf > mLutLut_merged_n47.vcf_stats
# copy the header
bcftools view -h mLutLut_merged_n47.vcf > mLutLut_merged_n47_nots.vcf
# use awk to exclude (!) transitions, from the vcf, and append onto the header
awk -F '\t' '! (($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C") || $0 ~ /^#/ )' mLutLut_merged_n47.vcf >> mLutLut_merged_n47_nots.vcf
# check its worked
bcftools stats mLutLut_merged_n47_nots.vcf > mLutLut_merged_n47_nots.vcf_stats

# convert back to plink, as input files for admixture
$PLINK --vcf mLutLut_merged_n47_nots.vcf --make-bed --out mLutLut_merged_n47_nots
# split files into modern/trophy
$PLINK --bfile mLutLut_merged_n47_nots --keep trophy_ids.txt --make-bed --out mLutLut_merged_n47_nots_trophy
$PLINK --bfile mLutLut_merged_n47_nots --keep modern_ids.txt --make-bed --out mLutLut_merged_n47_nots_modern

 ```

## Summary of 06_admixture/admixture.sh

Aim: Run admixture on modern samples, and project on historical samples.

``` 
#!/bin/bash

WDIR=/path/06_admixture
cd $WDIR

IN_PATH=/path/05_formatting
PLINK_OUT=mLutLut_merged_n47_nots_modern
VCF_3=mLutLut_merged_n47_nots_trophy

# Running admixture in a loop
for K in `seq 1 8`;
do
	# run admixture on modern data
	/path/admixture/dist/admixture_linux-1.3.0/admixture \
		${IN_PATH}/${PLINK_OUT}.bed $K --cv -j4 
	# Use learned allele frequencies as (fixed) input to next step
	cp ${PLINK_OUT}.${K}.P ${VCF_3}.${K}.P.in
	# Run projection admixture
	/path/admixture/dist/admixture_linux-1.3.0/admixture \
		-P ${IN_PATH}/${VCF_3}.bed $K --cv -j4
	# Combine for plotting
	cat ${PLINK_OUT}.${K}.Q > combined_bi.${K}.Q
	cat ${VCF_3}.${K}.Q >> combined_bi.${K}.Q
done

cat ${IN_PATH}/${PLINK_OUT}.fam > combined.fam
cat ${IN_PATH}/${VCF_3}.fam >> combined.fam

 ```

## Summary of 07_pca/pca.sh

Aim: Run PCA using modern data, and project on historical samples. 

``` 
#!/bin/bash

module load R/4.0.0
PLINK=/path/plink/plink
SMARTPCA="/path/eigensoft/EIG-7.2.1/src/eigensrc/smartpca"
CONVERTF="/path/eigensoft/EIG-7.2.1/src/convertf"

WDIR=/path/07_pca
cd $WDIR

IN_PATH=/path/05_formatting
MERGED=mLutLut_merged_n47_nots

# run R code to adjust *bim input file
Rscript --vanilla --slave R_code1.R

# copy so all the plink files match
cp ${IN_PATH}/${MERGED}.bed ${MERGED}_R.bed
cp ${IN_PATH}/${MERGED}.fam ${MERGED}_R.fam
# need col 6 to be 1 for every individual to be included in analyses
awk '{print $1, $2, $3, $4, $5, "1"}' ${MERGED}_R.fam > ${MERGED}_R.fam_temp
mv ${MERGED}_R.fam_temp ${MERGED}_R.fam
# write parameter file
echo "genotypename: ${MERGED}_R.bed
snpname: ${MERGED}_R.bim
indivname: ${MERGED}_R.fam
outputformat: EIGENSTRAT
genotypeoutname: merged.eigenstratgeno
snpoutname: merged.snp
indivoutname: merged.ind
familynames: YES
pordercheck: NO" > merged.par
# convert to eigenstrat format
$CONVERTF -p merged.par
# edit the merged.ind file in R
Rscript --vanilla --slave R_code2.R
# list with all reference population/samples to estimate PCs
echo  "East
North
SWEng
Wales" > modern_estimatePCs.txt
## another parameters file, specifying the projection
echo "genotypename: merged.eigenstratgeno
snpname: merged.snp
indivname: merged.ind
evecoutname: merged_projected.evec
evaloutname: merged_projected.eval
familynames: YES
numoutevec: 4
numthreads: 4
pordercheck: NO
#poplistname: modern_estimatePCs.txt
lsqproject: YES" > modern_pca_project.par
## run smartpca
$SMARTPCA -p modern_pca_project.par

 ```

## Summary of 08_privatesnps/fullscript.sh 

Aim: Calculate the number of SNPs private to each stronghold population. 

``` 
#!/bin/bash
module load bcftools/1.14
module load vcftools/0.1.16
module load htslib/1.10.2

WDIR=/path/08_privatesnps
cd $WDIR

BASE=mLutLut_merged_n47_nots
VCF=/path/05_formatting/${BASE}.vcf

## 1 - Split main vcf into population level files ##
EAST=${BASE}_east.vcf.gz
SWENG=${BASE}_sweng.vcf.gz
NORTH=${BASE}_north.vcf.gz
WALES=${BASE}_wales.vcf.gz
LUT01=${BASE}_LUT01.vcf.gz
LUT02=${BASE}_LUT02.vcf.gz

bcftools view $VCF --samples-file pop_east -Oz -o $EAST
bcftools index $EAST

bcftools view $VCF --samples-file pop_sweng -Oz -o $SWENG
bcftools index $SWENG

bcftools view $VCF --samples-file pop_north -Oz -o $NORTH
bcftools index $NORTH

bcftools view $VCF --samples-file pop_wales -Oz -o $WALES
bcftools index $WALES

bcftools view $VCF --samples LUT01_LUT01 -Oz -o $LUT01
bcftools index $LUT01

bcftools view $VCF --samples LUT02_LUT02 -Oz -o $LUT02
bcftools index $LUT02


## 2 - Remove empty lines (vcf only contains lines for variants found in that pop) ##
EAST1=${BASE}_east_noemptylines.vcf
SWENG1=${BASE}_sweng_noemptylines.vcf
NORTH1=${BASE}_north_noemptylines.vcf
WALES1=${BASE}_wales_noemptylines.vcf
LUT011=${BASE}_LUT01_noemptylines.vcf
LUT021=${BASE}_LUT02_noemptylines.vcf

vcf-subset -e $EAST > $EAST1
bgzip $EAST1
tabix -p vcf ${EAST1}.gz

vcf-subset -e $SWENG > $SWENG1
bgzip $SWENG1
tabix -p vcf ${SWENG1}.gz

vcf-subset -e $WALES > $WALES1
bgzip $WALES1
tabix -p vcf ${WALES1}.gz

vcf-subset -e $NORTH > $NORTH1
bgzip $NORTH1
tabix -p vcf ${NORTH1}.gz

vcf-subset -e $LUT01 > $LUT011
bgzip $LUT011
tabix -p vcf ${LUT011}.gz

vcf-subset -e $LUT02 > $LUT021
bgzip $LUT021
tabix -p vcf ${LUT021}.gz


## 3 - vcf compare to compare snps found in any combination of input files ##
	# output written to jobid.out
vcf-compare ${EAST1}.gz ${SWENG1}.gz ${WALES1}.gz ${NORTH1}.gz ${LUT011}.gz ${LUT021}.gz

 ```

## Summary of 02_mitobim/LUT01/mitobim.sh

Aim: Use Mitobim to assemble the whole mitochondrial genome for each sample seperately. (Repeated for LUT02)

``` 
#!/bin/bash

module load samtools/1.17
module load fastp/v0.20

IND=LUT01
BAM_PATH=/path/01_paleomix/${IND}.mLutLut1.2.bam

# pull out the mitochondrial scaffold from the paleomix output
samtools view $BAM_PATH LR822067.1 -o ${IND}.mLutLut1.2.mt.bam
samtools index ${IND}.mLutLut1.2.mt.bam

samtools sort -n ${IND}.mLutLut1.2.mt.bam | \
	samtools fastq -0 ${IND}.mLutLut1.2.mt.0.fq.gz

gunzip ${IND}.mLutLut1.2.mt.0.fq.gz
mv ${IND}.mLutLut1.2.mt.0.fq ${IND}.mLutLut1.2.mt.0.fastq

## Write mito genome ##
grep -A 10000 ^">LR822067.1" \
	/path/ref/GCA_902655055.2_mLutLut1.2_genomic.fna > \
	GCA_902655055.2_mLutLut1.2_genomic_mt.fa

# Write configuration file
echo -e "\n#manifest file for basic mapping assembly with illumina data using MIRA 4

project = $IND

job=genome,mapping,accurate

parameters = -DI:trt=/tmp/$IND -NW:mrnl=0 -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no 

readgroup
is_reference
data = GCA_902655055.2_mLutLut1.2_genomic_mt.fa
strain = $IND

readgroup = reads
data = ${IND}.mLutLut1.2.mt.0.fastq
technology = solexa
strain = $IND\n" > manifest_$IND.conf

## Run mira ##
mkdir -p /tmp/$IND
export LC_ALL=C
/path/mitobim_version/MITObim-master/docker/external_software/mira_4.0.2/mira manifest_$IND.conf
rm -rf /tmp/$IND

## Run MITObim.pl 
/path/mitobim_version/MITObim-master/MITObim.pl \
	--mirapath /path/mitobim_version/MITObim-master/docker/external_software/mira_4.0.2/ \
	-start 1 -end 10 \
	-sample $IND -ref GCA_902655055.2_mLutLut1.2_genomic_mt.fa \
       	-readpool ${IND}.mLutLut1.2.mt.0.fastq \
	-maf ${IND}_assembly/${IND}_d_results/${IND}_out.maf \
	--NFS_warn_only &> log

 ```

## Summary of 02_mitobim/mafft_align/mafft.sh

Aim: Aligning whole mitochondrial genome sequences using mafft. 

``` 
#!/bin/bash

cd /path/02_mitobim/mafft_align

# copy final fasta file from mitobim files
cp /path/02_mitobim/LUT01/iteration2/LUT01-GCA_902655055.2_mLutLut1.2_genomic_mt.fa-it2_noIUPAC.fasta .
cp /path/02_mitobim/LUT02/iteration2/LUT02-GCA_902655055.2_mLutLut1.2_genomic_mt.fa-it2_noIUPAC.fasta .
# rename files
mv LUT01-GCA_902655055.2_mLutLut1.2_genomic_mt.fa-it2_noIUPAC.fasta LUT01_mitobim.fasta
mv LUT02-GCA_902655055.2_mLutLut1.2_genomic_mt.fa-it2_noIUPAC.fasta LUT02_mitobim.fasta
# manually change header of each fasta sequence '>LUT01_mitobim' (vi ...)

# combine context sequences with newly generated sequences
cat UK_seq.fasta >> combined_sequences.fasta
cat LUT01_mitobim.fasta >> combined_sequences.fasta
cat LUT02_mitobim.fasta >> combined_sequences.fasta

# run through mafft to align
module load mafft-7.481-gcc-8.5.0-b6srufn
mafft --auto combined_sequences.fasta > combined_sequences.align

 ```

