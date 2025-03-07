\# NGS Analysis Workflow

*\# Default Batch Script –* 

*\#\!/bin/bash* 

*\#SBATCH \--job-name= (--JOB NAME--)* 

*\#SBATCH \--mail-type=ALL* 

*\#SBATCH \--mail-user=\[email\]*

*\#SBATCH \--qos=privileged* 

*\#SBATCH \-c 2* 

*\#SBATCH \--mem 2048* 

*\#SBATCH \-t 2:00:00* 

*\#SBATCH \-o (--JOB NAME--)\_log.out* 

*\#SBTACH \-e (--JOB NAME--)\_log.err* 

*echo "job started"* 

*\--JOB GOES HERE--* 

*echo "job ended"* 

\# Log In to HPC

$ ssh \[email\] 

\# Create a new directory 

$ mkdir NGS

$ cd NGS

\# Copy project file into your own directory and unzip it 

$ scp /global/project/\[locationID\]/Project\_NGS/NGSbams.tar.gz 

/global/home/\[personalID\]/NGS/ 

$ scp /global/project/\[locationID\]/Project\_NGS/ref22.tar.gz 

/global/home/\[personalID\]/NGS/ 

$ tar \-xzf NGSbams.tar.gz 

$ tar \-xzf ref22.tar.gz 

\# Make another directory and copy other files into it 

$ cd NGSbams 

$ mkdir other\_files 

$ scp 

/global/home/\[personalID\]/NGS/home/\[locationID\]/seqshop/example/ref22/dbsnp\_135.b37 .chr22.vcf.gz /global/home/\[personalID\]/NGS/NGSbams/other\_files/ $ scp 

/global/home/\[personalID\]/NGS/home/\[locationID\]/seqshop/example/ref22/dbsnp\_135.b37 .chr22.vcf.gz.tbi /global/home/\[personalID\]/NGS/NGSbams/other\_files/  
$ scp 

/global/home/\[personalID\]/NGS/home/\[locationID\]/seqshop/example/ref22/human.g1k.v37 .chr22.dict /global/home/\[personalID\]/NGS/NGSbams/other\_files/ $ scp 

/global/home/\[personalID\]/NGS/home/\[locationID\]/seqshop/example/ref22/human.g1k.v37 .chr22.fa /global/home/\[personalID\]/NGS/NGSbams/other\_files/ 

$ scp 

/global/home/\[personalID\]/NGS/home/\[locationID\]/seqshop/example/ref22/human.g1k.v37 .chr22.fa.fai /global/home/\[personalID\]/NGS/NGSbams/other\_files/ $ scp /global/home/\[personalID\]/NGS/GATK/Module5\_rev/other\_files/snpEff/ /global/home/\[personalID\]/NGS/NGSbams/other\_files/ 

$ scp /global/project/\[locationID\]/SKATphenofile1.txt 

/global/home/\[personalID\]/NGS/NGSbams/ 

\# Load samtools 

$ module load samtools 

$ samtools flagstat HG01248.recal.bam 

\# Sort .bam files 

$ nano sort\_alignment.sh 

*\#batch script JOB GOES HERE* 

module load samtools 

for f in \*.bam; do samtools sort \-@ 4 \-m 5G \-T /tmp/${f}.sort \-o ${f}.sort.bam ${f}; done $ sbatch sort\_alignment.sh 

\# Remove duplicates 

$ nano removedups.sh 

*\#batch script JOB GOES HERE* 

module load samtools 

for f in \*.sort.bam; do filename="${f%%.\*}"; samtools rmdup \-s ${f} 

${filename}.sort.rmdup.bam; done 

$ sbatch removedups.sh 

\# Index rmdup.bam files 

$ nano index\_bam.sh 

*\#batch script JOB GOES HERE* 

module load samtools  
for f in \*.rmdup.bam; do samtools index ${f}; done 

$ sbatch index\_bam.sh 

\# Variant calling 

$ nano variantcall.sh 

*\#batch script JOB GOES HERE* 

module load java 

module load gatk/4.1.2.0 

for f in \*.sort.rmdup.bam; do java \-Xmx2g \-jar 

$EBROOTGATK/gatk-package-4.1.2.0-local.jar HaplotypeCaller \-R 

other\_files/human.g1k.v37.chr22.fa \-I ${f} \-O ${f}.vcf \-L 22:36000000-37000000; done $ sbatch variantcall.sh 

$ less HG00551.sort.rmdup.bam.vcf 

\# Filter .vcf files 

$ nano filtervariants.sh 

*\#batch script JOB GOES HERE* 

module load java 

module load gatk/4.1.2.0 

for f in \*.vcf; do filename="${f%%.\*}"; java \-jar 

$EBROOTGATK/gatk-package-4.1.2.0-local.jar VariantFiltration \-R 

other\_files/human.g1k.v37.chr22.fa \--variant ${f} \-O 

${filename}.sort.rmdup.bam.filter.vcf \--filter-expression "QD \< 2.0" \--filter-expression "FS \> 200.0" \--filter-expression "MQ \< 40.0" \--filter-name "QDFilter" \--filter-name "FSFilter" \--filter-name "MQFilter"; done 

$ sbatch filtervariants.sh 

$ less HG00551.sort.rmdup.bam.filter.vcf 

\# Annotate using dbSNP 

$ nano dbsnpann.sh 

*\#batch script JOB GOES HERE* 

module load java 

module load gatk/4.1.2.0  
for f in \*.filter.vcf; do filename="${f%%.\*}"; java \-Xmx2g \-jar 

$EBROOTGATK/gatk-package-4.1.2.0-local.jar VariantAnnotator \-R 

other\_files/human.g1k.v37.chr22.fa \--dbsnp other\_files/dbsnp\_135.b37.chr22.vcf.gz \--variant ${f} \-O ${filename}.sort.rmdup.bam.filter.dbsnp.vcf \-L 22:36000000-37000000; done 

$ sbatch dbsnpann.sh 

$ less HG00551.sort.rmdup.bam.filter.dbsnp.vcf 

\# Annotate using snpEff 

$ nano snpeff.sh 

*\#batch script JOB GOES HERE* 

*( \--mem 8000* for processing) 

module load snpeff/4.3t 

for f in \*.dbsnp.vcf; do filename="${f%%.\*}"; java \-jar other\_files/snpEff/snpEff.jar eff \-c other\_files/snpEff/snpEff.config \-v \-no-intergenic \-i vcf \-o vcf hg19 ${f} \> ${filename}.sort.rmdup.bam.filter.dbsnp.snpeff.vcf; done 

$ sbatch snpeff.sh 

$ less HG00551.sort.rmdup.bam.filter.dbsnp.snpeff.vcf 

\# Merge vcf files and filter by MAF 

$ nano mergevcfs.sh 

*\#batch script JOB GOES HERE* 

module load samtools 

module load vcftools 

module load tabix 

module load nixpkgs/16.09 

for file in $(ls \*.dbsnp.vcf); do bgzip ${file}; done 

for file in $(ls \*.vcf.gz); do tabix ${file}; done 

vcf-merge \*.vcf.gz \> merged\_all\_variants.vcfvcftools 

vcftools \--vcf merged\_all\_variants.vcfvcftools \--max-maf 0.05 \--recode \--out merged\_all\_variants\_MAF0.05 

$ sbatch mergevcfs.sh

\# Did all the processing up to this point, but decided to use the Skat file that was provided to us for our final analysis. We decided not to filter by MAF to keep more variants and do common rare variant analysis. 

$ scp /global/project/\[locationID\]/Project\_NGS/SKATngs\_merged\_62subjects.vcf /global/home/\[personalID\]/NGS/NGSbams/ 

\# SnpEff annotations 

$ sbatch SKATngs\_snpeff.sh 

*\#batch script JOB GOES HERE* 

java \-jar other\_files/snpEff/snpEff.jar eff \-c other\_files/snpEff/snpEff.config \-v \-no-intergenic \-i vcf \-o vcf hg19 SKATngs\_merged\_62subjects.vcf \> 

SKATngs\_merged\_62subjects.snpeff.vcf 

\# Convert to plink 

$sbatch vcftoplink.sh 

*\#batch script JOB GOES HERE* 

module load plink/1.9b\_4.1-x86\_64 

plink \--vcf SKATngs\_merged\_62subjects.snpeff.vcf \--pheno 

SKATphenofile1.txt \--pheno-name PHENO \--allow-no-sex \--make-bed \-out SKATngs\_merged\_62subjects 

\# Use results from SnpEff, create SKAT.SetID 

$ nano snpsift\_setid.sh 

*\#batch script JOB GOES HERE* 

java \-jar other\_files/snpEff/SnpSift.jar filter "(ANN\[\*\].IMPACT has 

'HIGH')|(ANN\[\*\].IMPACT has 'MODERATE')|(ANN\[\*\].IMPACT has 

'LOW')|(ANN\[\*\].EFFECT has 'MODIFIER')" SKATngs\_merged\_62subjects.snpeff.vcf \> functionalVariants.vcf 

java \-jar other\_files/snpEff/SnpSift.jar extractFields \-e "." functionalVariants.vcf "ANN\[0\].GENE" ID \> SKAT.SetID 

$ sbatch snpsift\_setid.sh

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\# R Code for SKAT Analysis \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\# library(SKAT) 

\# Create SSD and info files 

File.Bed \<- "./SKATngs\_merged\_62subjects.bed" 

File.Bim \<- "./SKATngs\_merged\_62subjects.bim" 

File.Fam \<- "./SKATngs\_merged\_62subjects.fam" 

File.SetID \<- "./SKAT.SetID" 

File.SSD \<- "./SKATngs\_merged\_62subjects.SSD" 

File.Info \<- "./SKATngs\_merged\_62subjects.Info" 

Generate\_SSD\_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info) 

\# Get phenotype file and subtract 1 to get 1's and 0's for the phenotype FAM \<- Read\_Plink\_FAM(File.Fam, Is.binary=FALSE) 

y \<- FAM$Phenotype \-1 

\# Open SSD file 

SSD.INFO \<- Open\_SSD(File.SSD, File.Info) 

\# Number of samples 

SSD.INFO$nSample 

\# Number of sets

SSD.INFO$nSets 

\# Set type to “D” for dichotomous 

obj \<- SKAT\_Null\_Model(y \~ 1, out\_type \= "D") 

\# Association test with phenotypes and SNP sets in the SSD file   
out \<- SKAT.SSD.All(SSD.INFO, obj) 

out 

\# Generate a table of results in the order of p values 

output.df \= out$results 

output.df \= output.df\[order(output.df$P.value),\] 

write.table(output.df, file="./SKAT\_SSD\_ALL.out.txt", col.names \= TRUE, row.names \= FALSE) output.df 

\# Association test for the combined effect of common and rare variants 

obj \<-  SKAT\_Null\_Model(y \~ 1, out\_type="D") 

out \<- SKAT\_CommonRare.SSD.All(SSD.INFO, ob) 

\# Write out results into a table 

output.df \= out$results 

output.df \= output.df\[order(output.df$P.value),\] 

write.table(output.df, file="./SKAT\_SSD\_ALL.out.txt", col.names \= TRUE, row.names \= FALSE) 

\# Close SSD file 

Close\_SSD()  
output.df