####################################################################################################################
####################################################################################################################
# Pre-process raw sequence data with BWA, Picard, and GATK -- qsub script.
# Author: Haiying Kong
# Last Modified: 7 November 2016
####################################################################################################################
####################################################################################################################
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=55GB
#PBS -l walltime=200:00:00

####################################################################################################################
# Decide number of thread.
n_thread=8

####################################################################################################################
####################################################################################################################
# Software:
Picard=/home/kong/lib/Picard-2.6.0/picard.jar
GATK=/home/kong/lib/GATK-3.6/GenomeAnalysisTK.jar

# Reference databases:
GATK_hg=/home/kong/Reference/GATK/Genome/human_g1k_v37.fasta
GATK_localrealignment_1000G=/home/kong/Reference/GATK/LocalRealignment/1000G_phase1.indels.b37.vcf
GATK_localrealignment_Mills=/home/kong/Reference/GATK/LocalRealignment/Mills_and_1000G_gold_standard.indels.b37.vcf
GATK_dbSNP138=/home/kong/Reference/GATK/SNP/dbsnp_138.b37.vcf
Intervals_V5=/home/kong/Reference/Agilent/SureSelectHumanAllExon_V5/S04380110_Covered_nochr.bed
Intervals_V5UTRs=/home/kong/Reference/Agilent/SureSelectHumanAllExon_V5UTRs/S04380219_Covered_nochr.bed
Intervals_V6r2=/home/kong/Reference/Agilent/SureSelectHumanAllExon_V6r2/S07604514_Covered_nochr.bed
Intervals_V6UTRsr2=/home/kong/Reference/Agilent/SureSelectHumanAllExon_V6UTRsr2/S07604624_Covered_nochr.bed

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Run BWA, Picard and GATK.

#########################################
# BWA-MEM:
mkfifo -m 0644 ${BWA_dir}/NamedPipe_${sample}.sam
bwa mem $GATK_hg ${Data_dir}/${sample}_R1.fastq.gz ${Data_dir}/${sample}_R2.fastq.gz -t ${n_thread} > ${BWA_dir}/NamedPipe_${sample}.sam &

#########################################
# Picard:
java -Xmx50G -Djava.io.tmpdir=tmp -jar ${Picard} SortSam     \
  INPUT=${BWA_dir}/NamedPipe_${sample}.sam OUTPUT=${Picard_dir}/${sample}.sorted.bam    \
  SORT_ORDER=coordinate CREATE_INDEX=true

rm ${BWA_dir}/NamedPipe_${sample}.sam

mkfifo -m 0644 ${Picard_dir}/NamedPipe_${sample}.dedupped.sam
java -Xmx50G -Djava.io.tmpdir=tmp -jar ${Picard} MarkDuplicatesWithMateCigar     \
  INPUT=${Picard_dir}/${sample}.sorted.bam OUTPUT=${Picard_dir}/NamedPipe_${sample}.dedupped.sam     \
  METRICS_FILE=${Picard_dir}/${sample}.metrics.txt    \
  REMOVE_DUPLICATES=false ASSUME_SORTED=true MINIMUM_DISTANCE=500 &

IFS='_' read -a RG <<< "$sample"

java -Xmx50G -Djava.io.tmpdir=tmp -jar ${Picard} AddOrReplaceReadGroups     \
  I=${Picard_dir}/NamedPipe_${sample}.dedupped.sam O=${Picard_dir}/${sample}.readgroups.bam     \
  RGLB=$sample RGPL=illumina RGPU=$RG[1] RGSM=$sample CREATE_INDEX=true

rm ${Picard_dir}/${sample}.sorted.bam
rm ${Picard_dir}/${sample}.sorted.bai
rm ${Picard_dir}/NamedPipe_${sample}.dedupped.sam
rm ${Picard_dir}/${sample}.metrics.txt

#########################################
# GATK_IBR:
java -Xmx50G -Djava.io.tmpdir=tmp -jar ${GATK} -T RealignerTargetCreator     \
  -I ${Picard_dir}/${sample}.readgroups.bam -o ${GATK_IBR_dir}/${sample}.realigner.intervals     \
  -R $GATK_hg -known $GATK_localrealignment_1000G -known $GATK_localrealignment_Mills     \
  --disable_auto_index_creation_and_locking_when_reading_rods
java -Xmx50G -Djava.io.tmpdir=tmp -jar ${GATK} -T IndelRealigner    \
  -I ${Picard_dir}/${sample}.readgroups.bam -o ${GATK_IBR_dir}/${sample}.realigned.bam    \
  -R $GATK_hg -known $GATK_localrealignment_1000G -known $GATK_localrealignment_Mills    \
  -targetIntervals ${GATK_IBR_dir}/${sample}.realigner.intervals    \
  --disable_auto_index_creation_and_locking_when_reading_rods

rm ${Picard_dir}/${sample}.readgroups.bam
rm ${Picard_dir}/${sample}.readgroups.bai
rm ${GATK_IBR_dir}/${sample}.realigner.intervals

#########################################
# GATK_BQSR:
java -Xmx50G -Djava.io.tmpdir=tmp -jar ${GATK} -T BaseRecalibrator -nct ${n_thread}     \
  -I ${GATK_IBR_dir}/${sample}.realigned.bam -o ${GATK_BQSR_dir}/${sample}.recal.table    \
  -R $GATK_hg -knownSites $GATK_localrealignment_1000G     \
  -knownSites $GATK_localrealignment_Mills -knownSites $GATK_dbSNP138    \
  --disable_auto_index_creation_and_locking_when_reading_rods
java -Xmx50G -Djava.io.tmpdir=tmp -jar ${GATK} -T PrintReads -nct ${n_thread}     \
  -I ${GATK_IBR_dir}/${sample}.realigned.bam -o ${GATK_BQSR_dir}/${sample}.recal.bam    \
  -R $GATK_hg -BQSR ${GATK_BQSR_dir}/${sample}.recal.table    \
  --disable_auto_index_creation_and_locking_when_reading_rods

rm ${GATK_IBR_dir}/${sample}.realigned.bam
rm ${GATK_IBR_dir}/${sample}.realigned.bai
rm ${GATK_BQSR_dir}/${sample}.recal.table

#########################################
# GATK_DOC:
java -Xmx50G -Djava.io.tmpdir=tmp -jar ${GATK} -T DepthOfCoverage    \
  -I ${GATK_BQSR_dir}/${sample}.recal.bam -o ${GATK_DOC_dir}/${sample}    \
  -R $GATK_hg -L ${!Intervals}     \
  -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0    \
  --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 --includeRefNSites --countType COUNT_FRAGMENTS    \
  --disable_auto_index_creation_and_locking_when_reading_rods

####################################################################################################################
####################################################################################################################
