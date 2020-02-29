####################################################################################################################
####################################################################################################################
# Pre-process raw sequence data with BWA, Picard, and GATK.
# Author: Haiying Kong
# Last Modified: 28 October 2016
####################################################################################################################
####################################################################################################################
#!/bin/bash

####################################################################################################################
####################################################################################################################
# Run the software on all batches.
batches=(ILSE1669_A22 ILSE530X_B20 ILSE63XX_C28 GATC_A5 GATC_B11 ServiceXS_A8)

for batch in ${batches[@]}
do

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
dir_name=/home/kong/Haiying/Projects/PrimaryMelanoma/${batch}

# Change to working directory.
temp_dir=${dir_name}/temp
if [ ! -d "${temp_dir}" ]
  then mkdir ${temp_dir}
fi
cd ${temp_dir}

####################################################################################################################
# Select reference database for intervals.
case $batch in
  ILSE1669_A22)
    Intervals=Intervals_V5UTRs;;
  ILSE530X_B20)
    Intervals=Intervals_V5UTRs;;
  ILSE63XX_C28)
    Intervals=Intervals_V6UTRsr2;;
  GATC_A5)
    Intervals=Intervals_V6r2;;
  GATC_B11)
    Intervals=Intervals_V6r2;;
  ServiceXS_A8)
    Intervals=Intervals_V5;;
esac

####################################################################################################################
# Define directories to save log files, error files, and intermediate results.
log_dir=${dir_name}/log/BWA_Picard_GATK
if [ ! -d "${dir_name}/log" ]
  then mkdir ${dir_name}/log
fi
if [ -d "${log_dir}" ]
  then rm -r ${log_dir}
fi
mkdir ${log_dir}

error_dir=${dir_name}/error/BWA_Picard_GATK
if [ ! -d "${dir_name}/error" ]
  then mkdir ${dir_name}/error
fi
if [ -d "${error_dir}" ]
  then rm -r ${error_dir}
fi
mkdir ${error_dir}

Data_dir=${dir_name}/Data

if [ ! -d "${dir_name}/Lock" ]
  then mkdir ${dir_name}/Lock
fi

BWA_dir=${dir_name}/Lock/BWA_MEM
if [ -d "${BWA_dir}" ]
  then rm -r ${BWA_dir}
fi
mkdir ${BWA_dir}

Picard_dir=${dir_name}/Lock/Picard
if [ -d "${Picard_dir}" ]
  then rm -r ${Picard_dir}
fi
mkdir ${Picard_dir}

GATK_IBR_dir=${dir_name}/Lock/GATK_IndelBasedRealignment
if [ -d "${GATK_IBR_dir}" ]
  then rm -r ${GATK_IBR_dir}
fi
mkdir ${GATK_IBR_dir}

GATK_BQSR_dir=${dir_name}/Lock/GATK_BaseQualityScoreRecalibration
if [ -d "${GATK_BQSR_dir}" ]
  then rm -r ${GATK_BQSR_dir}
fi
mkdir ${GATK_BQSR_dir}

GATK_DOC_dir=${dir_name}/Lock/GATK_DepthOfCoverage
if [ -d "${GATK_DOC_dir}" ]
  then rm -r ${GATK_DOC_dir}
fi
mkdir ${GATK_DOC_dir}

####################################################################################################################
# Get data file names.
cd ${Data_dir}
data_files=($(ls *.fastq.gz))

####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Get sample names.
samples=($(echo ${data_files[@]%_R*.fastq.gz} | tr ' ' '\n' | sort -u | tr '\n' ' '))

####################################################################################################################
# Run pipeline on all samples.
####################################################################################################################
for sample in ${samples[@]}
do
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_BWA_Picard_GATK_${sample}    \
    -v sample=${sample},Data_dir=${Data_dir},BWA_dir=${BWA_dir},Picard_dir=${Picard_dir},GATK_IBR_dir=${GATK_IBR_dir},GATK_BQSR_dir=${GATK_BQSR_dir},GATK_DOC_dir=${GATK_DOC_dir},temp_dir=${temp_dir},Intervals=${Intervals}   \
    /home/kong/Haiying/Projects/PrimaryMelanoma/BatchCode_qsub/PreProcessing/BWA_Picard_GATK_NamedPipe_qsub.sh
done

####################################################################################################################
####################################################################################################################

done

####################################################################################################################
####################################################################################################################
