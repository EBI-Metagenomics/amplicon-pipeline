#!/usr/bin/env bash


######################################################################################
##
## FIRST SECTION:
## Here we define the requirements we need to run our job: CPUs, RAM, run time, ...
##
######################################################################################

#BSUB -J FETCH_RAW_READS_MULTIPLE # job name
#BSUB -W 04:00      # wall-clock time (hrs:mins)
#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q production   # queue
#BSUB -e /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/logs/asv-test/download_survey_datasets_error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/logs/asv-test/download_survey_datasets_output.%J  # output file name in which %J is replaced by the job ID



######################################################################################
##
## SECOND SECTION:
## Here we define the modules we need to be loaded in order to run our job.
## In this case, we need to load the R module, ... but first we make sure
## we start off with a sane, blank environment so we purge all modules
##
######################################################################################

module purge

mitload fetchtool

######################################################################################
##
## THIRD SECTION:
## Here we define the executable/binary/script we want LSF to run, our real data
## crunching section. Remember, in our case: we want to run an R script that
## calculates factorial of 8 and our R script is called "r.script".
##
## As you can see, we specify the FULL path to the script and we redirect the output
## to a file.
##
######################################################################################

datasets="/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/filtered_datasets.tsv"

{
    read # skip first line
    while IFS= read -r line
    do
        arr=($line)
        study="${arr[0]}"
        run="${arr[1]}"

        if [ -d "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/data/$study" ] 
        then
            echo "Study $study was already downloaded, skipping." 
        else
            sh fetch-reads-tool.sh -v -p $study -ru $run -d /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/data
            
            if [ $? == "0" ]
            then
                echo "${study} ${run}" >> /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/data/success.txt
            else
                echo "${study} ${run}" >> /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/data/failed.txt
            fi
        fi
        
    done
} < "$datasets"