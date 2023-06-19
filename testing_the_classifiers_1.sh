#!/bin/bash

#SBATCH --job-name=testing_classifiers
#SBATCH --output=/storage/users/cferguson/jobs/final_run/testing_classifiers_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --array=1-24
#SBATCH --mem-per-cpu=3GB
#SBATCH --cpus-per-task=10

############ software
# impoting modules that I'll need
module load apps/anaconda-4.7.12.tcl

# allowing the script to access my conda enviroments
eval "$(conda shell.bash hook)"

############ defining variables
# setting the base directory
BaseDir=/storage/users/cferguson/rotation_2_data/Data_for_classifier_assessment/sym_data

# defining witch simulated dataset to use
data=SNPS_and_homopolymers

# pulling out the accession/filenames
file_names=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${BaseDir}/../file_names.txt )

# grabbing the names of the folders
folder_names=$( ls ${BaseDir}/sym_data_1/SNPs )

# grabbing the simulated data filenames
sym_names=$( ls ${BaseDir} )

echo ${file_names}

########### kraken2
conda activate kraken2

echo kraken2

database=/storage/users/cferguson/Classification_database/kraken2

for y in ${sym_names[@]}; do
   report=/storage/users/cferguson/rotation_2_data/classified_data/${y}/${data}/kraken2/report
   out_reads=/storage/users/cferguson/rotation_2_data/classified_data/${y}/${data}/kraken2/reads

   for x in ${folder_names[@]}; do
      echo ${y}
      echo ${x}

      file1=${BaseDir}/${y}/${data}/${x}/${file_names}

      time kraken2 --db ${database} ${file1} --threads 20 --report ${report}/${x}/${file_names}.out --classified-out ${out_reads}/${x}/${file_names}.out

   done

done

conda deactivate

########### kaiju
conda activate kaiju

echo kaiju

database=/storage/users/cferguson/Classification_database/kaiju

for y in ${sym_names[@]}; do
   report=/storage/users/cferguson/rotation_2_data/classified_data/${y}/${data}/kaiju/report
   out_reads=/storage/users/cferguson/rotation_2_data/classified_data/${y}/${data}/kaiju/reads

   for x in ${folder_names[@]}; do
      echo ${y}
      echo ${x}

      file1=${BaseDir}/${y}/${data}/${x}/${file_names}

      time kaiju -z 20  -t ${database}/nodes.dmp -f ${database}/viruses/kaiju_db_viruses.fmi -i ${file1} -o ${out_reads}/${x}/${file_names}.out

   done

done

conda deactivate

########### centrifuge
conda activate centrifuge

echo centrifuge

database=/storage/users/cferguson/Classification_database/centrifuge

for y in ${sym_names[@]}; do
   report=/storage/users/cferguson/rotation_2_data/classified_data/${y}/${data}/centrifuge/report
   out_reads=/storage/users/cferguson/rotation_2_data/classified_data/${y}/${data}/centrifuge/reads

   for x in ${folder_names[@]}; do
      echo ${y}
      echo ${x}

      file1=${BaseDir}/${y}/${data}/${x}/${file_names}

      time centrifuge -x ${database}/abv -f ${file1} -S ${out_reads}/${x}/${file_names}.out --report-file ${report}/${x}/${file_names}.out -p 20

   done

done

conda deactivate


