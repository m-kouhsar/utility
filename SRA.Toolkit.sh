#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=72:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

if [ $# -eq 0 ]
then
  echo "Error: No accission supplied!"
  exit 1
elif [ $# -eq 1 ]
then
  Acc_list_file=$1
  Out_Dir=./
elif [ $# -gt 1 ]
then
  Acc_list_file=$1
  Out_Dir=$2
fi

while read Acc
do
  echo "Downloading $Acc"
  fastq-dump -A $Acc -O $Out_Dir
done < $Acc_list_file

echo "All samples were downloaded to $Out_Dir"
