#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=10:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=0-19

#########################################################################
InputDir=/lustre/projects/Research_Project-191391/Project_10396/11_fastp_trimmed
OutDir=/lustre/projects/Research_Project-191391/Morteza/miRNA/Results/Project.11008.V0304.NorCog/flash.10396
Thread=15
min_overlap=10
max_overlap=150

#########################################################################
Samples=(${InputDir}/*R1*)

Num_samp=${#Samples[@]}

denom_2=$(( SLURM_ARRAY_TASK_COUNT / 2 ))
window_size=$(( ( Num_samp + denom_2 ) / SLURM_ARRAY_TASK_COUNT ))

lower=$(( SLURM_ARRAY_TASK_ID * window_size ))

Samples_1=(${Samples[@]:${lower}:${window_size}})

if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]
then
    Samples_1=(${Samples[@]:$lower})
fi

echo Input directory: $InputDir
echo Output directory: $OutDir
echo Total number of samples: $Num_samp
echo Start array index: $SLURM_ARRAY_TASK_MIN
echo End array index : $SLURM_ARRAY_TASK_MAX
echo numer of arrays: $SLURM_ARRAY_TASK_COUNT
echo current array index: $SLURM_ARRAY_TASK_ID
echo Number of samples in current array: ${#Samples_1[@]}

###########################################################################
if [ ${#Samples_1[@]} != 0 ]
then
    mkdir -p $OutDir
    cd $OutDir
    j=0
    for name in ${Samples_1[@]}
    do
        j=$((j+1))
        R1=$name
        R2=${name/R1/R2}
        name1=$( basename $R1)
        name1=${name1/R1/R12}
        echo "**********************************************************************************************************"
        echo "                 Merging Sample $j: $name1"
        echo "**********************************************************************************************************"
        flash -m $min_overlap -M $max_overlap $R1 $R2 -o ${name1}.flash -t $Thread 2>&1 | tee ${name1}.flash.log
    done
else
    echo "There is no samples in the current array!"
fi


