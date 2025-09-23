#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=01:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#########################################################################

# liftover plink binary genotype data
# 

# ========= User input =========
INPUT_PREFIX=$1                         # e.g. yourdata (without .bed/.bim/.fam)
CHAIN_FILE=$2                           # UCSC Chain file (.chain or .chain.gz)
OUTPUT_PREFIX=$3                        # Output files prefix

# ========= Step 1. Extract BIM info =========
echo "[INFO] Extracting SNP positions from BIM..."
awk '{print "chr"$1, $4-1, $4, $2}' OFS="\t" ${INPUT_PREFIX}.bim > ${INPUT_PREFIX}.SNP.Pos.bed

# ========= Step 2. Run UCSC liftOver =========
echo "[INFO] Running liftOver..."
liftOver ${INPUT_PREFIX}.SNP.Pos.bed $CHAIN_FILE ${INPUT_PREFIX}.LiftOver.SNP.Pos.bed ${INPUT_PREFIX}.unmapped.bed

# ========= Step 3. Rebuild BIM file =========
echo "[INFO] Rebuilding BIM file with hg38 coordinates..."
awk 'NR==FNR {newpos[$4]=$2; chr[$4]=$1; next} 
     {if ($2 in newpos) {$1=chr[$2]; sub(/^chr/,"",$1); $4=newpos[$2]} print}' \
     ${INPUT_PREFIX}.LiftOver.SNP.Pos.bed ${INPUT_PREFIX}.bim > ${OUTPUT_PREFIX}.bim

# ========= Step 4. Copy BED and FAM =========
cp ${INPUT_PREFIX}.bed ${OUTPUT_PREFIX}.bed
cp ${INPUT_PREFIX}.fam ${OUTPUT_PREFIX}.fam

# ========= Step 5. Report =========
echo "[DONE] Liftover complete."
echo "Output files: ${OUTPUT_PREFIX}.bed/.bim/.fam"
echo "Unmapped SNPs saved in: ${INPUT_PREFIX}.unmapped.bed"

