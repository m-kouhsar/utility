#!/bin/bash
# PlinkBinaryLiftOver.sh
# Usage:
# ./PlinkBinaryLiftOver.sh input_prefix chain_file output_prefix [ref_fasta]
# Example:
# ./PlinkBinaryLiftOver.sh mydata hg19ToHg38.over.chain.gz mydata_hg38 hg38.fa

set -e

# --- Input Parameters ---
INPUT_PREFIX=$1      # Original PLINK prefix (e.g., mydata)
CHAIN_FILE=$2        # UCSC chain file (e.g., hg19ToHg38.over.chain.gz)
OUTPUT_PREFIX=$3     # New PLINK prefix (e.g., mydata_hg38)
REF_FASTA=$4         # Target reference genome fasta (for allele matching)

# --- Function to display usage ---
show_usage() {
    echo "Usage:"
    echo "./PlinkBinaryLiftOver.sh input_prefix chain_file output_prefix ref_fasta"
    echo "Example:"
    echo "./PlinkBinaryLiftOver.sh mydata hg19ToHg38.over.chain.gz mydata_hg38 hg38.fa"
}

# --- Check the number of arguments ---
# The script requires at least 3 arguments, and accepts an optional 4th.
if [ "$#" -lt 4 ]; then
    echo "[ERROR] Incorrect number of arguments provided." >&2
    echo ""
    show_usage
    exit 1
fi

# --- File & Directory Setup ---
OUT_DIR=$(dirname "$OUTPUT_PREFIX")
IN_FILE=$(basename "$INPUT_PREFIX")

LOV_IN="${OUT_DIR}/${IN_FILE}.SNP.txt"
LOV_OUT="${OUT_DIR}/${IN_FILE}.SNP.lifted.txt"
LOV_REJ="${OUT_DIR}/${IN_FILE}.SNP.rejected.txt"
NEW_BIM="${OUT_DIR}/${IN_FILE}.lifted.bim"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# --- Main Workflow ---

echo "[INFO] Step 1: Convert .bim to BED-like format for liftOver..."
# Improved AWK command: Checks if 'chr' prefix exists. Assumes chain file uses 'chr' notation.
# If your chain file uses '1' instead of 'chr1', you would need to reverse this logic.
awk 'BEGIN{OFS="\t"} {
    chrom = $1
    if (chrom !~ /^chr/) {
        # Handle MT -> chrM conversion if necessary
        if (chrom == "MT") {
            chrom = "chrM"
        } else {
            chrom = "chr"chrom
        }
    }
    print chrom, $4-1, $4, $2
}' "${INPUT_PREFIX}.bim" > "$LOV_IN"

echo "[INFO] Step 2: Running UCSC liftOver..."
liftOver "$LOV_IN" "$CHAIN_FILE" "$LOV_OUT" "$LOV_REJ"

if [ $? -ne 0 ]; then
    echo "[ERROR] liftOver command failed. Please check your inputs and chain file."
    exit 1
fi

echo "[INFO] Step 3: Remove rejected SNPs from original .bim..."
# Create a new .bim file containing only SNPs that were successfully lifted
awk 'NR==FNR{a[$4]; next} !($2 in a)' "$LOV_REJ" "${INPUT_PREFIX}.bim" > "$NEW_BIM"

echo "[INFO] Step 4: Update positions in .bim using lifted SNPs..."
# Update the coordinates in the new .bim file
awk 'NR==FNR{a[$4]=$2+1; next} {if($2 in a) $4=a[$2]; print}' "$LOV_OUT" "$NEW_BIM" > "${NEW_BIM}.tmp"
mv "${NEW_BIM}.tmp" "$NEW_BIM"

echo "[INFO] Step 5: Rebuild PLINK binary using lifted coordinates and checking the alleles with reference genome..."
# Create a list of SNPs to keep and a map file for their new positions
awk '{print $2}' "$NEW_BIM" > "${OUT_DIR}/lifted_snps.txt"
awk '{print $2, $4}' "$NEW_BIM" > "${OUT_DIR}/lifted_map.txt"

# Use plink2 to create the new fileset
plink2 --bfile "$INPUT_PREFIX" \
       --extract "${OUT_DIR}/lifted_snps.txt" \
       --update-map "${OUT_DIR}/lifted_map.txt" \
       --fa "$REF_FASTA" \
       --rm-dup force-first \
       --sort-vars \
       --make-pgen \
       --out "$OUTPUT_PREFIX"


plink2 \
    --pfile "${OUTPUT_PREFIX}" \
    --make-bed \
    --out "${OUTPUT_PREFIX}"


# --- Cleanup and Summary ---
echo "[INFO] Cleaning up intermediate files..."
rm "$LOV_IN"
rm "$NEW_BIM"
rm "${OUT_DIR}/lifted_snps.txt"
rm "${OUT_DIR}/lifted_map.txt"
rm "${OUTPUT_PREFIX}.pvar"
rm "${OUTPUT_PREFIX}.psam"
rm "${OUTPUT_PREFIX}.pgen"


echo ""
echo "[INFO] LiftOver complete!"
echo "----------------------------------------"
echo "Lifted SNPs: $(wc -l < "$LOV_OUT")"
echo "Rejected SNPs: $(( $(wc -l < "$LOV_REJ") / 2 ))"
echo "Log of lifted SNPs: ${LOV_OUT}"
echo "Log of rejected SNPs: ${LOV_REJ}"
echo "New PLINK binary prefix: $OUTPUT_PREFIX"
echo "----------------------------------------"
