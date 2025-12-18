#!/bin/bash

# =====================
# === Calculate PGS ===
# =====================

weightFile="${1}"
bfile="${2}"        # merged .bed dataset (prefix only)
fileType="${3}"     # now always "bed"
idCol="${4}"
a1Col="${5}"
effectCol="${6}"
maf="${7}"
threads=${8}
outFile="${9}"

echo $'\n'"--- Calculate PGS | Settings ---"
echo "weightFile: ${weightFile}"
echo "bfile:      ${bfile}"
echo "idCol:      ${idCol}"
echo "a1Col:      ${a1Col}"
echo "effectCol:  ${effectCol}"
echo "maf:        ${maf}"
echo "threads:    ${threads}"
echo "outFile:    ${outFile}"$'\n'

# create target dir
targetDir="$(dirname "${outFile}")"
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# STEP 1 — Create weight file (remove SNPs with zero effect)
echo "Creating filtered weight file..."
awk -v cols="${idCol},${a1Col},${effectCol}" '
    BEGIN { ncols=split(cols,colnames,",") }
    NR==1 {
        for(i=1;i<=ncols;i++){
            for(j=1;j<=NF;j++){
                if ($j == colnames[i]) colnums[i]=j
            }
        }
    }
    $colnums[3] != 0 {
        print $colnums[1], $colnums[2], $colnums[3]
    }' <(gunzip -c "${weightFile}") > "${outFile}.weights"

# STEP 2 — Compute PGS using single merged BED dataset
echo "Running PLINK2 scoring on merged dataset..."

plink2 \
    --bfile "${bfile}" \
    --maf "${maf}" \
    --score "${outFile}.weights" 1 2 3 header cols=+scoresums \
    --threads "${threads}" \
    --out "${outFile}"

# Output will be:  outFile.sscore

# STEP 3 — Clean up
echo "Cleaning up."
rm -f "${outFile}.weights"

chmod 770 "${outFile}"*
echo "--- Completed: Calculate PGS ---"

