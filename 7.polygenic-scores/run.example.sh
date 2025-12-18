#!/bin/bash

# requirements:
#   - plink2 must be installed and accessible as 'plink2'

## to do: remember to also make pgs script executable with chmod +x pgs.predict.sh

# example script for computing polygenic scores using a merged .bed dataset
bfile="/Users/vb506/Documents/SBayesRC-scripts-and-weights/GenR/QC-bfiles/your_data_all_chr_merged"   # <-- path to merged individual level data in .bed/.bim/.fam format (without extension)
idCol="Name"
a1Col="A1"
effectCol="A1Effect"
maf=0.00
threads=50

for gwas in brainageFactor han jawinski kaufmann leonardsen smith wen; do
    outFile="/Users/vb506/Documents/SBayesRC-scripts-and-weights/GenR/pgs/pgs.sbayesrc.${gwas}"
    weightFile="/Users/vb506/Documents/SBayesRC-scripts-and-weights/GenR/pgsweights/sbayesrc.${gwas}.snpRes.gz"

    ./pgs.predict.sh \
        "${weightFile}" \
        "${bfile}" \
        "bed" \
        "${idCol}" \
        "${a1Col}" \
        "${effectCol}" \
        "${maf}" \
        "${threads}" \
        "${outFile}"
done

