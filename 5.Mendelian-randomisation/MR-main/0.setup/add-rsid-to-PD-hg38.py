import pandas as pd

# Load GP2 file
gp2 = pd.read_csv("GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt", sep="\t")

# Load VCF, skipping ## lines
vcf = pd.read_csv("Kaviar-160204-Public-hg38-trim.vcf", sep="\t", comment='#', header=None,
                  names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])

# Create merge key
gp2['key'] = gp2['chromosome'].astype(str) + ":" + gp2['base_pair_position'].astype(str)
vcf['key'] = vcf['CHROM'].astype(str) + ":" + vcf['POS'].astype(str)

# Merge
merged = gp2.merge(vcf[['key','ID']], on='key', how='left')

# Save
merged.to_csv("GP2_PD_with_Kaviar_IDs.txt", sep="\t", index=False)

# Remove duplicate IDs, keeping the first occurrence
merged_unique = merged.drop_duplicates(subset='ID', keep='first')
merged_unique.to_csv("GP2_PD_with_Kaviar_IDs_unique.txt.gz", sep="\t", index=False, compression='gzip')
