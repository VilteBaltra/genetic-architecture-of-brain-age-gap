### ENIGMA brain age scripts

This folder provides scripts to derive Han et al. (2021; PMID: 32424236) brain age gap (BAG) and to perform a genome-wide association study (GWAS) on it. There are two main components: phenotypic scripts and genetic (GWAS) scripts.

---

1. Phenotypic scripts
- Folder: `ENIGMA-brainage-local`  
- Purpose: Derives Han's brain age gap phenotype.  
- Instructions: See the `instructions` folder within `ENIGMA-brainage-local`.  

---

2. Genetic scripts
- Folder: `ENIGMA-brainage-GWAS-local`  
- Purpose: Runs the BAGHan ENIGMA GWAS.  
- Instructions: Follow `0.BrainAgeGWAS_genetic_analysis_instructions_20241216.docx`.  

---

**Notes**
- Both folders are approximately 1 GB as they include Singularity containers.  
- For UK Biobank–specific brain age estimation and GWAS code, see: [pjawinski/enigma_brainage](https://github.com/pjawinski/enigma_brainage)
