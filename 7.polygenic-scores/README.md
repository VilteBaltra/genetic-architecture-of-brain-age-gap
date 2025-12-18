### Polygenic score scripts

This folder contains scripts for deriving polygenic scores in Generation R and runnning trio genetic analysis.  

UK Biobank–specific code for brain age estimation, GWAS, and PheWAS (including SBayesRC weights) is available at [pjawinski/enigma_brainage](https://github.com/pjawinski/enigma_brainage).

#### Scripts

- `run.example.sh` – example demonstrating how to run `pgs.predict.sh`  
- `pgs.predict.sh` – applies weights using PLINK2 
- `trio_GenR.R` – runs validation and trio genetic analysis (standardises polygenic scores and residualises them for genetic PCs, then runs multiple linear regressions where child, parent, and trio PGSs predict multiple outcomes)

#### Usage

1. Update file paths in `run.example.sh`  
2. Make the main script executable:

```bash
chmod +x pgs.predict.sh
```

3. Run to generate PGS scores

```bash
./run.example.sh 
```

4. Run trio genetic analysis (trio script by Esther Walton)
