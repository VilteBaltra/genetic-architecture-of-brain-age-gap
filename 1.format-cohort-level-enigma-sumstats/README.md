### Scripts overview

`1.loop-script-format-Rsq.R` & `1.loop-script-format-info.R`
- Adds rsIDs if missing  
- Adds INFO/RSQ columns (differs per cohort, hence, separate scripts)  
- Performs basic quality control 
- Reads in helper functions from the `functions` folder

`2.info-filter-0.6.R`
- Loops through different cohorts  
- Filters SNPs with INFO or RSQ > 0.6

`3.manhattan-plot-per-cohort.R`
- Generates Manhattan plots for each cohort individually

---

**Notes**
- Ensure the `functions` folder is available before running scripts  
- Scripts should be executed in the order listed  
- R environment and required packages must be installed
