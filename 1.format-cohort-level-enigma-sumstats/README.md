This folder contains scripts that help to format GWAS summary
statistics.  
  
\# The first scripts add rsIDs if missing, add INFO/rsq column, and run
basic QC  
1.loop-script-format-Rsq.R and 1.loop-script-format-info.R \# reads in
functions from ‘functions’ folder  
  
\# The second script loops through different cohorts and filters for
info or Rsq \> 0.6  
2.info-filter-0.6.R  
  
\# The third script creates a Manhattan plot for each cohort
individually  
3.manhattan-plot-per-cohort.R  
