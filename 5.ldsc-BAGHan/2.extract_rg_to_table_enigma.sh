
#!/bin/bash

# This script extracts rg values to a table for each file formatted as brainageFactor-cor-*.log 

cd /campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/genetic-corr/output2/

# output header
echo -e "Trait\tGeneticCorrelation\tSE\tZscore\tPvalue" > results.tsv

# loop through each file
for file in enigma-cor-*.log; do
  # Extract trait name (e.g., adhd, asd, etc.)
  trait=$(echo "$file" | sed -E 's/enigma-cor-(.*)\.log/\1/')

  # extract values
  gc_line=$(grep "Genetic Correlation:" "$file")
  z_line=$(grep "Z-score:" "$file")
  p_line=$(grep "^P:" "$file")

  # parse values
  gc=$(echo "$gc_line" | awk '{print $3}')
  se=$(echo "$gc_line" | sed -E 's/.*\((.*)\)/\1/')
  z=$(echo "$z_line" | awk '{print $2}')
  p=$(echo "$p_line" | awk '{print $2}')

  # append to file
  echo -e "${trait}\t${gc}\t${se}\t${z}\t${p}" >> results.tsv
done
