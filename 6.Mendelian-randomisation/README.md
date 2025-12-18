### Mendelian Randomisation 

This folder contains code for **bidirectional Mendelian Randomization (MR)** analyses between the **BAG GWASs** (BAG factor + 6 individual BAG GWASs) and external traits.

---

#### Methods included:

- Main MR: Inverse-Variance Weighted
- Sensitivity MR: Weighted Median, MR-Egger  
- Cochran’s Q (heterogeneity)
- MR-Egger intercept (directional pleiotropy)
- Leave-one-out analyses
- Steiger tests and filtering
- MRlap (accounts for sample overlap, winner’s curse, and weak instruments)

All analyses are run **in both directions**:
- BAG → trait  
- Trait → BAG

#### Instrument Selection

- p < 5e-8 or p < 5e-6 if few genome-wide significant hits present
- LD clumping: r² < 0.001, 10 Mb window  
- Ambiguous palindromic SNPs removed

---

#### Notes

- All GWAS are of primarily European ancestry  
- Results reflect genetically proxied lifelong effects

