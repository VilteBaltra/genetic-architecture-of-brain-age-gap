# Manhattan plot
library(qqman)

# load data
df <- fread("brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz")
df <- as.data.frame(df)

png("my_manhattan_plot_20251124.png", width = 3400, height = 2100, res = 300)

# Manhattan plot with color
manhattan(df, chr="CHR", 
          bp="BP", 
          p="Pval_Estimate", 
          snp="SNP", 
          col=c('#ab4762', '#ee8a82'), 
          genomewideline = FALSE,
          suggestiveline = FALSE,
          cex = 1.2,            # point size
          cex.axis = 1.5,       # axis font size
          cex.lab = 1.4,
          ylim = c(0, 25))       


abline(h = -log10(5e-8), col = "black", lty = 1)        # genome-wide significance
abline(h = -log10(1e-5), col = "black", lty = 2)        # suggestive significance


dev.off()
