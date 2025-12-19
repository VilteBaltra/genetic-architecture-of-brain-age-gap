#!/usr/bin/env Rscript
# written by Philippe Jawinksi

# ==============================================
# === create phesant output table and figure ===
# ==============================================
#setwd("/Users/vb506/Documents/Documents-Nav-PC/manuscript-files-updated-July29/Philippe_output_enigma_brainage")

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=12) {
  stop(paste0('expected 12 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
trait = args[1] # trait="gap_gm"
phesantResults = args[2] # phesantResults="results/phewas/brainageFactor/phesant.output/results-combined.txt"
ncovs = as.numeric(args[3]) # ncovs=34
targetDir = args[4] # targetDir="results/phewas/brainageFactor/my_plots/"
imaging = args[5] # TRUE or FALSE, imaging = FALSE
multipleTesting = args[6] # 'fdr' or 'bonferroni' or 'both', multipleTesting = 'fdr'
ylim = as.numeric(args[7]) # upper y axis limit, ylim = 30
ybreaks = as.numeric(args[8]) # y axis breaks, ybreaks = 5
width = as.numeric(args[9]) # plot width (8.7 inch), width = 8.7
height = as.numeric(args[10]) # plot height (5.3 inch), height = 5.3
repel_nudge_y = as.numeric(args[11]) # plot height (5.3 inch), repel_nudge_y = 10
desat = args[12] # TRUE or FALSE, desat = FALSE

message(paste0('\n--- Settings ---',
        '\ntrait: ', trait,
        '\nphesantResults: ', phesantResults,
        '\nncovs: ', ncovs,
        '\ntargetDir: ', targetDir,
        '\nmultipleTesting: ', multipleTesting,
        '\nimaging: ', imaging,
        '\nylim: ', ylim,
        '\nybreaks: ', ybreaks,
        '\nwidth: ', width,
        '\nheight: ', height,
        '\nrepel_nudge_y: ', repel_nudge_y,
        '\ndesat: ', desat,'\n'))

# attach packages to current R session
for (pkg in c('data.table', 'dplyr', 'ggplot2', 'plotly', 'ggrepel')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# get phesant output 'combined results'
df = data.frame(fread(phesantResults, sep = '\t', header = T, quote = ""))

# remove imaging variables and those without pvalue
if (imaging == FALSE) {
  df = df[df$Cat2_Title!= "Imaging",]
}
df = df[!is.na(df$pvalue),]

# remove quotes in 'description'
df$description = gsub('"' , '', df$description, fixed = TRUE)

# calculate z-value (and se)
df$z = qnorm(df$p/2)*-sign(df$beta)
df$z[df$pvalue==0] = df[df$pvalue==0,] %>% 
                     mutate(z = beta/(abs(lower-upper)/(qnorm(1-(0.05/2))*2))) %>%
                     pull(., var = 'z')
df$se = abs(df$lower-df$upper)/(qnorm(1-(0.05/2))*2)
df$se[is.na(df$se)] = df$beta[is.na(df$se)]/df$z[is.na(df$se)]
df$se[df$beta==-999]=-999

# determine total sample size per phenotype
df$ntotal = df$n %>%
  gsub(pattern = ".*[(]", replacement = "") %>% 
  gsub(pattern = "[)].*", replacement = "") %>%
  gsub(pattern = ".*[/]", replacement = "") %>%
  as.numeric()

# calculate rho and absolute value of rho
df = df %>%
  mutate(rho = sqrt(z^2/(z^2 + (ntotal-2-ncovs)))*sign(beta)) %>%
  mutate(rho_abs = abs(rho))

# set pvalue to 3e-308 if 0
df$pvalue[df$pvalue == 0] = 3e-308

# make custom categories
df$custom_category = df$Cat3_Title
df$custom_category[df$Cat2_Title=='Imaging'] = 'Imaging'
df$custom_category[df$Cat2_Title=='Local environment'] = 'Lifestyle and environment'
df$custom_category[df$Cat2_Title=='Physical activity measurement'] = 'Physical activity'
df$custom_category[df$Cat2_Title=='Cognitive function online'] = 'Cognitive function'
df$custom_category[df$Cat2_Title=='Diet by 24-hour recall'] = 'Diet by 24-hour recall'
df$custom_category[df$Cat2_Title=='Digestive health'] = 'Digestive health'
df$custom_category[df$Cat2_Title=='Mental health'] = 'Mental health'
df$custom_category[df$Cat2_Title=='Work environment'] = 'Work environment'
df$custom_category[df$Cat2_Title=='Baseline characteristics'] = 'Sociodemographics'
df$custom_category[df$Cat2_Title=='Cognitive function'] = 'Cognitive function'
df$custom_category[df$Cat2_Title=='Physical measures'] = 'Physical measures'
df$custom_category[df$Cat3_Title=='Physical activity measurement'] = 'Physical activity'
df$custom_category[df$Cat3_Title=='Baseline characteristics'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Reception'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Employment'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Summary Administration'] = 'Hospital Inpatient - Administration'
df$custom_category[df$Cat3_Title=='Summary Operations'] = 'Operations'
df$custom_category[df$Cat3_Title=='Summary Maternity'] = 'Maternity and sex-specific factors'
df$custom_category[df$Cat3_Title=='Summary Psychiatric'] = 'Mental health'
df$custom_category[df$Cat3_Title=='Sex-specific factors'] = 'Maternity and sex-specific factors'
df$custom_category[df$Cat3_Title=='Health and medical history'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Medical conditions'] = 'Medical history and conditions'
df$custom_category[df$Cat2_Title=='Algorithmically-defined outcomes'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Psychosocial factors'] = 'Mental health'
df$custom_category[df$Cat3_Title=='Family history'] = 'Family history and early life factors'
df$custom_category[df$Cat3_Title=='Early life factors'] = 'Family history and early life factors'
df$custom_category[df$Cat3_Title=='Death Register'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=="Death register"] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Blood assays'] = 'Blood and urine assays'
df$custom_category[df$Cat3_Title=='Urine assays'] = 'Blood and urine assays'
df$custom_category[df$Cat3_Title=='Ongoing characteristics'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Process durations'] = 'Sociodemographics'

# calculate fdr
df$fdr = p.adjust(df$pvalue, method = 'BH')

# ----------------------------------
# --- draw phewas manhattan plot ---
# ----------------------------------

# create new data frame for plot
dfplot = df
dfplot = dfplot[order(dfplot$custom_category, dfplot$varName),]
ncategories = length(unique(dfplot$custom_category))

# set x axis positions for trait-phenotype associations
k = 0
for (category in unique(dfplot$custom_category)) {
    dfplot$pos[dfplot$custom_category == category] = ppoints(sum(dfplot$custom_category == category)) + k
    k = k +1
}

# set position for x axis breaks (centered)
axis.set = dfplot %>% 
  group_by(custom_category) %>% 
  summarize(center = mean(pos))

# set ylim significance threshold
if (multipleTesting == 'fdr') {
  sig = max(dfplot$pvalue[dfplot$fdr < 0.05])
} else if (multipleTesting == 'bonferroni' | multipleTesting == 'both') {
  sig = 0.05/nrow(dfplot)
}

# add variables for plotly tooltip
dfplot$variable = dfplot$description
dfplot$`field id` = dfplot$varName
dfplot$p = sprintf("%.6f", as.numeric(dfplot$pval))
dfplot$p[as.numeric(dfplot$pvalue) < 0.000001] = sprintf("%.2g", as.numeric(dfplot$pvalue[as.numeric(dfplot$pvalue) < 0.000001]))
dfplot$r = sprintf("%.6f", as.numeric(dfplot$rho))

# set pvalue of< ylim to ylim
dfplot$pvalue[dfplot$pvalue < 10^-ylim] = 10^-ylim

# get top hits per category for annotation
top = data.table(dfplot[order(dfplot$pvalue),], key = 'custom_category')
top = top[, head(.SD, 1), by=custom_category]
top = top[top$pvalue < sig,]

# set plot colors 
gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
desatfunc = function(cols, sat=0.8) {
  X = diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

# set.seed(1234)
# colcolors = gg_color_hue(ncategories)
# colcolors = colcolors[sample(ncategories, replace = F)]
# colcolors = c(tail(colcolors, 1), head(colcolors, -1)) # Shift colors one step to the right
colcolors <- ifelse(seq_len(ncategories) %% 2 == 1, '#ab4762', '#ee8a82')  # was  "#007C7C", "#40BFBF"
if (desat==TRUE) {
  fillcolors = desatfunc(colcolors, sat = 0.25)
} else {
  fillcolors = colcolors
}

# create manhattan plot
phewas_plot = ggplot2::ggplot(data = dfplot[dfplot$pvalue >= sig,], aes(textvariable = variable, textfield = `field id`, textn = n, textr = r, textp = p)) +
  geom_point(aes(y=pos, x=-log10(pvalue), fill=as.factor(custom_category)), alpha = 1, size = 1.5, stroke = 0, shape = 21) +
  geom_point(data = dfplot[dfplot$pvalue < sig,], aes(y=pos, x=-log10(pvalue), color=as.factor(custom_category)), fill = NA,  alpha = 1, size = 1.5, stroke = 0, shape = 16) +
  scale_fill_manual(values = fillcolors) +
  scale_color_manual(values = colcolors[unique(dfplot$custom_category) %in% dfplot$custom_category[dfplot$pvalue < sig]]) +
  scale_y_continuous(trans = "reverse", expand = c(0, 0), limits = c(max(axis.set$center) + 1, -0.5), label = unique(dfplot$custom_category), breaks = axis.set$center) + # label = axis.set$CHR % label = c(1:22, 'X', 'Y MT', '') label = c(1:18,'', 20, '', 22, 'X', 'Y MT', '')
  scale_x_continuous(expand = c(0,0), limits = c(0,ylim), breaks = seq(0,ylim,ybreaks)) +
  geom_vline(xintercept = -log10(sig), color = "black", linetype = "solid", size = 0.25) + 
  geom_segment(aes(y=min(axis.set$center),yend=max(axis.set$center),x=0,xend=0), colour = "black", size = 0.25) +
  geom_segment(aes(x = 0, xend = ylim, y = max(axis.set$center) + 1, yend = max(axis.set$center) + 1), size = 0.25) +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", size = 0.25),
    axis.ticks.length=unit(.15, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # axis.text.x = element_text(size = 9, vjust = 0.5, hjust = 0, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    # axis.text.y = element_text(angle = 0, size = 9, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
    # axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 5, b = 0, l = 5)),
    # axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    
    axis.text.x = element_text(size = 17, vjust = 0.5, hjust = 0, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(angle = 0, size = 17, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
    axis.title.y = element_text(size = 17, margin = margin(t = 0, r = 5, b = 0, l = 5)),
    axis.title.x = element_text(size = 17, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    
    # # VB chunk below instead of above for bigger font 
    # axis.text.x = element_text(size = 14, vjust = 0.5, hjust = 0, margin = margin(t = 5)),
    # axis.text.y = element_text(size = 14, vjust = 0.5, margin = margin(r = 3)),
    # axis.title.x = element_text(size = 18, margin = margin(t = 10)),
    # axis.title.y = element_text(size = 18, margin = margin(r = 5)),
    
    plot.margin=unit(c(0.25,1,0,-0.5),"cm"))

# # VB added
# phewas_plot = phewas_plot + theme_bw(base_size = 14)

# add fdr line
if (multipleTesting == 'both') {
  phewas_plot = phewas_plot + geom_vline(xintercept = -log10(max(dfplot$pvalue[dfplot$fdr < 0.05])), color = "black", linetype = "dashed", size = 0.25)
}

# draw ggplot
phewas_ggplot = phewas_plot +
    labs(y = "", x = expression("-log"[10]*"("*italic(p)*")")) # + 
    # geom_text_repel(data = top, aes(x=pos, y=-log10(pvalue), label=description),
    #                 max.time = 5,
    #                 max.iter = 100000,
    #                 box.padding = 0.5,
    #                 nudge_x = 0.5,
    #                 nudge_y = repel_nudge_y, # 10
    #                 direction = 'both',
    #                 size = 2.5,
    #                 angle = 0,
    #                 hjust = 0,
    #                 segment.size = 0.25,
    #                 segment.color = "grey10",
    #                 segment.square  = T,
    #                 segment.squareShape  = 1,
    #                 segment.angle     = -90,
    #                 segment.curvature = 1,
    #                 segment.ncp = 1)

# draw plotly
phewas_plotly = plotly::ggplotly(phewas_plot + labs(y = ""), tooltip = c("textvariable", "textfield", "textn", "textr", "textp"), dynamicTicks = FALSE, width = width*100, height = height*100)
phewas_plotly = phewas_plotly %>% config(displayModeBar = F) %>% layout(xaxis=list(fixedrange=TRUE)) %>% layout(yaxis=list(fixedrange=TRUE))

# ==========================================
# === saves figure and output data frame ===
# ==========================================

# save ggplot as png
message('(1/3) Saving ggplot as png.')
system(paste0('mkdir -p ', targetDir))
png(width = width, height = height, units = "in", res = 600, filename = paste0(targetDir, '/phesant_big_font.png'))
phewas_ggplot
invisible(dev.off())
system(paste0('chmod 770 ', targetDir,'/phesant_big_font.png'))

# save plotly as interactive html
# change height and width: phewas_plotly = phewas_plotly %>% layout(height = 800, width = 1200)
message('(2/3) Saving interactive plotly as html.')
htmlwidgets::saveWidget(phewas_plotly, paste0(targetDir, '/phesant_big_font.html'), title = paste0('PheWAS ', trait), selfcontained = TRUE)
system(paste0('chmod 770 ', targetDir,'/phesant_big_font.html'))
system(paste0('rm -rf ',targetDir,'/phesant_files_big_font'))
      
# export result table
message('(3/3) Saving plain text file with results summary.')
output = df[,c('description', 'varName', 'custom_category', 'Path', 'varType', 'resType', 'n',  'beta', 'lower', 'upper', 'se', 'pvalue', 'fdr', 'z', 'ntotal', 'rho')]
write.table(output, file = paste0(targetDir,'/phesant.summary_big_font.txt'), sep = '\t', quote = F, row.names = F, col.names = T)
system(paste0('chmod 770 ', targetDir,'/phesant.summary_big_font.txt'))
message('-- Creating PHESANT output table and figure completed. ---')
