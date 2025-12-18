prepare.files=function(df,cov){
  #col_cov=c('AGE','SEX','DX')
  col_cov=c('SUBJID','AGE','SEX','DX',grep("SITE", cov, value=T))
  
  col_MRI=c('ICV',
            'Mvent','Mthal','Mcaud','Mput','Mpal','Mhippo','Mamyg','Maccumb',
            'M_bankssts_thickavg','M_caudalanteriorcingulate_thickavg','M_caudalmiddlefrontal_thickavg',
            'M_cuneus_thickavg','M_entorhinal_thickavg','M_fusiform_thickavg','M_inferiorparietal_thickavg',
            'M_inferiortemporal_thickavg','M_isthmuscingulate_thickavg','M_lateraloccipital_thickavg','M_lateralorbitofrontal_thickavg',
            'M_lingual_thickavg','M_medialorbitofrontal_thickavg','M_middletemporal_thickavg','M_parahippocampal_thickavg',
            'M_paracentral_thickavg','M_parsopercularis_thickavg','M_parsorbitalis_thickavg','M_parstriangularis_thickavg',
            'M_pericalcarine_thickavg','M_postcentral_thickavg','M_posteriorcingulate_thickavg','M_precentral_thickavg','M_precuneus_thickavg',
            'M_rostralanteriorcingulate_thickavg','M_rostralmiddlefrontal_thickavg','M_superiorfrontal_thickavg','M_superiorparietal_thickavg',
            'M_superiortemporal_thickavg','M_supramarginal_thickavg','M_frontalpole_thickavg','M_temporalpole_thickavg','M_transversetemporal_thickavg',
            'M_insula_thickavg','M_bankssts_surfavg','M_caudalanteriorcingulate_surfavg','M_caudalmiddlefrontal_surfavg','M_cuneus_surfavg',
            'M_entorhinal_surfavg','M_fusiform_surfavg','M_inferiorparietal_surfavg','M_inferiortemporal_surfavg','M_isthmuscingulate_surfavg',
            'M_lateraloccipital_surfavg','M_lateralorbitofrontal_surfavg','M_lingual_surfavg','M_medialorbitofrontal_surfavg','M_middletemporal_surfavg',
            'M_parahippocampal_surfavg','M_paracentral_surfavg','M_parsopercularis_surfavg','M_parsorbitalis_surfavg','M_parstriangularis_surfavg',
            'M_pericalcarine_surfavg','M_postcentral_surfavg','M_posteriorcingulate_surfavg','M_precentral_surfavg','M_precuneus_surfavg',
            'M_rostralanteriorcingulate_surfavg','M_rostralmiddlefrontal_surfavg','M_superiorfrontal_surfavg','M_superiorparietal_surfavg',
            'M_superiortemporal_surfavg','M_supramarginal_surfavg','M_frontalpole_surfavg','M_temporalpole_surfavg','M_transversetemporal_surfavg','M_insula_surfavg', 'AGE', 'SUBJID')
  

  cols=c(col_cov,col_MRI)
  cols[!cols %in% names(data)]  
  
  data <- data[,cols]
  
  
  # Persons <10 or >75 years of age are excluded due to scarce datapoints in those age
  # extremes. Age variable is floored.
  data <- subset(data, AGE>9 & AGE<76)
  # data$AGE <- floor(data$AGE) # VB removed floor as this rounded age to lower number (removed decimals)
  data$AGE <- data$AGE
  
  # Calculate missings. If a subject has >10% missing values it is likely the segmentations of
  # correlated structures are not reliable. In this case, exclude participant.
  # The models use a standard scaler and median imputation for missing values.
  dim(data)
  num_nans <- rowSums(is.na(data[,col_MRI]))
  
  
  #new_data <- data[num_nans < 7.7,]
  new_data <- data[num_nans < 0.1*dim(data[,col_MRI])[2],] # VB: changed data[,col_MRI])[1] to data[,col_MRI])[2] to select 77 FS features (as row above)
  data <- new_data
  dim(data)
  
  # Partition the dataset into males and females and cases and controls.
  # Write separate .csv files for males and females
  # Upload to the PHOTON platform
  # . Female model ??? www.photon-ai.com/brainage/PZ46WUMIRP
  # . Male ??? www.photon-ai.com/brainage/WU6003FYD7
  
  males <- subset(data, SEX==0)
  males_raw <- males[,col_MRI]

  # Ensure the 'rawfiles' directory exists
  if (!dir.exists("rawfiles")) {
      dir.create("rawfiles")
  }
  # Save the CSV file in the 'rawfiles' directory
  write.csv(males_raw, file = "rawfiles/males_raw.csv", row.names = FALSE)
  
  females <- subset(data, SEX==1)
  females_raw <- females[,col_MRI]
  
  # Save the CSV file in the 'rawfiles' directory
  write.csv(females_raw, file = "rawfiles/females_raw.csv", row.names = FALSE)
 
  return(list(males=males,females=females))
}
