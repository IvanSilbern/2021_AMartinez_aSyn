local({
  
  library(data.table)
  
  if(!dir.exists("tables")) dir.create("tables")
  
  df1 <- fread("temp\\PeptQuant_average_Hippoc.txt")
  df2 <- fread("temp\\PeptQuant_average_Str_OB_SN.txt")
  lod <- fread("temp\\Linear_range_aSyn15N.txt")
  
  # Add region and enzyme
  df1[, Region := "Hippocampus"]
  df1[, Enzyme := "Trypsin"]
  df1[grepl("^K", Peptide), Enzyme := "LysN"]
  
  # rename mice
  df1[Condition == "HU-P3", Condition := "TM1"]
  df1[Condition == "HU-TR", Condition := "TM2"]
  df1[Condition == "P3",    Condition := "TP"]
  df1[Condition == "TR",    Condition := "D119"]
  
  # combine tables
  df <- rbind(df1, df2)
  
  # set significance level at 2 digits for the mean
  df[, aSyn_conc_mean := signif(aSyn_conc_mean, 2)]
  df[, aSyn_conc_sd   := signif(aSyn_conc_sd, 1)]
  
  # add starting point of the linear range (log2 Intensity)
  df <- merge(df, lod[, c("NAME", "linear_range")], by.x = "Peptide", by.y = "NAME", all.x = TRUE)
  
  # add aSyn specificities
  df[Peptide == "EGVVHGVATVAEK", aSyn_specificity := "TM1;TM2;D119"]
  df[Peptide == "EGVVHGVATVPEK", aSyn_specificity := "TM1"]
  df[Peptide == "EGVVHGVTTVAEK", aSyn_specificity := "WT"]
  df[Peptide == "EQVTNVGGAVVTGVTAVAQK", aSyn_specificity := "WT;TM1;TM2;D119"]
  df[Peptide == "EQVTNVGGAVVTGVTPVAQK", aSyn_specificity := "TP"]
  df[Peptide == "KNEEGAPQEGILEDMPVD", aSyn_specificity := "D119"]
  df[Peptide == "TVEGAGNIAAATGFVK", aSyn_specificity := "WT"]
  df[Peptide == "TVEGAGSIAAATGFVK", aSyn_specificity := "TM1;TM2;TP;D119"]
  df[Peptide == "KEGVVHGVATVAE", aSyn_specificity := "TM1;TM2;D119"]
  df[Peptide == "KEGVVHGVTTVAE", aSyn_specificity := "WT"]
  
  # remove carry over in WT
  df <- df[!(Peptide == "TVEGAGSIAAATGFVK" & Condition == "WT")]
  
  # rename columns
  names(df) <- c("Peptide", "Mouse", "log2_Intensity_mean", "log2_Intensity_sd",
                 "aSyn_fmol_ug_mean", "aSyn_fmol_ug_sd", "Region", "Enzyme", "log2_Intensity_LLOQ",
                 "aSyn_specificity")
  df <- df[, c("Region", "Enzyme", "Peptide", "Mouse", "aSyn_specificity",
               "log2_Intensity_mean", "log2_Intensity_sd",
               "log2_Intensity_LLOQ", "aSyn_fmol_ug_mean", "aSyn_fmol_ug_sd")]
  
  fwrite(df, "tables\\aSyn_quant.txt", sep = "\t")
  
  })