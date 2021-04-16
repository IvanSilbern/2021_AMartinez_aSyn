local({
  
  # analyze results from Striatum, Olfactory Bulb and Substantia Nigra
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  
  # spiked in amounts of aSyn per 100ug of total protein in ng and pmol
  # * for WT Biological replicates 14, 15, 16 only 100 ng were spiked in instead of 200 ng
  # respective corrections are made below
  
  hsyn_molw <- 14460
  hsyn_ng   <- 200
  hsyn_pmol <- 1000 * hsyn_ng / hsyn_molw 
  
  if(!dir.exists("temp")) dir.create("temp")
  
  ff <- list.files("Skyline_results", pattern = "\\.tsv")
  ff <- ff[grepl("Striatum_OlfactBulb_SubstNigr", ff)]
  
  df <- data.table()
  for(i in seq_along(ff)){
    
    df <- rbind(df, fread(paste0("Skyline_results/", ff[i]), check.names = TRUE))
  
  }
  
  # remove test samples with pure proteins
  df <- df[!grepl("_pure", df$Replicate)]
  
  # remove not identified transitions
  df <- df[!is.na(Peak.Rank)]
  
  # use only rank 1-6 transitions
  df <- df[Peak.Rank %in% 1:6]
  
  # sum transition intensities per Peptide / Replicate / Precursor.Charge / Isotope.Label.Type
  df_sum <- df[, list(Area = sum(Area)), by = c("Peptide", "Replicate", "Precursor.Charge", "Isotope.Label.Type")]
  
  # Extract sample information
  
  # technical (injection) replicate
  df_sum[, Tech_Replicate := str_match(Replicate, "_(\\d+)$")[, 2]]
  df_sum[is.na(Tech_Replicate), Tech_Replicate := "01"]
  
  # is a blank sample?
  df_sum[, blank := grepl("blank", Replicate)]
 
  # remove blanks
  df_sum <- df_sum[blank == FALSE]
  
  # Condition
  df_sum[, Region := str_match(Replicate, "^([^_]+)")[, 2]]
  df_sum[, Enzyme := str_match(Replicate, "^[^_]+_([^_]+)")[, 2]]
  df_sum[, Condition := str_match(Replicate, "^[^_]+_[^_]+_([^_]+)")[, 2]]
  df_sum[, Bio_Replicate := str_match(Replicate, "^[^_]+_[^_]+_[^_]+_([^_]+)_")[, 2]]
 
  # for EQVTNVGGAVVTGVTAVAQK and EQVTNVGGAVVTGVTPVAQK
  # keep only triply charged species
  
  df_sum <- df_sum[-which(Peptide %in% c("EQVTNVGGAVVTGVTAVAQK", "EQVTNVGGAVVTGVTPVAQK") & Precursor.Charge == 2)]
  
  # Average Technical replicates
  df_sum <- df_sum[, list(Area = mean(Area)), by = c("Region", "Enzyme", "Condition",
                                                     "Bio_Replicate",  "Peptide", "Precursor.Charge",
                                                     "Isotope.Label.Type")]
  
  # log2-transform
  df_sum[, Area := log2(Area)]
  
  # split in heavy and light peptides
  df_h <- df_sum[Isotope.Label.Type == "heavy"]
  df_l <- df_sum[Isotope.Label.Type == "light"]
  
  # derive normalization factors based on heavy peptides
  df_h[, norm_factor := 2^(Area - median(Area)), by = c("Peptide")]
  
  # heavy peptide that is used for normalization
  df_l[Peptide == "EGVVHGVATVAEK", hpept := "EGVVHGVATVAEK"]
  df_l[Peptide == "EGVVHGVTTVAEK", hpept := "EGVVHGVATVAEK"]
  df_l[Peptide == "EGVVHGVATVPEK", hpept := "EGVVHGVATVAEK"]
  df_l[Peptide == "TVEGAGNIAAATGFVK", hpept := "TVEGAGSIAAATGFVK"]
  df_l[Peptide == "TVEGAGSIAAATGFVK", hpept := "TVEGAGSIAAATGFVK"]
  df_l[Peptide == "EQVTNVGGAVVTGVTAVAQK", hpept := "EQVTNVGGAVVTGVTAVAQK"]
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK", hpept := "EQVTNVGGAVVTGVTAVAQK"]
  df_l[Peptide == "KNEEGAPQEGILEDMPVD", hpept := "KEGVVHGVATVAE"]
  
  # add heavy peptide areas
  names(df_h)[names(df_h) == "Area"] <- "hpept_area"
  df_l <- merge(df_l, 
                df_h[, c("Region", "Enzyme", "Condition", "Bio_Replicate", "Peptide", "Precursor.Charge", "hpept_area", "norm_factor")],
                by.x = c("Region", "Enzyme", "Condition", "Bio_Replicate", "hpept",   "Precursor.Charge"),
                by.y = c("Region", "Enzyme", "Condition", "Bio_Replicate", "Peptide", "Precursor.Charge"),
                all.x = TRUE)
  
  # log2 ratio light/heavy
  df_l[, log2FC_lh := Area - hpept_area]
  
  # calculate approx. aSyn concentration (only for peptides that have EXACT heavy-counter part)
  # in pmol per 100 ug of total protein 
  df_l[hpept == Peptide, aSyn_conc := 2^log2FC_lh * ..hsyn_pmol]
  
  # convert asyn concentration to fmol per 1 ug of total protein
  df_l[, aSyn_conc := 10*aSyn_conc]
  
  # normalize peptide intensities by the factor
  df_l[, Area_norm := log2(2^Area / norm_factor)]
  
  # Summarize data per Biological Replicate
  df_l[, Area_mean := mean(Area_norm), by = c("Region", "Enzyme", "Condition", "Peptide")]
  df_l[, Area_sd   := sd(Area_norm),   by = c("Region", "Enzyme", "Condition", "Peptide")]
  
  df_l[, aSyn_conc_mean := mean(aSyn_conc), by = c("Region", "Enzyme", "Condition", "Peptide")]
  df_l[, aSyn_conc_sd   := sd(aSyn_conc),   by = c("Region", "Enzyme", "Condition", "Peptide")]
  
  
  # indirectly estimate concentrations of "EQVTNVGGAVVTGVTPVAQK" and "KNEEGAPQEGILEDMPVD" in TM1 and TM2
  # set peptide concentrations equal to TVEGAGSIAAATGFVK concentration in TP an D119
  # use calculate peptide concentration in TP and D119 based on intensity fold changes
  
  # estimated aSyn concentration in TP based on the TVEGAGSIAAATGFVK peptide
  p3 <- df_l[Peptide == "TVEGAGSIAAATGFVK" & Condition == "TP"]
  
  # assume the same concentration for EQVTNVGGAVVTGVTPVAQK in TP
  temp <- merge(df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "TP",
                     -c("aSyn_conc", "aSyn_conc_mean", "aSyn_conc_sd")],
                p3[, c("Region", "Enzyme", "Condition", "Bio_Replicate",
                       "aSyn_conc", "aSyn_conc_mean", "aSyn_conc_sd")],
                by = c("Region", "Enzyme", "Condition", "Bio_Replicate"))
  
  # based on the intensity fold change calculate EQVTNVGGAVVTGVTPVAQK in TM1
  temp <- merge(temp[, c("Region", "Enzyme", "Condition", "Bio_Replicate",
                         "Area_norm", "Area_mean", "Area_sd",
                         "aSyn_conc", "aSyn_conc_mean", "aSyn_conc_sd")],
                df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "TM1",
                     c("Region", "Enzyme", "Bio_Replicate",
                       "Area_norm", "Area_mean", "Area_sd")],
                by = c("Region", "Enzyme", "Bio_Replicate"),
                suffixes = c(".TP", ".TM1"))
  
  temp[, fc := 2^(Area_norm.TM1 - Area_norm.TP)]
  temp[, aSyn_conc.TM1 := fc*aSyn_conc]
  
  temp[, aSyn_conc_mean.TP  := mean(aSyn_conc), by = c("Region", "Enzyme")]
  temp[, aSyn_conc_mean.TM1 := mean(aSyn_conc.TM1), by = c("Region", "Enzyme")]
  
  temp[, fc_mean := Area_mean.TM1 - Area_mean.TP]
  temp[, fc_sd   := sqrt(Area_sd.TM1^2 + Area_sd.TP^2)]
  
  temp[, aSyn_conc_sd.TM1 := 2^(fc_mean + fc_sd)*aSyn_conc_mean.TP - aSyn_conc_mean.TM1]
  
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "TM1", aSyn_conc      := temp$aSyn_conc.TM1]
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "TM1", aSyn_conc_mean := temp$aSyn_conc_mean.TM1]
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "TM1", aSyn_conc_sd   := temp$aSyn_conc_sd.TM1]
  
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "TP", aSyn_conc      := temp$aSyn_conc]
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "TP", aSyn_conc_mean := temp$aSyn_conc_mean.TP]
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "TP", aSyn_conc_sd   := temp$aSyn_conc_sd]
  
  # repeat the same for KNEEGAPQEGILEDMPVD
  
  # estimated aSyn concentration in TR based on the TVEGAGSIAAATGFVK peptide
  tr <- df_l[Peptide == "TVEGAGSIAAATGFVK" & Condition == "D119" & Enzyme == "Trypsin"]
  
  # assume the same concentration for EQVTNVGGAVVTGVTPVAQK in P3
  temp <- merge(df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "D119",
                     -c("aSyn_conc", "aSyn_conc_mean", "aSyn_conc_sd")],
                tr[, c("Region", "Condition", "Bio_Replicate",
                       "aSyn_conc", "aSyn_conc_mean", "aSyn_conc_sd")],
                by = c("Region", "Condition", "Bio_Replicate"))
  
  # based on the intensity fold change calculate EQVTNVGGAVVTGVTPVAQK in HU-P3
  temp <- merge(temp[, c("Region", "Enzyme", "Condition", "Bio_Replicate",
                         "Area_norm", "Area_mean", "Area_sd",
                         "aSyn_conc", "aSyn_conc_mean", "aSyn_conc_sd")],
                df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "TM2",
                     c("Region", "Enzyme", "Bio_Replicate",
                       "Area_norm", "Area_mean", "Area_sd")],
                by = c("Region", "Enzyme", "Bio_Replicate"),
                suffixes = c(".D119", ".TM2"))
  
  temp[, fc := 2^(Area_norm.TM2 - Area_norm.D119)]
  temp[, aSyn_conc.TM2 := fc*aSyn_conc]
  
  temp[, aSyn_conc_mean.D119  := mean(aSyn_conc), by = c("Region", "Enzyme")]
  temp[, aSyn_conc_mean.TM2 := mean(aSyn_conc.TM2), by = c("Region", "Enzyme")]
  
  temp[, fc_mean := Area_mean.TM2 - Area_mean.D119]
  temp[, fc_sd   := sqrt(Area_sd.TM2^2 + Area_sd.D119^2)]
  
  temp[, aSyn_conc_sd.TM2 := 2^(fc_mean + fc_sd)*aSyn_conc_mean.D119 - aSyn_conc_mean.TM2]
  
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "TM2", aSyn_conc      := temp$aSyn_conc.TM2]
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "TM2", aSyn_conc_mean := temp$aSyn_conc_mean.TM2]
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "TM2", aSyn_conc_sd   := temp$aSyn_conc_sd.TM2]
  
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "D119", aSyn_conc      := temp$aSyn_conc]
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "D119", aSyn_conc_mean := temp$aSyn_conc_mean.D119]
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "D119", aSyn_conc_sd   := temp$aSyn_conc_sd]
  
  # write result tables
  
  fwrite(df_l, "temp\\PeptQuant_Str_OB_SN.txt", sep = "\t")
  fwrite(df_l[, .SD[1], by = c("Peptide", "Region", "Enzyme", "Condition"), .SDcols = c("Area_mean",
                                                                                        "Area_sd",
                                                                                        "aSyn_conc_mean",
                                                                                        "aSyn_conc_sd")],
                              "temp\\PeptQuant_average_Str_OB_SN.txt",
         sep = "\t")
  
  })
  