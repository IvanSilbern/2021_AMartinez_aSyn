local({
  
  # analyze results from Hippocampus
  
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
  
  tryp <- fread("Skyline_results\\20180403_Hippocamp_aSyn_SRM_trypsin_TransitionResults.tsv", check.names = TRUE)
  lysn <- fread("Skyline_results\\20180522_Hippocamp_aSyn_SRM_LysN_TransitionResults.tsv", check.names = TRUE)
  
  # combine tables
  df <- rbind(tryp, lysn)
  
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
  df_sum[, blank := grepl("_bl$", Replicate)]
  
  # Condition
  df_sum[, Condition := str_match(Replicate, "^\\d+_([^_]+)")[, 2]]
  df_sum[grepl("_LysN_", Replicate), Condition := str_match(Replicate, "^([^_]+)")[, 2]]
  
  # Biological Replicate
  df_sum[, Bio_Replicate := str_match(Condition, "^(\\d+)-")[, 2]]
  df_sum[is.na(Bio_Replicate), Bio_Replicate := "00"]
  
  # Condition
  df_sum[, Condition := gsub("^\\d+-", "", Condition)]
  
  # remove blanks
  df_sum <- df_sum[blank == FALSE]
  
  # Average Technical replicates
  df_sum <- df_sum[, list(Area = mean(Area)), by = c("Condition", "Bio_Replicate",  "Peptide", "Precursor.Charge", "Isotope.Label.Type")]
  
  # correct for less spike-in amount
  df_sum[Condition == "WT" & Bio_Replicate %in% c("14", "15", "16") & Isotope.Label.Type == "heavy", Area := Area * 2 ]
  
  
  # log2-transform
  df_sum[, Area := log2(Area)]
  
  # split in heavy and light peptides
  df_h <- df_sum[Isotope.Label.Type == "heavy"]
  df_l <- df_sum[Isotope.Label.Type == "light"]
  
  # derive normalization factors based on heavy peptides
  df_h[, norm_factor := 2^(Area - median(Area)), by = c("Peptide")]
  
  # heavy peptide that is used for normalization
  df_l[Peptide == "EGVVAAAEK", hpept := "EGVVAAAEK"]
  df_l[Peptide == "EGVLYVGSK", hpept := "EGVLYVGSK"]
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
                df_h[, c("Condition", "Bio_Replicate", "Peptide", "Precursor.Charge", "hpept_area", "norm_factor")],
                by.x = c("Condition", "Bio_Replicate", "hpept",   "Precursor.Charge"),
                by.y = c("Condition", "Bio_Replicate", "Peptide", "Precursor.Charge"),
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
  df_l[, Area_mean := mean(Area_norm), by = c("Condition", "Peptide")]
  df_l[, Area_sd   := sd(Area_norm),   by = c("Condition", "Peptide")]
  
  df_l[, aSyn_conc_mean := mean(aSyn_conc), by = c("Condition", "Peptide")]
  df_l[, aSyn_conc_sd   := sd(aSyn_conc),   by = c("Condition", "Peptide")]
  
  
  # indirectly estimate concentrations of "EQVTNVGGAVVTGVTPVAQK" and "KNEEGAPQEGILEDMPVD" in Hu-P3 and Hu-TR
  # set peptide concentrations equal to TVEGAGSIAAATGFVK concentration in P3 an TR
  # use calculate peptide concentration in P3 and TR based on intensity fold changes
  
  # estimated aSyn concentration in P3 based on the TVEGAGSIAAATGFVK peptide
  p3 <- df_l[Peptide == "TVEGAGSIAAATGFVK" & Condition == "P3"]
  
  # assume the same concentration for EQVTNVGGAVVTGVTPVAQK in P3
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "P3", aSyn_conc := p3$aSyn_conc]
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "P3", aSyn_conc_mean := p3$aSyn_conc_mean]
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "P3", aSyn_conc_sd   := p3$aSyn_conc_sd]
  p3 <- df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "P3"]
  
  # based on the intensity fold change calculate EQVTNVGGAVVTGVTPVAQK in HU-P3
  fc <- 2^(df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "HU-P3"]$Area_norm - p3$Area_norm)
  conc <- fc * p3$aSyn_conc
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "HU-P3", aSyn_conc := conc]
  
  # compute mean EQVTNVGGAVVTGVTPVAQK concentration in Hu-P3
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "HU-P3", aSyn_conc_mean := mean(aSyn_conc)]
  
  # use error propagation to caculate sd of the log2 fold change and estimated aSyn concentration
  fc_mean <-  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "HU-P3"]$Area_mean[1] - p3$Area_mean[1]
  fc_sd  <- sqrt(df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "HU-P3"]$Area_sd[1]^2 + p3$Area_sd[1]^2)
  aSyn_mean <-  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "HU-P3"]$aSyn_conc_mean[1]
  aSyn_sd <- 2^(fc_mean + fc_sd)*p3$aSyn_conc_mean[1] - aSyn_mean 
  
  df_l[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition == "HU-P3", aSyn_conc_sd := aSyn_sd]
  
  # repeat the same for KNEEGAPQEGILEDMPVD
  
  # estimated aSyn concentration in TR based on the TVEGAGSIAAATGFVK peptide
  tr <- df_l[Peptide == "TVEGAGSIAAATGFVK" & Condition == "TR"]
  
  # assume the same concentration for EQVTNVGGAVVTGVTPVAQK in P3
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "TR", aSyn_conc := tr$aSyn_conc]
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "TR", aSyn_conc_mean := tr$aSyn_conc_mean]
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "TR", aSyn_conc_sd   := tr$aSyn_conc_sd]
  tr <- df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "TR"]
  
  # based on the intensity fold change calculate EQVTNVGGAVVTGVTPVAQK in HU-P3
  fc <- 2^(df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "HU-TR"]$Area_norm - tr$Area_norm)
  conc <- fc * tr$aSyn_conc
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "HU-TR", aSyn_conc := ..conc]
  
  # compute mean EQVTNVGGAVVTGVTPVAQK concentration in Hu-P3
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "HU-TR", aSyn_conc_mean := mean(aSyn_conc)]
  
  # use error propagation to caculate sd of the log2 fold change and estimated aSyn concentration
  fc_mean <-  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "HU-TR"]$Area_mean[1] - tr$Area_mean[1]
  fc_sd  <- sqrt(df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "HU-TR"]$Area_sd[1]^2 + tr$Area_sd[1]^2)
  aSyn_mean <-  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "HU-TR"]$aSyn_conc_mean[1]
  aSyn_sd <- 2^(fc_mean + fc_sd)*tr$aSyn_conc_mean[1] - aSyn_mean 
  
  df_l[Peptide == "KNEEGAPQEGILEDMPVD" & Condition == "HU-TR", aSyn_conc_sd := aSyn_sd]
  
  fwrite(df_l, "temp\\PeptQuant_Hippoc.txt", sep = "\t")
  fwrite(df_l[, .SD[1], by = c("Peptide", "Condition"), .SDcols = c("Area_mean",
                                                                    "Area_sd",
                                                                    "aSyn_conc_mean",
                                                                    "aSyn_conc_sd")],
          "temp\\PeptQuant_average_Hippoc.txt",
         sep = "\t")
  
  })