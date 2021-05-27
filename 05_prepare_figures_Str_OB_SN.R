local({
  
  # prepare figures for alpha-synuclein peptide concentrations in Striatum, Olfactory bulb and Substantia Nigra
  
  library(data.table)
  library(ggplot2)
  
  if(!dir.exists("plots")) dir.create("plots")
  
  # define colors:
  col_tr   <- rgb(0, 255, 0, maxColorValue = 255)
  #col_hutr <- rgb(0, 0, 255, maxColorValue = 255)
  col_hutr <- "#993300"
  col_hup3 <- rgb(0, 0, 255, maxColorValue = 255)
  col_p3   <- rgb(255, 0, 0, maxColorValue = 255)
  col_wt   <- rgb(0, 0, 0, maxColorValue = 255)
  
  df     <- fread("temp\\PeptQuant_Str_OB_SN.txt")
  df_sum <- fread("temp\\PeptQuant_average_Str_OB_SN.txt")
  
  regions <- c("Striatum", "OlfactoryBulb", "SubstantiaNigra")
  
  # round concentration estimates to 2 significant digits for mean and 1 for sd
  df[, aSyn_conc_mean_round := signif(aSyn_conc_mean, 2)]
  df[, aSyn_conc_sd_round := signif(aSyn_conc_sd, 1)]
  
  df_sum[, aSyn_conc_mean_round := signif(aSyn_conc_mean, 2)]
  df_sum[, aSyn_conc_sd_round   := signif(aSyn_conc_sd, 1)]
  
  # EQVTNVGGAVVTGVTAVAQK
  temp <- df[Peptide == "EQVTNVGGAVVTGVTAVAQK" & Condition %in% c("WT", "TM1", "TM2", "D119")]
  temp[, Condition := factor(Condition, levels = c("WT", "TM1", "TM2", "D119"))]
  
  # add P3 condition (peptide not detected)
  for(i in seq_along(regions)){
    
    temp <- rbind(temp, temp[1])
    temp[nrow(temp), Region := ..regions[i]]
    temp[nrow(temp), Condition := "TP"]
    temp[nrow(temp), Area_norm := 12.5]
    temp[nrow(temp), aSyn_conc := 0]  
    
    }
  
  
  temp[, Condition := factor(Condition, levels = c("WT", "TM1", "TM2", "TP", "D119"))]
  temp[, Region := factor(Region, levels = regions)]
  temp <- temp[order(Region, Condition)]
  
  temp_sum <- df_sum[Peptide == "EQVTNVGGAVVTGVTAVAQK" & Condition %in% c("WT", "TM1", "TM2", "D119")]
  temp_sum[, Region := factor(Region, levels = regions)]
  temp_sum[, mean_sd := unlist(Map(function(x, y){
    
    substitute(mean %+-% sd, list(mean = round(x, 2),  sd = round(y, 2)))  
    
  }, x = aSyn_conc_mean_round, y = aSyn_conc_sd_round))]
  
  # add P3 condition (peptide not detected
  for(i in seq_along(regions)){
    
    temp_sum <- rbind(temp_sum, temp_sum[1])
    temp_sum[nrow(temp_sum), Region := ..regions[i]]
    temp_sum[nrow(temp_sum), Condition := "TP"]
    temp_sum[nrow(temp_sum), Area_mean := 12.5]
    temp_sum[nrow(temp_sum), aSyn_conc_mean_round := 0]
    temp_sum[nrow(temp_sum), aSyn_conc_sd_round := 0]
    temp_sum[nrow(temp_sum), mean_sd := substitute("not~detected")]
    
  }
  
  temp_sum[, Condition := factor(Condition, levels = c("WT", "TM1", "TM2", "TP", "D119"))]
  temp_sum <- temp_sum[order(Region, Condition)]
  
  g <- ggplot(data = temp, aes(x = Condition,
                               y = aSyn_conc,
                               color = Condition))
  g <- g + scale_color_manual(values = c("WT" = col_wt, "TM1" = col_hup3, "TM2" = col_hutr, "TP" = col_p3, "D119" = col_tr ))
  g <- g + facet_wrap(~Region)
  g <- g + geom_boxplot(fill = NA, width = 0.6, middle = temp_sum$aSyn_conc_mean_round,
                        lower  = temp_sum$aSyn_conc_mean_round - temp_sum$aSyn_conc_sd_round,
                        upper = temp_sum$aSyn_conc_mean_round + temp_sum$aSyn_conc_sd_round,
                        ymax = NA, ymin = NA,
                        outlier.color = NA)
  g <- g + geom_point(data = temp[Condition != "TP"], position = position_jitter(width = 0.2, height = 0, seed = 998), size = 3, alpha = 0.8)
  g <- g + scale_x_discrete(labels = c(expression(paste(italic(alpha*"Syn"^{"WT"}))),
                                       expression(paste(italic(alpha*"Syn"^{"tm1"}))),
                                       expression(paste(italic(alpha*"Syn"^{"tm2"}))),
                                       expression(paste(italic(alpha*"Syn"^{"TP"}))),
                                       expression(paste(italic(alpha*"Syn"^{Delta*"119"})))
                                       ))
  g <- g + geom_text(data = temp_sum, aes(label = mean_sd,
                                          y = aSyn_conc_mean_round + 1.5*(aSyn_conc_sd_round + 2)),
                     size = 3, fontface = "bold", nudge_y = 3, parse = TRUE)
    g <- g + ylab(expression(paste({alpha*"Syn, "}, "fmol/"*{mu*"g"}))) + xlab("\n61 EQVTNVGGAVVTGVTAVAQK 80 (Mm/Hu)")
  g <- g + theme_bw()
  g <- g + guides(color = FALSE)
  print(g)
  
  pdf("plots\\SupplFig_EQVTNVGGAVVTGVTAVAQK_Str_OB_SN.pdf", width = 10, height = 5)
  print(g)
  dev.off()
  
  # TVEGAGSIAAATGFVK
  
  temp <- df[Peptide == "TVEGAGSIAAATGFVK" & Condition %in% c("TM1", "TM2", "TP", "D119")]
  temp[, Condition := factor(Condition, levels = c("TM1", "TM2", "TP", "D119"))]
  temp[, Region := factor(Region, levels = regions)]
  temp <- temp[order(Region, Condition)]
  
  temp_sum <- df_sum[Peptide == "TVEGAGSIAAATGFVK" & Condition %in% c("TM1", "TM2", "TP", "D119")]
  temp_sum[, Condition := factor(Condition, levels = c("TM1", "TM2", "TP", "D119"))]
  temp_sum[, Region := factor(Region, levels = regions)]
  temp_sum <- temp_sum[order(Region, Condition)]
  temp_sum[, mean_sd := unlist(Map(function(x, y){
    
    substitute(mean %+-% sd, list(mean = round(x, 2),  sd = round(y, 2)))  
    
  }, x = aSyn_conc_mean_round, y = aSyn_conc_sd_round))]
  
  g <- ggplot(data = temp, aes(x = Condition, y = aSyn_conc, color = Condition))
  g <- g + facet_wrap(~Region)
  g <- g + geom_boxplot(fill = NA, width = 0.6, middle = temp_sum$aSyn_conc_mean_round,
                        lower  = temp_sum$aSyn_conc_mean_round - temp_sum$aSyn_conc_sd_round,
                        upper = temp_sum$aSyn_conc_mean_round + temp_sum$aSyn_conc_sd_round,
                        ymax = NA, ymin = NA,
                        outlier.color = NA)
  g <- g + geom_point(position = position_jitter(width = 0.2, height = 0, seed = 998), size = 3, alpha = 0.8)
  g <- g + scale_color_manual(values = c("TM1" = col_hup3, "TM2" = col_hutr, "TP" = col_p3, "D119" = col_tr))
  g <- g + scale_x_discrete(labels = c(expression(paste(italic(alpha*"Syn"^{"tm1"}))),
                                       expression(paste(italic(alpha*"Syn"^{"tm2"}))),
                                       expression(paste(italic(alpha*"Syn"^{"TP"}))),
                                       expression(paste(italic(alpha*"Syn"^{Delta*"119"})))
                                       ))
  g <- g + geom_text(data = temp_sum, aes(label = mean_sd,
                                          y = aSyn_conc_mean_round + 1.5*aSyn_conc_sd_round),
                     size = 3, fontface = "bold", nudge_y = 3, parse = TRUE)
  g <- g + ylab(expression(paste({alpha*"Syn, "}, "fmol/"*{mu*"g"})))
  g <- g + xlab("\n81 TVEGAGSIAAATGFVK 96 (Hu)")
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(face = "bold"))
  g <- g + guides(color = FALSE)
  print(g)
  
  pdf("plots\\SupplFig_TVEGAGSIAAATGFVK_Str_OB_SN.pdf", width = 10, height = 5)
  print(g)
  dev.off()
  
  # EQVTNVGGAVVTGVTPVAQK
  
  temp <- df[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition %in% c("TM1", "TP")]
  temp[, Condition := factor(Condition, levels = c("TM1", "TP"))]
  temp[, Region := factor(Region, levels = regions)]
  temp <- temp[order(Region, Condition)]
  
  temp_sum <- df_sum[Peptide == "EQVTNVGGAVVTGVTPVAQK" & Condition %in% c("TM1", "TP")]
  temp_sum[, Condition := factor(Condition, levels = c("TM1", "TP"))]
  temp_sum[, Region := factor(Region, levels = regions)]
  temp_sum <- temp_sum[order(Region, Condition)]
  temp_sum[, mean_sd := unlist(Map(function(x, y){
    
    substitute(mean %+-% sd, list(mean = round(x, 1),  sd = round(y, 1)))  
    
  }, x = aSyn_conc_mean_round, y = aSyn_conc_sd_round))]
  
  
  
  g <- ggplot(data = temp, aes(x = Condition, y = aSyn_conc, color = Condition))
  g <- g + facet_wrap(~Region)
  g <- g + geom_boxplot(fill = NA, width = 0.6, middle = temp_sum$aSyn_conc_mean_round,
                        lower  = temp_sum$aSyn_conc_mean_round - temp_sum$aSyn_conc_sd_round,
                        upper = temp_sum$aSyn_conc_mean_round + temp_sum$aSyn_conc_sd_round,
                        ymax = NA, ymin = NA,
                        outlier.color = NA)
  g <- g + geom_point(position = position_jitter(width = 0.2, height = 0, seed = 998), size = 3, alpha = 0.8)
  g <- g + scale_color_manual(values = c("TM1" = col_hup3, "TP" = col_p3))
  g <- g + scale_x_discrete(labels = c(expression(paste(italic(alpha*"Syn"^{"tm1"}))),
                                       expression(paste(italic(alpha*"Syn"^{"TP"})))
  ))
  g <- g + geom_text(data = temp_sum, aes(label = mean_sd,
                                          y = aSyn_conc_mean_round + 1.5*aSyn_conc_sd_round),
                     size = 3, fontface = "bold", nudge_y = 0.2, parse = TRUE)
  g <- g + ylab(expression(paste({alpha*"Syn, "}, "fmol/"*{mu*"g"})))
  g <- g + xlab("\n61 EQVTNVGGAVVTGVTPVAQK 80 (TP)")
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(face = "bold"))
  g <- g + guides(color = FALSE)
  print(g)
  
  pdf("plots\\SupplFig_EQVTNVGGAVVTGVTPVAQK_Str_OB_SN.pdf", width = 10, height = 5)
  print(g)
  dev.off()
  
  # KNEEGAPQEGILEDMPVD
  
  temp <- df[Peptide == "KNEEGAPQEGILEDMPVD" & Condition %in% c("TM2", "D119")]
  temp[, Condition := factor(Condition, levels = c("TM2", "D119"))]
  temp[, Region := factor(Region, levels = regions)]
  temp <- temp[order(Region, Condition)]
  
  temp_sum <- df_sum[Peptide == "KNEEGAPQEGILEDMPVD" & Condition %in% c("TM2", "D119")]
  temp_sum <- temp_sum[order(Condition)]
  temp_sum[, Condition := factor(Condition, levels = c("TM2", "D119"))]
  temp_sum[, Region := factor(Region, levels = regions)]
  temp_sum <- temp_sum[order(Region, Condition)]
  temp_sum[, mean_sd := unlist(Map(function(x, y){
    
    substitute(mean %+-% sd, list(mean = round(x, 1),  sd = round(y, 1)))  
    
  }, x = aSyn_conc_mean_round, y = aSyn_conc_sd_round))]

  
  g <- ggplot(data = temp, aes(x = Condition, y = aSyn_conc, color = Condition))
  g <- g + facet_wrap(~Region)
  g <- g + geom_boxplot(fill = NA, width = 0.6, middle = temp_sum$aSyn_conc_mean_round,
                        lower  = temp_sum$aSyn_conc_mean_round - temp_sum$aSyn_conc_sd_round,
                        upper = temp_sum$aSyn_conc_mean_round + temp_sum$aSyn_conc_sd_round,
                        ymax = NA, ymin = NA,
                        outlier.color = NA) 
  g <- g + geom_point(position = position_jitter(width = 0.2, height = 0, seed = 998), size = 3, alpha = 0.8)
  g <- g + scale_color_manual(values = c("TM2" = col_hutr, "D119" = col_tr))
  g <- g + scale_x_discrete(labels = c(expression(paste(italic(alpha*"Syn"^{"tm2"}))),
                                       expression(paste(italic(alpha*"Syn"^{Delta*"119"})))
  ))
  g <- g + geom_text(data = temp_sum, aes(label = mean_sd,
                                          y = aSyn_conc_mean_round + 1.5*aSyn_conc_sd_round),
                     size = 3, fontface = "bold", nudge_y = 0.2, parse = TRUE)
  g <- g + ylab(expression(paste({alpha*"Syn, "}, "fmol/"*{mu*"g"})))
  g <- g + xlab(expression(paste("\n102 KNEEGAPQEGILEDMPVD 119 (", alpha*"Syn"^{Delta*"119"},")")))
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(face = "bold"))
  g <- g + guides(color = FALSE)
  print(g)
  
  pdf("plots\\SupplFig_KNEEGAPQEGILEDMPVD_Str_OB_SN.pdf", width = 6, height = 5)
  print(g)
  dev.off()
  
  })
  