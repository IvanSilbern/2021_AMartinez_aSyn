local({
  
  # check linearity of peptide quantification
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(MSstats)
  
  set.seed(100)
  
  lm_eqn <- function(slope, intercept){
    
    eq <- substitute(italic(log[2](Intensity)) == b %.% italic(log[10](fmol))* + italic(a), 
                     list(b  = format(intercept, digits = 3, nsmall = 3),
                          a  = format(slope, digits = 3, nsmall = 3)))
    as.character(as.expression(eq)) 
    
  }
  
  if(!dir.exists("temp")) dir.create("temp")
  if(!dir.exists("plots")) dir.create("plots")
  
  df <- fread("Skyline_results\\20171227_aSyn15N_Dilutions_TransitionResults.tsv", check.names = TRUE)
  df[, Dilution := str_match(Replicate, "^(\\d+)")[, 2]]
  df[, Tech_Replicate := str_match(Replicate, "-(\\d+)$")[, 2]]
  df <- df[Peak.Rank %in% 1:6]
  df <- df[, list(Area = sum(Area)), by = c("Peptide", "Precursor.Charge", "Dilution", "Tech_Replicate")]
  df[, Area := log2(Area)]
  df[, Dilution := factor(Dilution)]
  df[, aSyn_fmol := 690/2^(as.integer(Dilution)-1)]
  df[is.na(Dilution), aSyn_fmol := 0]
  
  df_ms <- df[, c("Peptide", "Area", "aSyn_fmol")]
  df_ms[, Area := 2^Area]
  names(df_ms) <- c("NAME", "INTENSITY", "CONCENTRATION")
  
  ms_mod <- vector("list", length = length(unique(df_ms$NAME)))
  
  for(i in seq_along(ms_mod)){
    
    ms_mod[[i]] <- linear_quantlim(df_ms[NAME == unique(df_ms$NAME)[i]])
  
  }
  
  ms_mod <- rbindlist(ms_mod)
  
  lod_ <- ms_mod[, list(yintercept = (2*LOD[1]*SLOPE[1] + INTERCEPT[1])), by = "NAME"]
    
  # MSstats misidentifies the lod_ for EGVVHGVATVAEK and EGVLYVGSK
  # use one identified based on the graph
  lod_[NAME == "EGVVHGVATVAEK", yintercept := 2^16]
  lod_[NAME == "EGVLYVGSK", yintercept := 2^19]
  lod_[NAME == "EGVVAAAEK", yintercept := 2^18]
  
  # remove EGVVAAEK due to higher variation to prevent overplotting
  # and since the results for EGVVAAAEK are similar to EGVLYVGSK
  lod_ <- lod_[NAME != "EGVVAAAEK"]
  df <- df[Peptide != "EGVVAAAEK"]
  
  lod_ <- lod_[order(NAME)]
  lod_[, linear_range := log2(yintercept)]
  
  # add brakes
  main_breaks <- sort(c(15, 20, 25, 14.4, 16.0, 19.0))
  col_breaks <- c("#00BFC4", "black", "#7CAE00", "#F8766D", "black", "black")
  
  
  g <- ggplot(df, aes(x = log10(aSyn_fmol), y = Area, color = Peptide, group = Peptide))
  g <- g + geom_point()
  g <- g + scale_y_continuous(breaks = main_breaks, minor_breaks = NULL)
  g <- g + ylab(expression(paste("log"[2], " Intensity")))
  g <- g + xlab(expression(paste("log"[10], " alpha-Synuclein, fmol")))
  g <- g + geom_hline(data = lod_, aes(yintercept = log2(yintercept), color = NAME))
  g <- g + theme_bw()
  g <- g + theme(axis.text.y = element_text(color = col_breaks))
  
  pdf("plots\\SupplFig2F_aSyn_peptides_linearity.pdf", width = 10)
  print(g)
  dev.off()
  
  fwrite(lod_, "temp\\Linear_range_aSyn15N.txt")
  fwrite(df, "temp\\Dilutions_aSyn15N.txt", sep = "\t")
  
  })
  