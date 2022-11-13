# Function for running correlations with each environmental variable in a dataframe with no NAs
# Procedure:
# Take an mctoolsr object
# 1. Take environmental dataframe, extract variable of interest, and remove it from df
# 2. Run correlations
# 3. Save df

multiple_correlations <- function(env_nona, var) {
  
  # Get dependent variable of interest
  dv <- env_nona %>%
    dplyr::select(var)
  
  # Get other BGC variables without methane
  bgc <- env_nona %>%
    dplyr::select(var, everything()) %>%
    dplyr::select(-1)
  
  # Dataframe to store the correlation results
  cor.df <- as.data.frame(matrix(data = NA, nrow = ncol(bgc), ncol = 7)) %>%
    set_names(c("Variable", "r", "PearsonP", "rho", "SpearmanP", "tau", "KendallP"))
  suppressWarnings(for (i in 1:ncol(bgc)) {
    cor.df$Variable[i] <- colnames(bgc)[i]
    cor.df$r[i] <- round(
      cor.test(bgc[,i], dv[,1], method = "pearson")$estimate[1], 
      digits = 2)
    cor.df$PearsonP[i] <-round(
      cor.test(bgc[,i], dv[,1], method = "pearson")$p.value[1], 
      digits = 7)
    cor.df$rho[i] <- round(
      cor.test(bgc[,i], dv[,1], method = "spearman")$estimate[1], 
      digits = 2)
    cor.df$SpearmanP[i] <-round(
      cor.test(bgc[,i], dv[,1], method = "spearman")$p.value[1], 
      digits = 7)
    cor.df$tau[i] <- round(
      cor.test(bgc[,i], dv[,1], method = "kendall")$estimate[1], 
      digits = 2)
    cor.df$KendallP[i] <-round(
      cor.test(bgc[,i], dv[,1], method = "kendall")$p.value[1], 
      digits = 7)
  })
  
  # Add adjusted p values to the dataframe
  cor.df <- cor.df %>% 
    mutate(PearsonPfdr = p.adjust(PearsonP, method = "fdr")) %>%
    mutate(SpearmanPfdr = p.adjust(SpearmanP, method = "fdr")) %>%
    mutate(KendallPfdr = p.adjust(KendallP, method = "fdr")) %>%
    mutate(PearsonPcut = factor(ifelse(PearsonPfdr < 0.05, 
                                       "Pfdr < 0.05", "Pfdr > 0.05"),
                                levels = c("Pfdr < 0.05","Pfdr > 0.05"))) %>%
    mutate(SpearmanPcut = factor(ifelse(SpearmanPfdr < 0.05,
                                        "Pfdr < 0.05", "Pfdr > 0.05"),
                                 levels = c("Pfdr < 0.05","Pfdr > 0.05"))) %>%
    mutate(KendallPcut = factor(ifelse(KendallPfdr < 0.05,
                                       "Pfdr < 0.05", "Pfdr > 0.05"),
                                levels = c("Pfdr < 0.05","Pfdr > 0.05")))
  cor.df
}
