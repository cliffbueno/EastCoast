# Function for running correlations with different sets of taxa, subset to a threshold
# Procedure:
# Take an mctoolsr object
# 1. Divide threshold by 100
# 2. Summarize taxonomy
# 3. Get taxa greater than defined threshold of community
# 4. Run correlations
# 5. Plot

meth_corr_by_bgc <- function(env_nona) {
  
  # Get methane
  meth <- env_nona$CH4_ug_m2_h
  
  # Get other BGC variables without methane
  bgc <- env_nona %>%
    dplyr::select(-CH4_ug_m2_h)
  
  # Dataframe to store the correlation results
  methcor <- as.data.frame(matrix(data = NA, nrow = ncol(bgc), ncol = 7)) %>%
    set_names(c("Variable", "r", "PearsonP", "rho", "SpearmanP", "tau", "KendallP"))
  suppressWarnings(for (i in 1:ncol(bgc)) {
    methcor$Variable[i] <- colnames(bgc)[i]
    methcor$r[i] <- round(
      cor.test(bgc[,i], meth, method = "pearson")$estimate[1], 
      digits = 2)
    methcor$PearsonP[i] <-round(
      cor.test(bgc[,i], meth, method = "pearson")$p.value[1], 
      digits = 7)
    methcor$rho[i] <- round(
      cor.test(bgc[,i], meth, method = "spearman")$estimate[1], 
      digits = 2)
    methcor$SpearmanP[i] <-round(
      cor.test(bgc[,i], meth, method = "spearman")$p.value[1], 
      digits = 7)
    methcor$tau[i] <- round(
      cor.test(bgc[,i], meth, method = "kendall")$estimate[1], 
      digits = 2)
    methcor$KendallP[i] <-round(
      cor.test(bgc[,i], meth, method = "kendall")$p.value[1], 
      digits = 7)
  })
  
  # Add adjusted p values to the dataframe
  methcor <- methcor %>% 
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
  
  # Melt the coefficients
  methcor_long_coef <- melt(methcor,
                            measure.vars = c("r", "rho", "tau"),
                            id.vars = c("Variable"))
  # Melt the p values
  methcor_long_pval <- 
    melt(methcor,
         measure.vars = c("PearsonPcut", "SpearmanPcut", "KendallPcut"),
         id.vars = c("Variable"))
  
  # Combine into a long dataframe for plotting
  methcor_long <- data.frame(Variable = methcor_long_coef$Variable, 
                             Test = methcor_long_coef$variable, 
                             Coefficient = methcor_long_coef$value,
                             P = methcor_long_pval$value) %>%
    mutate(Test = dplyr::recode(Test,
                                r = "Pearson",
                                rho = "Spearman",
                                tau = "Kendall"))
  
  # Plot all of the coefficients (order by rho)
  meth_corr_plot <- ggplot(methcor_long, aes(reorder(Variable, abs(Coefficient), max), 
                                             Coefficient, shape = Test, colour = P)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 2, alpha = 0.9) +
    labs(x = NULL,
         y = "Correlation coefficient",
         shape = "Test",
         colour = "Significance") +
    scale_colour_manual(values = c("#F8766D", "#619CFF")) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = c(1,0),
          legend.justification = c(1,0),
          legend.background = element_blank(),
          legend.margin = margin(-0.2,0.2,0.1,0.1, unit="cm"),
          axis.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 6))
  
  meth_corr_plot
}
