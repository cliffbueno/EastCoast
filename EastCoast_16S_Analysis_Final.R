# East Coast/SF 16S final data analysis
# by Cliff Bueno de Mesquita, Tringe Lab, JGI, Fall 2022
# Samples from SF Bay, Delaware River, 
# Use just freshwater to oligohaline samples
# Just control and +ASW treatment, or transplants to oligohaline site
# SF: freshwater, oligohaline
# DE: freshwater, oligohaline
# DE: freshwater transplanted to oligohaline
# DE: control, +ASW
# NC: control, +ASW
# SC: control, +ASW
# Assess site, salinity, depth


#### 1. Setup ####
library(plyr) # Data manipulation
library(tidyverse) # Data manipulation
library(mctoolsr) # Microbial analyses
library(RColorBrewer) # Colors
library(vegan) # Multivariate analyses
library(indicspecies) # Indicator species
library(car) # Stats
library(FSA) # SE
library(magrittr) # Set names
library(PMCMRplus) # Stats
library(readxl) # Excel
library(writexl) # Excel
library(plotly) # Interactive plots
library(ggmap) # Maps
library(ggsn) # Maps
library(multcomp) # Tukey HSD and significance letters
library(emmeans) # Tukey HSD and significance letters
library(scales) # View colors
library(cowplot) # Multipanels
library(qvalue) # q values for indicator species
library(reshape2) # melt
library(gridExtra) # graphs
library(grid) # graphs
library(cowplot) # graphs
library(ggpubr) # graphs
library(ggExtra) # graphs
library(ggh4x) # graphs
library(dendextend) # graphs
library(corrplot) # correlation plots
library(pheatmap) # pretty heatmap
library(gplots) # color ramps

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
`%notin%` <- Negate(`%in%`)

# Guild subsetting module from other repository
source("~/Documents/GitHub/SF_microbe_methane/modules/3_OTU_subsetting_modules_v.0.4_strip.r")

# Correlations
source("~/Documents/GitHub/EastCoast/meth_corr_by_taxonomy.R")
source("~/Documents/GitHub/EastCoast/meth_corr_by_bgc.R")
source("~/Documents/GitHub/EastCoast/multiple_correlations.R")

# Plotting
source("~/Documents/GitHub/EastCoast/cliffplot_taxa_bars.R")

# Repository path
setwd("~/Documents/GitHub/EastCoast/")

# Wyatt Hartman's guild color palette
# Note that MeOB don't exist in this dataset, so removed
# But ANME do exist in this dataset, so add
# Extra methanogen guilds added so colors added too
Guild_cols <- read.table("~/Documents/GitHub/SF_microbe_methane/data/colors/Guild_color_palette.txt",
                         sep='\t') %>%
  dplyr::select(Guild, G_index, color) %>%
  set_names(c("Guild", "Index", "color")) %>%
  mutate(Index = rev(Index)) %>%
  add_row(Guild = "ANME", Index = 10, color = "#836FFF") %>%
  add_row(Guild = "CH4_me", Index = 16, color = "#FDC086") %>%
  add_row(Guild = "CH4_mix", Index = 17, color = "#FFFF99") %>%
  filter(Guild != "MeOB") %>%
  arrange(Index)



# Input, filter, and rarefy
input_filt <- readRDS("input_filt_comb_wBGC.rds")
input_filt$map_loaded <- input_filt$map_loaded %>%
  mutate(Estuary = factor(Estuary,
                          levels = c("Waccamaw", "Alligator",
                                     "Delaware", "SF")),
         Salt = "NA",
         Depth = recode_factor(input_filt$map_loaded$Depth,
                                                   " 0-5" = "0-5",
                                                   " 5-15" = "5-15",
                                                   "0.02" = "0-5",
                                                   "0.025" = "0-5",
                                                   "0.1" = "5-15",
                                                   "0.12" = "5-15",
                                                   "0.125" = "5-15"),
         Field = "Lab")

# Add column for if unmanipulated field sample
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "SF") {
    input_filt$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Delaware" & 
      input_filt$map_loaded$Site[i] == "Soil") {
    input_filt$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Alligator" & 
      input_filt$map_loaded$Site[i] == "Soil") {
    input_filt$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Waccamaw" & 
      input_filt$map_loaded$Detail[i] == "Control") {
    input_filt$map_loaded$Field[i] <- "Field"
  }
}

# Add column for salinity class
# SF
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Site[i] == "Sandmound" | input_filt$map_loaded$Site[i] == "West Pond") {
    input_filt$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Site[i] == "Mayberry" | input_filt$map_loaded$Site[i] == "Browns") {
    input_filt$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Site[i] == "Joice" | input_filt$map_loaded$Site[i] == "Rush Ranch") {
    input_filt$map_loaded$Salt[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Site[i] == "Goodyear" | input_filt$map_loaded$Site[i] == "White Slough" |
      input_filt$map_loaded$Site[i] == "Tolay" | input_filt$map_loaded$Site[i] == "China Camp" |
      input_filt$map_loaded$Site[i] == "Muzzi") {
    input_filt$map_loaded$Salt[i] <- "Polyhaline"
  }
}

# SC
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Waccamaw" &
      input_filt$map_loaded$Detail[i] == "Control") {
    input_filt$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Waccamaw" &
      input_filt$map_loaded$Detail[i] == "Freshwater amended") {
    input_filt$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Waccamaw" &
      input_filt$map_loaded$Detail[i] == "Saltwater amended") {
    input_filt$map_loaded$Salt[i] <- "Oligohaline"
  }
}

# NC
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Alligator" &
      input_filt$map_loaded$Detail[i] == "5ppt ASW") {
    input_filt$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Alligator" &
      input_filt$map_loaded$Detail[i] == "SW_noSO4") {
    input_filt$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Alligator" &
      input_filt$map_loaded$Detail[i] == "DI_ctrl") {
    input_filt$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Alligator" &
      input_filt$map_loaded$Detail[i] == "SO4 amended") {
    input_filt$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Estuary[i] == "Alligator" &
      input_filt$map_loaded$Detail[i] == "Field Reference") {
    input_filt$map_loaded$Salt[i] <- "Freshwater"
  }
}

# Delaware River
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Info[i] == "TFM1_source" | 
      input_filt$map_loaded$Info[i] == "TFM2_source" |
      input_filt$map_loaded$Info[i] == "TFM1@TFM2" |
      input_filt$map_loaded$Info[i] == "TFM1@TFM2-40cm" |
      input_filt$map_loaded$Info[i] == "Freshwater" |
      input_filt$map_loaded$Info[i] == "Freshwater 4wk" |
      input_filt$map_loaded$Info[i] == "Freshwater 7wk" |
      input_filt$map_loaded$Info[i] == "Freshwater 12wk" |
      input_filt$map_loaded$Info[i] == "TFM2_source") {
    input_filt$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Info[i] == "OligoHal_source" |
      input_filt$map_loaded$Info[i] == "TFM1@OligoHal" |
      input_filt$map_loaded$Info[i] == "TFM1@OligoHal-40cm" |
      input_filt$map_loaded$Info[i] == "5 ppt ASW 4wk" |
      input_filt$map_loaded$Info[i] == "5 ppt ASW 7wk" |
      input_filt$map_loaded$Info[i] == "5 ppt ASW 12wk") {
    input_filt$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(input_filt$map_loaded)) {
  if (input_filt$map_loaded$Info[i] == "MesoHal_source" |
      input_filt$map_loaded$Info[i] == "TFM1@MesoHal" |
      input_filt$map_loaded$Info[i] == "TFM1@MesoHal-40cm") {
    input_filt$map_loaded$Salt[i] <- "Mesohaline"
  }
}

input_filt$map_loaded <- input_filt$map_loaded %>%
  unite("EstSalt", c(Estuary, Salt), sep = "_", remove = F) %>%
  mutate(EstSalt = factor(EstSalt,
                          levels = c("SF_Freshwater", "Alligator_Freshwater",
                                     "Delaware_Freshwater", "Waccamaw_Freshwater",
                                     "SF_Oligohaline", "Alligator_Oligohaline",
                                     "Delaware_Oligohaline", "Waccamaw_Oligohaline",
                                     "SF_Mesohaline", "Delaware_Mesohaline",
                                     "SF_Polyhaline")),
         Salt = factor(Salt,
                       levels = c("Freshwater", "Oligohaline", 
                                  "Mesohaline", "Polyhaline")))

# Filter freshwater and oligohaline samples
frol <- filter_data(input_filt,
                    filter_cat = "Salt",
                    keep_vals = c("Freshwater", "Oligohaline"))

# Filter just Control and +ASW
frol <- filter_data(frol,
                    filter_cat = "Detail",
                    filter_vals = c("SW_noSO4", "SO4 amended"))

# Filter just final time points
frol <- filter_data(frol,
                    filter_cat = "Detail",
                    filter_vals = c("Freshwater", "Freshwater 4wk", "Freshwater 7wk",
                                    "5 ppt ASW 4wk", "5 ppt ASW 7wk"))

# Filter to keep just the more flooded "organs" (-40 cm elevation)
# Makes sense because Alligator data are from flooded samples
# Waccamaw samples get 10-30 cm of flooding
frol <- filter_data(frol,
                    filter_cat = "Detail",
                    filter_vals = c("TFM1@TFM2", "TFM1@OligoHal")) # n = 157

# Filter out unwanted NC samples
frol <- filter_data(frol,
                    filter_cat = "sampleID",
                    filter_vals = c("TL_nw_d1_DI_ctrl_AF1", "TL_nw_d1_DI_ctrl_AF3", 
                                    "TL_nw_d1_DI_ctrl_AF4", "TL_nw_d1_ASW_noS_BF3",
                                    "TL_nw_d1_ASW_noS_BF4", "TL_nw_d1_ASW_noS_BF5"))

# Filter out DE field samples that were separate and don't have CH4
frol <- filter_data(frol,
                    filter_cat = "Detail",
                    filter_vals = c("TFM1_source", "TFM2_source", "OligoHal_source"))

# Filter out NC field samples that don't have CH4
frol <- filter_data(frol,
                    filter_cat = "Detail",
                    filter_vals = c("Field Reference"))

# Rarefy at 31264
sort(colSums(frol$data_loaded))
mean(colSums(frol$data_loaded))
se(colSums(frol$data_loaded))
# Original depth 131844.9 ± 3373.033
# Drop Sandmound_TuleB_D1 (1296 reads)
set.seed(530)
frol <- single_rarefy(frol, 31264) # Now n = 133
sort(colSums(frol$data_loaded))

# OTU Richness
frol$map_loaded$rich <- specnumber(frol$data_loaded, 
                                              MARGIN = 2)

# Shannon diversity
frol$map_loaded$shannon <- diversity(frol$data_loaded, 
                                                index = "shannon", 
                                                MARGIN = 2)

# Save
# saveRDS(frol, "frol.rds")



#### 2. Combined ####
frol <- readRDS("frol.rds")
frol$map_loaded$EstSalt <- factor(frol$map_loaded$EstSalt,
                                  levels = c("Alligator_Oligohaline", "Alligator_Freshwater",
                                            "Delaware_Oligohaline", "Delaware_Freshwater", 
                                            "SF_Oligohaline", "SF_Freshwater",
                                            "Waccamaw_Oligohaline", "Waccamaw_Freshwater"))
table(frol$map_loaded$Estuary)

#### _Alpha ####
leveneTest(frol$map_loaded$rich ~ frol$map_loaded$Estuary) # Homogeneous
leveneTest(frol$map_loaded$rich ~ frol$map_loaded$Salt) # Homogeneous
leveneTest(frol$map_loaded$rich ~ frol$map_loaded$Depth) # Homogeneous
m <- aov(rich ~ Estuary + Salt + Depth, data = frol$map_loaded)
Anova(m, type = "II") # All sig.
m <- aov(rich ~ EstSalt, data = frol$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "EstSalt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(frol$map_loaded$rich)+(max(frol$map_loaded$rich)-min(frol$map_loaded$rich))/20)

leveneTest(frol$map_loaded$rich ~ frol$map_loaded$Estuary) # Almost homogeneous
leveneTest(frol$map_loaded$shannon ~ frol$map_loaded$Salt) # Homogeneous
leveneTest(frol$map_loaded$rich ~ frol$map_loaded$Depth) # Homogeneous
m1 <- aov(shannon ~ Estuary + Salt + Depth, data = frol$map_loaded)
Anova(m1, type = "II") # All sig
m1 <- aov(shannon ~ EstSalt, data = frol$map_loaded)
shapiro.test(m1$residuals) # Not normal
summary(m1)
t1 <- emmeans(object = m1, specs = "EstSalt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(frol$map_loaded$shannon)+(max(frol$map_loaded$shannon)-min(frol$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- frol$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
png("FinalFigs/Figure1.png", width = 7, height = 4.5, units = "in", res = 300)
ggplot(alpha_long, aes(EstSalt, value)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Salt)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, 
              aes(fill = Depth, shape = Estuary, colour = Salt)) +
  geom_text(data = label_df, aes(EstSalt, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site",
       colour = "Salinity",
       fill = "Depth (cm)",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(shape = guide_legend(order = 1,),
         colour = guide_legend(order = 2),
         fill = guide_legend(override.aes = list(shape = c(16, 1)), order = 3)) +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0.25, 0, 0, -0.1, unit = "cm"),
        legend.key.size = unit(0.4, "cm"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 25, "pt"))
dev.off()



#### _Beta  ####
# Bray-Curtis dissimilarity
bc <- calc_dm(frol$data_loaded)

# Stats
set.seed(100)
adonis2(bc ~ frol$map_loaded$Estuary + frol$map_loaded$Salt + frol$map_loaded$Depth) # Salt and Depth sig
anova(betadisper(bc, frol$map_loaded$Estuary)) # Dispersion not homogeneous
anova(betadisper(bc, frol$map_loaded$Salt)) # Dispersion not homogeneous
anova(betadisper(bc, frol$map_loaded$Depth)) # Dispersion homogeneous

# Get variables
env <- frol$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, Salinity_ppt_all)
env_nona <- na.omit(env) # n = 110

# Envfit
pcoa <- cmdscale(bc, k = nrow(frol$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.8
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("CH4", "Salinity"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
frol$map_loaded$Axis01 <- scores(pcoa)[,1]
frol$map_loaded$Axis02 <- scores(pcoa)[,2]

# Get legend
micro.hulls <- ddply(frol$map_loaded, c("Salt"), find_hull)
g_forleg <- ggplot(frol$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt), fill = NA,
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(fill = Depth, colour = Salt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salinity",
       fill = "Depth (cm)",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(shape = guide_legend(order = 1,),
         colour = guide_legend(order = 2),
         fill = guide_legend(override.aes = list(shape = c(16, 1)), order = 3)) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0.25, 0, 0, -0.1, unit = "cm"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
leg <- get_legend(g_forleg)

# Plot
micro.hulls <- ddply(frol$map_loaded, c("EstSalt"), find_hull)
g <- ggplot(frol$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = EstSalt), fill = NA,
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(fill = Depth, colour = Salt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salt",
       fill = "Depth (cm)",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red", "blue", "red", "blue", "red",
                                 "blue", "red", "blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(shape = guide_legend(order = 1,),
         colour = guide_legend(order = 2),
         fill = guide_legend(override.aes = list(shape = c(16, 1)), order = 3)) +
  theme_bw() +  
  ggtitle("(a) Combined") +
  theme(legend.position = "none",
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0.25, 0, 0, -0.1, unit = "cm"),
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt")) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 2, color = "black")
g

# Ordistep
comm_nona <- as.data.frame(t(frol$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova #

# Make other PCoAs individually
sf <- filter_data(frol,
                  filter_cat = "Estuary",
                  keep_vals = "SF")
bc_sf <- calc_dm(sf$data_loaded)
env_sf <- sf$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, 
                Salinity_ppt_all, SO4_mgL, DOC_mgL, Fe_mgL, Mn_mgL, Cu_mgL, Zn_mgL,
                sed_pH, sed_NH4_mgL, sed_NO3_mgL, sed_PO4_mgL, sed_Cl_mgL, sed_SO4_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_Bulk_dens,
                sed_Fe_mgL, sed_Mn_mgL, sed_Cu_mgL, sed_Zn_mgL)
env_nona_sf <- na.omit(env_sf) # n = 30 (just D2)
pcoa_sf <- cmdscale(bc_sf, k = nrow(sf$map_loaded) - 1, eig = T)
set.seed(100)
ef_sf <- envfit(pcoa_sf, env_sf, permutations = 999, na.rm = TRUE)
ef_sf
ordiplot(pcoa_sf)
plot(ef_sf, p.max = 0.05, cex = 0.5)
manual_factor_sf <- 0.3
vec.df_sf <- as.data.frame(ef_sf$vectors$arrows*sqrt(ef_sf$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_sf,
         Dim2 = Dim2 * manual_factor_sf) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_sf$vectors$pvals < 0.05) %>%
  filter(variables != "sed_Cl_mgL") %>%
  mutate(shortnames = c("CH4", "CO2", "Salinity", "DOC", "Fe", "Mn",
                        "Zn", "sed_pH", "sed_NH4", "sed_NO3", "sed_PO4", "C", 
                        "N", "C:N", "BD", "sed_Fe", "sed_Mn", "sed_Cu", "sed_Zn"))
pcoaA1 <- round((eigenvals(pcoa_sf)/sum(eigenvals(pcoa_sf)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_sf)/sum(eigenvals(pcoa_sf)))[2]*100, digits = 1)
sf$map_loaded$Axis01 <- scores(pcoa_sf)[,1]
sf$map_loaded$Axis02 <- scores(pcoa_sf)[,2]
micro.hulls <- ddply(sf$map_loaded, c("Salt"), find_hull)
g2 <- ggplot(sf$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt), fill = NA,
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = Salt, fill = Depth), shape = 24) +
  geom_segment(data = vec.df_sf,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = subset(vec.df_sf, shortnames != "CO2" & shortnames != "CH4"),
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 2, color = "black") +
  geom_text(data = subset(vec.df_sf, shortnames == "CO2"),
            aes(x = Dim1, y = Dim2 - 0.03, label = shortnames),
            size = 2, color = "black") +
  geom_text(data = subset(vec.df_sf, shortnames == "CH4"),
            aes(x = Dim1 - 0.04, y = Dim2, label = shortnames),
            size = 2, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  theme_bw() +
  ggtitle("(b) SF") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g2

de <- filter_data(frol,
                  filter_cat = "Estuary",
                  keep_vals = "Delaware")
detra <- filter_data(de,
                     filter_cat = "Site",
                     keep_vals = "Soil mesocosm")
bc_detra <- calc_dm(detra$data_loaded)
env_detra <- detra$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, N2O_ug_m2_h, 
                Salinity_ppt_all, 
                NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL,
                Fe_mgL, Porosity, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d)
env_nona_detra <- na.omit(env_detra) # n = 8
pcoa_detra <- cmdscale(bc_detra, k = nrow(detra$map_loaded) - 1, eig = T)
set.seed(100)
ef_detra <- envfit(pcoa_detra, env_detra, permutations = 999, na.rm = TRUE)
ef_detra
ordiplot(pcoa_detra)
plot(ef_detra, p.max = 0.05, cex = 0.5)
manual_factor_detra <- 0.64
vec.df_detra <- as.data.frame(ef_detra$vectors$arrows*sqrt(ef_detra$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_detra,
         Dim2 = Dim2 * manual_factor_detra) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_detra$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "SR"))
pcoaA1 <- round((eigenvals(pcoa_detra)/sum(eigenvals(pcoa_detra)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_detra)/sum(eigenvals(pcoa_detra)))[2]*100, digits = 1)
detra$map_loaded$Axis01 <- scores(pcoa_detra)[,1]
detra$map_loaded$Axis02 <- scores(pcoa_detra)[,2]
micro.hulls <- ddply(detra$map_loaded, c("Salt"), find_hull)
g3 <- ggplot(detra$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt), fill = NA,
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = Salt, fill = Depth), shape = 23,
             show.legend = F) +
  geom_segment(data = vec.df_detra,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_detra,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 2, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  theme_bw() +
  ggtitle("(c) Delaware (field exp.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g3

sc <- filter_data(frol,
                  filter_cat = "Estuary",
                  keep_vals = "Waccamaw")
bc_sc <- calc_dm(sc$data_loaded)
env_sc <- sc$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, CH4_pw_air_ppmv,
                CH4_pot_umol_gdw_h, CO2_pot_umol_gdw_h,
                N2_umol_m2_h,	SOD_umol_m2_h, NO3_umol_m2_h,	NH4_umol_m2_h, SRP_umol_m2_h, DON_umol_m2_h,
                Conductivity_uS_cm, Salinity_ppt_all, DIC_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_per_org, sed_per_inorg, pH)
env_nona_sc <- na.omit(env_sc) # n = 15
pcoa_sc <- cmdscale(bc_sc, k = nrow(sc$map_loaded) - 1, eig = T)
set.seed(100)
ef_sc <- envfit(pcoa_sc, env_sc, permutations = 999, na.rm = TRUE)
ef_sc
ordiplot(pcoa_sc)
plot(ef_sc, p.max = 0.05, cex = 0.5)
manual_factor_sc <- 0.2
vec.df_sc <- as.data.frame(ef_sc$vectors$arrows*sqrt(ef_sc$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_sc,
         Dim2 = Dim2 * manual_factor_sc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_sc$vectors$pvals < 0.05) %>%
  filter(variables != "Conductivity_uS_cm") %>%
  mutate(shortnames = c("CO2", "N2", "SOD", "Salinity", "N", "C:N", "pH"))
pcoaA1 <- round((eigenvals(pcoa_sc)/sum(eigenvals(pcoa_sc)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_sc)/sum(eigenvals(pcoa_sc)))[2]*100, digits = 1)
sc$map_loaded$Axis01 <- scores(pcoa_sc)[,1]
sc$map_loaded$Axis02 <- scores(pcoa_sc)[,2]
micro.hulls <- ddply(sc$map_loaded, c("Salt"), find_hull)
g4 <- ggplot(sc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt), fill = NA,
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = Salt, fill = Depth), shape = 21,
             show.legend = F) +
  geom_segment(data = vec.df_sc,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_sc,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 2, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  theme_bw() +
  ggtitle("(d) Waccamaw (field exp.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g4

# Note: in initial analysis, thought there was error, but looks like maybe there isn't
# Just not clear separation, no hard evidence of any sample label error
deinc <- filter_data(de,
                     filter_cat = "Site",
                     keep_vals = "Soil incubation")
bc_deinc <- calc_dm(deinc$data_loaded)
env_deinc <- deinc$map_loaded %>%
  dplyr::select(Salinity_ppt_all, 
                NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL,
                Fe_mgL, Porosity, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d)
env_nona_deinc <- na.omit(env_deinc) # n = 8
pcoa_deinc <- cmdscale(bc_deinc, k = nrow(deinc$map_loaded) - 1, eig = T)
set.seed(100)
ef_deinc <- envfit(pcoa_deinc, env_deinc, permutations = 999, na.rm = TRUE)
ef_deinc
ordiplot(pcoa_deinc)
plot(ef_deinc, p.max = 0.05, cex = 0.5)
manual_factor_deinc <- 0.5
vec.df_deinc <- as.data.frame(ef_deinc$vectors$arrows*sqrt(ef_deinc$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_deinc,
         Dim2 = Dim2 * manual_factor_deinc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_deinc$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "SR"))
pcoaA1 <- round((eigenvals(pcoa_deinc)/sum(eigenvals(pcoa_deinc)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_deinc)/sum(eigenvals(pcoa_deinc)))[2]*100, digits = 1)
deinc$map_loaded$Axis01 <- scores(pcoa_deinc)[,1]
deinc$map_loaded$Axis02 <- scores(pcoa_deinc)[,2]
micro.hulls <- ddply(deinc$map_loaded, c("Salt"), find_hull)
g5 <- ggplot(deinc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt), fill = NA,
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = Salt, fill = Depth), shape = 23,
             show.legend = F) +
#  geom_segment(data = vec.df_deinc,
#               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
#               arrow = arrow(length = unit(0.5, "cm")),
#               colour = "gray", alpha = 0.5,
#               inherit.aes = FALSE) + 
#  geom_text(data = vec.df_deinc,
#            aes(x = Dim1, y = Dim2, label = shortnames),
#            size = 2, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  theme_bw() +
  ggtitle("(e) Delaware (lab exp.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g5

nc <- filter_data(frol,
                  filter_cat = "Estuary",
                  keep_vals = "Alligator")
bc_nc <- calc_dm(nc$data_loaded)
env_nc <- nc$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, N2O_ug_m2_h, Salinity_ppt_all, 
                TOC_mgL, TN_mgL, NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL, Br_mgL, NO3_mgL,
                DIN_mgL, DON_mgL, pH)
env_nona_nc <- na.omit(env_nc) # n = 8
pcoa_nc <- cmdscale(bc_nc, k = nrow(nc$map_loaded) - 1, eig = T)
set.seed(100)
ef_nc <- envfit(pcoa_nc, env_nc, permutations = 999, na.rm = TRUE)
ef_nc
ordiplot(pcoa_nc)
plot(ef_nc, p.max = 0.05, cex = 0.5)
manual_factor_nc <- 0.33
vec.df_nc <- as.data.frame(ef_nc$vectors$arrows*sqrt(ef_nc$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_nc,
         Dim2 = Dim2 * manual_factor_nc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_nc$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("CH4", "CO2", "N2O", "Salinity", "NH4", "Br"))
pcoaA1 <- round((eigenvals(pcoa_nc)/sum(eigenvals(pcoa_nc)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_nc)/sum(eigenvals(pcoa_nc)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(pcoa_nc)[,1]
nc$map_loaded$Axis02 <- scores(pcoa_nc)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Salt"), find_hull)
g6 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt), fill = NA,
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = Salt, fill = Depth), shape = 22,
             show.legend = F) +
  geom_segment(data = vec.df_nc,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_nc,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 2, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  theme_bw() +
  ggtitle("(f) Alligator (lab exp.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g6

pcoa_plots <- plot_grid(g, g2, g3, g4, g5, g6, ncol = 2)
fig2 <- plot_grid(pcoa_plots, leg, rel_widths = c(0.85, 0.15))
fig2
png("FinalFigs/Figure2.png", width = 8, height = 6, units = "in", res = 300)
fig2
dev.off()


#### _Taxa ####
cliffplot_taxa_bars(input = frol, level = 1, variable = "EstSalt")
cliffplot_taxa_bars(input = frol, level = 2, variable = "EstSalt")
cliffplot_taxa_bars(input = frol, level = 3, variable = "EstSalt")
cliffplot_taxa_bars(input = frol, level = 4, variable = "EstSalt")
cliffplot_taxa_bars(input = frol, level = 5, variable = "EstSalt")
cliffplot_taxa_bars(input = frol, level = 6, variable = "EstSalt")
cliffplot_taxa_bars(input = frol, level = 9, variable = "EstSalt")

# Phyla, all samples
tax_sum_phyla <- summarize_taxonomy(input = frol, level = 2, report_higher_tax = F)
frol$map_loaded$sampleID <- rownames(frol$map_loaded)
barsP <- plot_taxa_bars(tax_sum_phyla,
                       frol$map_loaded,
                       "sampleID",
                       num_taxa = 12,
                       data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., frol$map_loaded, by = c("group_by" = "sampleID"))
facet_names <- c("Waccamaw" = "Waccamaw",
                 "Alligator" = "Alligator",
                 "Delaware" = "Delaware",
                 "SF" = "SF",
                 "Freshwater" = "Fresh",
                 "Oligohaline" = "Oligo")
phy <- ggplot(barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_nested(~ Estuary + Salt, space = "free", scales = "free_x", 
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 25, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 25, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
phy_leg <- get_legend(phy)
phy <- phy + theme(legend.position = "none")

# Guilds all samples
tax_sum_guilds <- summarize_taxonomy(input = frol, level = 9, report_higher_tax = F)
frol$map_loaded$sampleID <- rownames(frol$map_loaded)
barsG <- plot_taxa_bars(tax_sum_guilds,
                       frol$map_loaded,
                       "sampleID",
                       num_taxa = 20,
                       data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., frol$map_loaded, by = c("group_by" = "sampleID"))
tallest_bar <- barsG %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
gui <- ggplot(barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Estuary + Salt, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 2, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_blank(),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 10, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 10, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
gui_leg <- get_legend(gui)
gui <- gui + theme(legend.position = "none")

pg <- plot_grid(phy, gui, ncol = 1)
ll <- plot_grid(phy_leg, gui_leg, ncol = 1)

plot_grid(pg, ll, rel_widths = c(0.75, 0.25))

png("Figure3_forPPT.png", width = 5.5, height = 6, units = "in", res = 300)
pg
dev.off()

# Stats
# Test effect of salinity class on guilds and top phyla
tax_sum_phyla_sc <- summarize_taxonomy(input = sc, level = 2, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsP$taxon)
res_phy_sc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla_sc,
                                          metadata_map = sc$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "Waccamaw",
         Level = "Phylum",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))
tax_sum_guilds_sc <- summarize_taxonomy(input = sc, level = 9, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsG$taxon)
res_gui_sc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_guilds_sc,
                                          metadata_map = sc$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "Waccamaw",
         Level = "Guild",
         Taxon = rownames(.)) %>%
  arrange(match(Taxon, Guild_cols$Guild)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_phyla_nc <- summarize_taxonomy(input = nc, level = 2, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsP$taxon)
res_phy_nc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla_nc,
                                          metadata_map = nc$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "Alligator",
         Level = "Phylum",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_guilds_nc <- summarize_taxonomy(input = nc, level = 9, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsG$taxon)
tax_sum_guilds_nc[15,] <- rep(0, ncol(tax_sum_guilds_nc))
tax_sum_guilds_nc[16,] <- rep(0, ncol(tax_sum_guilds_nc))
tax_sum_guilds_nc[17,] <- rep(0, ncol(tax_sum_guilds_nc))
rownames(tax_sum_guilds_nc)[15:17] <- c("Anamx", "ANME", "CH4_ac")
res_gui_nc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_guilds_nc,
                                          metadata_map = nc$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "Waccamaw",
         Level = "Guild",
         Taxon = rownames(.)) %>%
  arrange(match(Taxon, Guild_cols$Guild)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_phyla_de <- summarize_taxonomy(input = de, level = 2, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsP$taxon)
res_phy_de <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla_de,
                                          metadata_map = de$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "Delaware",
         Level = "Phylum",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_guilds_de <- summarize_taxonomy(input = de, level = 9, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsG$taxon)
res_gui_de <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_guilds_de,
                                          metadata_map = de$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "Delaware",
         Level = "Guild",
         Taxon = rownames(.)) %>%
  arrange(match(Taxon, Guild_cols$Guild)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_phyla_sf <- summarize_taxonomy(input = sf, level = 2, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsP$taxon)
res_phy_sf <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla_sf,
                                          metadata_map = sf$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "SF",
         Level = "Phylum",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_guilds_sf <- summarize_taxonomy(input = sf, level = 9, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsG$taxon)
res_gui_sf <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_guilds_sf,
                                          metadata_map = sf$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "SF",
         Level = "Guild",
         Taxon = rownames(.)) %>%
  arrange(match(Taxon, Guild_cols$Guild)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

taxa_results <- rbind(res_phy_sc, res_phy_nc, res_phy_de, res_phy_sf,
                      res_gui_sc, res_gui_nc, res_gui_de, res_gui_sf)

res_phy <- res_phy_sc %>%
  mutate(WA = SaltEffect) %>%
  dplyr::select(Level, Taxon, WA) %>%
  mutate(AL = res_phy_nc$SaltEffect,
         DE = res_phy_de$SaltEffect,
         SF = res_phy_sf$SaltEffect)

res_gui <- res_gui_sc %>%
  mutate(WA = SaltEffect) %>%
  dplyr::select(Level, Taxon, WA) %>%
  mutate(AL = res_gui_nc$SaltEffect,
         DE = res_gui_de$SaltEffect,
         SF = res_gui_sf$SaltEffect)

res <- rbind(res_phy, res_gui)
write_xlsx(res, "taxa_results.xlsx", format_headers = F)
brewer.pal(12, "Paired")


#### _Indicators ####
# MULTIPATT (list ASVs associated with each group)
set.seed(1202)
frol_mp <- multipatt(t(frol$data_loaded), 
                     frol$map_loaded$Salt, 
                     func = "r.g", 
                     control = how(nperm=999),
                     max.order = 1)
frol_mp_results <- frol_mp$sign %>%
  mutate(q.value = qvalue(frol_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(frol_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.45) %>%
  left_join(., frol$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(frol_mp_results)) {
  if (frol_mp_results$s.Freshwater[i] == 1) {
    frol_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(frol_mp_results)) {
  if (frol_mp_results$s.Oligohaline[i] == 1) {
    frol_mp_results$Group[i] <- "Oligohaline"
  }
}
table(frol_mp_results$Group)
frol_asv <- summarize_taxonomy(frol, level = 8, report_higher_tax = F)
frol_asv_all <- data.frame("RelAbundance" = round(rowMeans(frol_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
frol_mp_corrs <- as.data.frame(frol_mp$str) %>%
  dplyr::select(1:length(levels(frol$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% frol_mp_results$ASV) %>%
  set_names(c("Freshwater", "Oligohaline", "ASV"))
# Add corrs and taxonomy
frol_mp_results <- frol_mp_results %>%
  left_join(., frol_asv_all, by = "ASV") %>%
  left_join(., frol_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(frol_mp_corrs)[1:2], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(frol_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
frol.hm.melted <- frol_mp_results %>%
  dplyr::select(taxon, names(frol_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
frol.hm <- ggplot(data = frol.hm.melted, 
                  aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(frol.hm.melted$taxon), labels = unique(frol.hm.melted$taxon),
                   limits = rev(levels(frol.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
frol.l <- get_legend(frol.hm)
frol.hm.clean <- frol.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)),
                                                                   size = 5),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
frol.bp.y <- ggplot(data = frol_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(frol_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")

top <- plot_grid(frol.hm.clean, frol.bp.y, ncol = 2, rel_widths = c(0.7, 0.3))
png("FinalFigs/Figure4.png", width = 5, height = 5, units = "in", res = 300)
plot_grid(top, frol.l, nrow = 2, rel_heights = c(0.87, 0.13))
dev.off()



#### _CH4 ####
# CH4 by salt and experiment
# Need to break into 5, add addition info
frol$map_loaded$Exp <- "NA"
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Delaware" & 
      frol$map_loaded$Site[i] == "Soil mesocosm") {
    frol$map_loaded$Exp[i] <- "Field experiment"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Delaware" & 
      frol$map_loaded$Site[i] == "Soil incubation") {
    frol$map_loaded$Exp[i] <- "Lab experiment"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "SF") {
    frol$map_loaded$Exp[i] <- "Observational"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Waccamaw") {
    frol$map_loaded$Exp[i] <- "Field experiment"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Alligator") {
    frol$map_loaded$Exp[i] <- "Lab experiment"
  }
}
frol$map_loaded$Exp <- factor(frol$map_loaded$Exp,
                              levels = c("Lab experiment", "Field experiment",
                                         "Observational"))
frol$map_loaded$Estuary2 <- factor(frol$map_loaded$Estuary,
                                   levels = c("Alligator", 
                                              "Delaware",
                                              "Waccamaw",
                                              "SF"))
# Stats
leveneTest(CH4_ug_m2_h ~ Salt, data = sf$map_loaded)
shapiro.test(sf$map_loaded$CH4_ug_m2_h)
t.test(CH4_ug_m2_h ~ Salt, data = sf$map_loaded) # Sig.
wilcox.test(CH4_ug_m2_h ~ Salt, data = sf$map_loaded) # Sig.

leveneTest(CH4_ug_m2_h ~ Salt, data = nc$map_loaded)
shapiro.test(nc$map_loaded$CH4_ug_m2_h)
t.test(CH4_ug_m2_h ~ Salt, data = nc$map_loaded) # Marg.
wilcox.test(CH4_ug_m2_h ~ Salt, data = nc$map_loaded) # Marg.

leveneTest(CH4_ug_m2_h ~ Salt, data = sc$map_loaded)
shapiro.test(sc$map_loaded$CH4_ug_m2_h)
t.test(CH4_ug_m2_h ~ Salt, data = sc$map_loaded) # Sig.
wilcox.test(CH4_ug_m2_h ~ Salt, data = sc$map_loaded) # Sig.

leveneTest(CH4_ug_m2_h ~ Salt, data = detra$map_loaded)
shapiro.test(detra$map_loaded$CH4_ug_m2_h)
t.test(CH4_ug_m2_h ~ Salt, data = detra$map_loaded) # NSD
wilcox.test(CH4_ug_m2_h ~ Salt, data = detra$map_loaded) # NSD

leveneTest(CH4_ug_m2_h ~ Salt, data = deinc$map_loaded)
shapiro.test(deinc$map_loaded$CH4_ug_m2_h)
t.test(CH4_ug_m2_h ~ Salt, data = deinc$map_loaded)
wilcox.test(CH4_ug_m2_h ~ Salt, data = deinc$map_loaded)

label_df2 <- data.frame("Exp" = c("Lab experiment", "Lab experiment", 
                                  "Lab experiment", "Lab experiment",
                                  "Field experiment", "Field experiment",
                                  "Field experiment", "Field experiment",
                                  "Observational", "Observational"),
                        "Estuary2" = c("Alligator", "Alligator",
                                       "Delaware", "Delaware",
                                       "Delaware", "Delaware",
                                       "Waccamaw", "Waccamaw",
                                       "SF", "SF"),
                        "y" = c(max(nc$map_loaded$CH4_ug_m2_h) + max(nc$map_loaded$CH4_ug_m2_h)/10,
                                max(nc$map_loaded$CH4_ug_m2_h) + max(nc$map_loaded$CH4_ug_m2_h)/10,
                              max(deinc$map_loaded$CH4_ug_m2_h) + max(deinc$map_loaded$CH4_ug_m2_h)/10,
                              max(deinc$map_loaded$CH4_ug_m2_h) + max(deinc$map_loaded$CH4_ug_m2_h)/10,
                              max(detra$map_loaded$CH4_ug_m2_h) + max(detra$map_loaded$CH4_ug_m2_h)/10,
                              max(detra$map_loaded$CH4_ug_m2_h) + max(detra$map_loaded$CH4_ug_m2_h)/10,
                              max(sc$map_loaded$CH4_ug_m2_h) + max(sc$map_loaded$CH4_ug_m2_h)/10,
                              max(sc$map_loaded$CH4_ug_m2_h) + max(sc$map_loaded$CH4_ug_m2_h)/10,
                              max(sf$map_loaded$CH4_ug_m2_h) + max(sf$map_loaded$CH4_ug_m2_h)/10,
                              max(sf$map_loaded$CH4_ug_m2_h) + max(sf$map_loaded$CH4_ug_m2_h)/10),
                        "x" = c("Freshwater", "Oligohaline", "Freshwater", "Oligohaline",
                                "Freshwater", "Oligohaline", "Freshwater", "Oligohaline",
                                "Freshwater", "Oligohaline"),
                        "label" = c("a", "b", "", "", "", "", "a", "b", "a", "b"))
label_df2$Exp <- factor(label_df2$Exp,
                              levels = c("Lab experiment", "Field experiment",
                                         "Observational"))
label_df2$Estuary2 <- factor(label_df2$Estuary2,
                                   levels = c("Alligator", 
                                              "Delaware",
                                              "Waccamaw",
                                              "SF"))
                        
png("FinalFigs/Figure5.png", width = 8, height = 4, units = "in", res = 300)
ggplot(frol$map_loaded, aes(Salt, CH4_ug_m2_h)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Salt)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, 
              aes(fill = Depth, shape = Estuary2, colour = Salt)) +
  geom_text(data = label_df2, aes(x, y, label = label), 
            size = 4, color = "black", inherit.aes = F) +
  labs(y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)"),
       colour = "Salinity",
       fill = "Depth (cm)",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(breaks = c("Waccamaw", "Alligator", "Delaware", "SF"), 
                     values = c(21, 22, 23, 24)) +
  guides(shape = guide_legend(order = 1,),
         colour = guide_legend(order = 2),
         fill = guide_legend(override.aes = list(shape = c(16, 1)), order = 3)) +
  facet_nested_wrap(~ Exp + Estuary2, scales = "free_y", nrow = 1) +
  theme_bw() +
  theme(legend.position = "right",
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0.25, 0, 0, -0.1, unit = "cm"),
        legend.key.size = unit(0.4, "cm"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()
  

# Correlations between CH4 and BGC and Guilds by each site
# Use the custom functions written to do this

# Get all corrs
sf_cor <- multiple_correlations(env_nona = env_nona_sf, var = "CH4_ug_m2_h") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sf_corG <- meth_corr_by_taxonomy(input = sf, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sf_corP <- meth_corr_by_taxonomy(input = sf, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sf_res <- rbind(sf_cor, sf_corG, sf_corP) %>%
  rename(SF_rho = rho,
         SF_sig = SpearmanPcut)

detra_cor <- multiple_correlations(env_nona = env_nona_detra, var = "CH4_ug_m2_h") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
detra_corG <- meth_corr_by_taxonomy(input = detra, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
detra_corP <- meth_corr_by_taxonomy(input = detra, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
de_res <- rbind(detra_cor, detra_corG, detra_corP) %>%
  rename(DE_rho = rho,
         DE_sig = SpearmanPcut)

nc_cor <- multiple_correlations(env_nona = env_nona_nc, var = "CH4_ug_m2_h") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
nc_corG <- meth_corr_by_taxonomy(input = nc, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
nc_corP <- meth_corr_by_taxonomy(input = nc, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
nc_res <- rbind(nc_cor, nc_corG, nc_corP) %>%
  rename(AL_rho = rho,
         AL_sig = SpearmanPcut)

sc_cor <- multiple_correlations(env_nona = env_nona_sc, var = "CH4_ug_m2_h") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sc_corG <- meth_corr_by_taxonomy(input = sc, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sc_corP <- meth_corr_by_taxonomy(input = sc, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sc_res <- rbind(sc_cor, sc_corG, sc_corP) %>%
  rename(WA_rho = rho,
         WA_sig = SpearmanPcut)

# Make df
# Prep. data frame
CH4_res <- data.frame("Variable" = c("CO2_ug_m2_h", "sed_per_org", "DOC_mgL", 
                                     "TOC_mgL", "sed_per_C", "sed_CN", 
                                     "sed_per_N", "DON_mgL", "DIN_mgL", "TN_mgL", 
                                     "NH4_mgL", "sed_NH4_mgL", "PO4_mgL", "sed_PO4_mgL",
                                     "pH", "Salinity_ppt_all", "NO3_mgL", "sed_NO3_mgL", 
                                     "N2O_ug_m2_h", "SRP_umol_m2_h", "SO4_mgL", 
                                     "sed_SO4_mgL", "Fe_mgL", "sed_Fe_mgL", "Mn_mgL", 
                                     "sed_Mn_mgL", "CH4_H2", "CH4_ac", "CH4_me", 
                                     "CH4_mix", "ANAMX", "AOA", "AOB", "NOB", "SRB_syn", 
                                     "Firmicutes", "ANME", "MOB_I", "MOB_II", "MOB_IIa", 
                                     "SRB", "FeRB"),
                      "Shortname" = c("CO2 Flux", "SOM", "DOC", "TOC", "C", "C:N", 
                                      "N", "DON", "DIN", "TN", "NH4", "NH4_s", "PO4",
                                      "PO4_s", "pH", "Salinity", "NO3", "NO3_s", 
                                      "N2O Flux", "SR", "SO4", "SO4_s", "Fe", "Fe_s", 
                                      "Mn", "Mn_s", "CH4_H2", "CH4_ac", "CH4_me",
                                      "CH4_mix", "ANAMX", "AOA", "AOB", "NOB",
                                      "SRB_syn", "Firmicutes", "ANME", "MOB_I",
                                      "MOB_II", "MOB_IIa", "SRB", "FeRB"),
                      "Type" = c("Flux", "Sediment", "Porewater", "Porewater", 
                                 "Sediment", "Sediment", "Sediment", "Porewater", 
                                 "Porewater", "Porewater", "Porewater", "Sediment", 
                                 "Porewater", "Sediment", "Porewater", "Porewater", 
                                 "Porewater", "Sediment", "Flux", "Porewater", 
                                 "Porewater", "Sediment", "Porewater", "Sediment", 
                                 "Porewater", "Sediment", "Sediment", "Sediment", 
                                 "Sediment", "Sediment", "Sediment", "Sediment", 
                                 "Sediment", "Sediment", "Sediment", "Sediment", 
                                 "Sediment", "Sediment", "Sediment", "Sediment", 
                                 "Sediment", "Sediment"),
                      "Prediction" = c("Positive", "Positive", "Positive", "Positive", 
                                       "Positive", "Negative", "Positive", "Positive", 
                                       "Positive", "Positive", "Positive", "Positive", 
                                       "Positive", "Positive", "Positive", "Negative", 
                                       "Negative", "Negative", "Negative", "Negative", 
                                       "Negative", "Negative", "Negative", "Negative", 
                                       "Negative", "Negative", "Positive", "Positive", 
                                       "Positive", "Positive", "Positive", "Positive", 
                                       "Positive", "Positive", "Positive", "Positive", 
                                       "Negative", "Negative", "Negative", "Negative", 
                                       "Negative", "Negative")) %>%
  left_join(., sf_res, by = "Variable") %>%
  left_join(., de_res, by = "Variable") %>%
  left_join(., nc_res, by = "Variable") %>%
  left_join(., sc_res, by = "Variable")

# Pretty heatmap
CH4_res_meta <- CH4_res %>%
  dplyr::select(Shortname, Type, Prediction, SF_sig, DE_sig, AL_sig, WA_sig) %>%
  mutate_if(is.character, as.factor)
CH4_res_mat <- CH4_res %>%
  column_to_rownames(var = "Shortname") %>%
  dplyr::select(SF_rho, DE_rho, AL_rho, WA_rho) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(CH4_res_mat), 
                       WA_sig = CH4_res_meta$WA_sig,
                       AL_sig = CH4_res_meta$AL_sig,
                       DE_sig = CH4_res_meta$DE_sig,
                       SF_sig = CH4_res_meta$SF_sig,
                       Prediction = CH4_res_meta$Prediction,
                       Type = CH4_res_meta$Type)
ann_colors <- list(Type = c(Flux = "#FFFF99", 
                            Porewater = "#A6CEE3", 
                            Sediment = "#B15928"),
                   Prediction = c(Positive = "red", 
                                  Negative = "blue"),
                   SF_sig = c(`Pfdr < 0.05` = "black", 
                              `Pfdr > 0.05` = "white"),
                   DE_sig = c(`Pfdr < 0.05` = "black", 
                              `Pfdr > 0.05` = "white"),
                   AL_sig = c(`Pfdr < 0.05` = "black", 
                              `Pfdr > 0.05` = "white"),
                   WA_sig = c(`Pfdr < 0.05` = "black", 
                              `Pfdr > 0.05` = "white"))
pheatmap(CH4_res_mat,
         legend = T,
         legend_breaks = c(-0.5, 0, 0.5, max(na.omit(CH4_res_mat))),
         legend_labels = c("-0.5", "0", "0.5", "rho\n"),
         main = "",
         #color = bluered(100),
         #border_color = NA,
         #na_col = "gray",
         scale = "none",
         angle_col = 315,
         fontsize = 10,
         fontsize_row = 10,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = c(14, 15, 26, 36),
         filename = "FinalFigs/Figure6.png",
         width = 5,
         height = 7)
dev.off()

show_col(brewer_pal(palette = "Paired")(12))
brewer_pal(palette = "Paired")(12)



#### 3. Lab Exps ####
# A closer look at these two experiments with contrasting methane/salinity results
# Highlight differences in phyla, guilds, methanogens, BGC
# Make Figure 6 multipanel figure
lab <- filter_data(frol,
                   filter_cat = "Site",
                   keep_vals = c("Soil incubation", "Soil Incubation")) # 23 samples
table(lab$map_loaded$Estuary)

lab$map_loaded <- lab$map_loaded %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Salt = recode_factor(Detail,
                              "Freshwater 12wk" = "Control",
                              "DI_ctrl" = "Control",
                              "5 ppt ASW 12wk" = "+ASW",
                              "5ppt ASW" = "+ASW"))
lab$map_loaded$EstSalt <- paste(lab$map_loaded$Estuary,
                                lab$map_loaded$Salt,
                                sep = "_")

# Phyla
lab_phyla <- summarize_taxonomy(lab, level = 2, report_higher_tax = F)
lab_barsP <- plot_taxa_bars(lab_phyla, lab$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
lab_legP <- get_legend(ggplot(lab_barsP, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)))
labP <- ggplot(lab_barsP, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.position = "none")
labP

# Guilds
lab_guilds <- summarize_taxonomy(lab, level = 9, report_higher_tax = F)
lab_barsGu <- plot_taxa_bars(lab_guilds, lab$map_loaded, "EstSalt",
                             num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
tallest_barG <- lab_barsGu %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
lab_legG <- get_legend(ggplot(lab_barsGu, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_barG$sum)/100, max(tallest_barG$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)))
labG <- ggplot(lab_barsGu, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none")
labG

# Methanobacteriaceae
lab_Methanobacteriaceae <- filter_taxa_from_input(lab,
                                                  taxa_to_keep = c("Methanobacteriaceae"),
                                                  at_spec_level = 5)
Methanobacteriaceae_wTax <- summarize_taxonomy(lab_Methanobacteriaceae, 
                                               level = 8, 
                                               report_higher_tax = T, 
                                               relative = FALSE)
rownames(Methanobacteriaceae_wTax) <- substring(rownames(Methanobacteriaceae_wTax), 62)
rownames(Methanobacteriaceae_wTax) <- substring(rownames(Methanobacteriaceae_wTax), 22)
rownames(Methanobacteriaceae_wTax) <- str_replace(rownames(Methanobacteriaceae_wTax), "ASV", "OTU")
nb.cols <- nrow(Methanobacteriaceae_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
lab_barsMethanobacteriaceae <- plot_taxa_bars(Methanobacteriaceae_wTax, 
                                              lab$map_loaded, 
                                              "EstSalt", 
                                              num_taxa = nrow(Methanobacteriaceae_wTax), 
                                              data_only = T) %>%
  mutate(mean_value = mean_value/26429) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
tallest_barM <- lab_barsMethanobacteriaceae %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
lab_legM <- get_legend(ggplot(lab_barsMethanobacteriaceae, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Methanobacteriaceae\nGenus; Species; OTU ID") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_barM$sum)/100, max(tallest_barM$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)))

labM <- ggplot(lab_barsMethanobacteriaceae, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Methanobacteriaceae\nGenus; Species; OTU ID") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none")

p <- plot_grid(labP, labG, labM, ncol = 1, rel_heights = c(0.35, 0.3, 0.35))
l <- plot_grid(lab_legP, lab_legG, lab_legM, ncol = 1, align = "v", axis = "b",
               rel_heights = c(0.4, 0.33, 0.5))
png("FinalFigs/Figure7.png", width = 5, height = 7, units = "in", res = 300)
plot_grid(p, l, ncol = 2, align = "h", rel_widths = c(0.57, 0.43))
dev.off()


