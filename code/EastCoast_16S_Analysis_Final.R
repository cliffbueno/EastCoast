# East Coast/SF 16S final data analysis
# by Cliff Bueno de Mesquita, Tringe Lab, JGI, Fall 2022 - Spring 2023
# Minor updates following reviewer comments, Spring 2024
# Samples from SF Bay, Delaware River, SC, NC
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
library(asbio) # Partial R2s
library(picante) # Trees
library(iCAMP) # NTI
library(lsr) # Cohen's d effect size
library(ape) # Phylogenetics
library(usmap) # Maps
library(ggrepel) # repel text

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
`%notin%` <- Negate(`%in%`)

# Guild subsetting module from other repository
source("~/Documents/GitHub/SF_microbe_methane/modules/3_OTU_subsetting_modules_v.0.4_strip.r")

# Correlations
source("~/Documents/GitHub/EastCoast/code/meth_corr_by_taxonomy.R")
source("~/Documents/GitHub/EastCoast/code/meth_corr_by_bgc.R")
source("~/Documents/GitHub/EastCoast/code/multiple_correlations.R")

# Plotting
source("~/Documents/GitHub/EastCoast/code/cliffplot_taxa_bars.R")
source("~/Documents/GitHub/EastCoast/code/plot_venn_diagram2.R")

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


#### _Process ####
# Input, filter, and rarefy (done once, now can input processed file at at 2. Combined)
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
#### _Start here ####
frol <- readRDS("data/frol.rds")
frol$map_loaded$EstSalt <- factor(frol$map_loaded$EstSalt,
                                  levels = c("Alligator_Oligohaline", "Alligator_Freshwater",
                                            "Delaware_Oligohaline", "Delaware_Freshwater", 
                                            "SF_Oligohaline", "SF_Freshwater",
                                            "Waccamaw_Oligohaline", "Waccamaw_Freshwater"))
# Add Lab/Field/Observational variable
frol$map_loaded$Exp <- "NA"
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Delaware" & 
      frol$map_loaded$Site[i] == "Soil mesocosm") {
    frol$map_loaded$Exp[i] <- "Field"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Delaware" & 
      frol$map_loaded$Site[i] == "Soil incubation") {
    frol$map_loaded$Exp[i] <- "Lab"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "SF") {
    frol$map_loaded$Exp[i] <- "Obs"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Waccamaw") {
    frol$map_loaded$Exp[i] <- "Field"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Alligator") {
    frol$map_loaded$Exp[i] <- "Lab"
  }
}
table(frol$map_loaded$Estuary)
frol$map_loaded$Exp <- factor(frol$map_loaded$Exp,
                              levels = c("Obs", "Field", "Lab"))
frol$map_loaded$Estuary2 <- factor(frol$map_loaded$Estuary,
                                   levels = c("SF", 
                                              "Waccamaw",
                                              "Delaware",
                                              "Alligator"))
frol$map_loaded$ExpEstSalt <- paste(frol$map_loaded$Exp, frol$map_loaded$EstSalt, sep = "_")
frol$map_loaded$ExpEstSalt <- factor(frol$map_loaded$ExpEstSalt,
                                  levels = c("Lab_Alligator_Oligohaline", "Lab_Alligator_Freshwater",
                                             "Lab_Delaware_Oligohaline", "Lab_Delaware_Freshwater", 
                                             "Field_Delaware_Oligohaline", "Field_Delaware_Freshwater",
                                             "Obs_SF_Oligohaline", "Obs_SF_Freshwater",
                                             "Field_Waccamaw_Oligohaline", "Field_Waccamaw_Freshwater"))
# Estuary with DE split into lab and field
frol$map_loaded$Estuary3 <- NA
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Delaware" & 
      frol$map_loaded$Site[i] == "Soil mesocosm") {
    frol$map_loaded$Estuary3[i] <- "DE Field"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Delaware" & 
      frol$map_loaded$Site[i] == "Soil incubation") {
    frol$map_loaded$Estuary3[i] <- "DE Lab"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "SF") {
    frol$map_loaded$Estuary3[i] <- "SF"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Waccamaw") {
    frol$map_loaded$Estuary3[i] <- "Waccamaw"
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$Estuary[i] == "Alligator") {
    frol$map_loaded$Estuary3[i] <- "Alligator"
  }
}
table(frol$map_loaded$Estuary3)
frol$map_loaded$Estuary3 <- factor(frol$map_loaded$Estuary3,
                                   levels = c("SF", "Waccamaw", "DE Field", "DE Lab", "Alligator"))
frol$map_loaded$Est3Salt <- paste(frol$map_loaded$Estuary3, frol$map_loaded$Salt, sep = "_")
frol$map_loaded$Est3Salt <- factor(frol$map_loaded$Est3Salt,
                                   levels = c("Alligator_Freshwater",
                                              "Alligator_Oligohaline",
                                              "DE Lab_Freshwater",
                                              "DE Lab_Oligohaline",
                                              "DE Field_Freshwater",
                                              "DE Field_Oligohaline",
                                              "SF_Freshwater",
                                              "SF_Oligohaline",
                                              "Waccamaw_Freshwater",
                                              "Waccamaw_Oligohaline"))
levels(frol$map_loaded$Est3Salt)

# Add MG:MT ratio to metadata
frol_guilds <- summarize_taxonomy(input = frol, level = 9, report_higher_tax = F) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sampleID = rownames(.)) %>%
  mutate(AO_NOB = (AOA + AOB) / NOB) %>%
  mutate(MG = CH4_ac + CH4_H2 + CH4_me + CH4_mix) %>%
  mutate(MT = MOB_I + MOB_II + MOB_IIa + ANME) %>%
  mutate(MG_MT = (CH4_ac + CH4_H2 + CH4_me + CH4_mix)/(ANME + MOB_I + MOB_II + MOB_IIa)) %>%
  left_join(., frol$map_loaded, by = "sampleID")

frol$map_loaded$MG_MT <- frol_guilds$MG_MT

# Add NC C, N, CN. Got data later so adding here!
nc_cn <- read_excel("data/Copy of CHN data.xls", sheet = 2) %>%
  slice(12:91) %>%
  dplyr::select(Name, `%N`, `%C`, `C:N`) %>%
  mutate(Name = gsub("bottom", "bot", Name)) %>%
  separate(Name, into = c("ID", "Depth"), sep = " ", remove = F) %>%
  mutate(Treatment = substr(ID, start = 1, stop = 1),
         Hydro = substr(ID, start = 2, stop = 2),
         Replicate = substr(ID, start = 3, stop = 3)) %>%
  filter(Hydro == "F") %>%
  mutate(Treatment = recode_factor(Treatment,
                                   "A" = "Control",
                                   "B" = "ASW",
                                   "C" = "ASW_noS",
                                   "D" = "SO4")) %>%
  filter(Treatment %in% c("Control", "ASW")) %>%
  arrange(desc(Depth), desc(Treatment)) %>%
  filter(Name %notin% c("BF1 top", "BF2 top", "AF1 top", "BF2 bot", "AF2 bot"))
rownames(frol$map_loaded)
frol$map_loaded$sed_per_C[47:61] <- nc_cn$`%C`
frol$map_loaded$sed_per_N[47:61] <- nc_cn$`%N`
frol$map_loaded$sed_CN[47:61] <- nc_cn$`C:N`

# Add Delaware lab pH data (Nat found much later on!)
de_ph <- read_xlsx("data/Delaware Lab Expt pH.xlsx", sheet = 2) %>%
  arrange(sampleID)
sum(rownames(frol$map_loaded)[9:16] != de_ph$sampleID) # good
frol$map_loaded$pH[9:16] <- de_ph$pH

# Add Delaware CH4 data
# Nat sent whole-core integrated values. Add duplicates to D1 and D2 like the other datasets
# When plotting, only plot one depth though.
# Don't use other lab biogeochem data because they are all integrated and not separated by depth
# Values would have to be repeated by depth and this is not accurate, better to leave out.
# And actually, can't even back convert from integrated to concentration because don't have datasheet.
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$sampleID[i] == "TS_FW_d1_12_1") {
    frol$map_loaded$CH4_ug_m2_h[i] <- 20303.49838
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$sampleID[i] == "TS_FW_d2_12_1") {
    frol$map_loaded$CH4_ug_m2_h[i] <- 20303.49838
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$sampleID[i] == "TS_FW_d1_12_2") {
    frol$map_loaded$CH4_ug_m2_h[i] <- 84770.97026
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$sampleID[i] == "TS_FW_d2_12_2") {
    frol$map_loaded$CH4_ug_m2_h[i] <- 84770.97026
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$sampleID[i] == "TS_SW_d1_12_1") {
    frol$map_loaded$CH4_ug_m2_h[i] <- 32976.95881
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$sampleID[i] == "TS_SW_d2_12_1") {
    frol$map_loaded$CH4_ug_m2_h[i] <- 32976.95881
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$sampleID[i] == "TS_SW_d1_12_2") {
    frol$map_loaded$CH4_ug_m2_h[i] <- 68221.38115
  }
}
for (i in 1:nrow(frol$map_loaded)) {
  if (frol$map_loaded$sampleID[i] == "TS_SW_d2_12_2") {
    frol$map_loaded$CH4_ug_m2_h[i] <- 68221.38115
  }
}

# Export this metadata
#write.csv(frol$map_loaded, "data/metadata_used.csv")

# Subsets
sf <- filter_data(frol,
                  filter_cat = "Estuary",
                  keep_vals = "SF")
de <- filter_data(frol,
                  filter_cat = "Estuary",
                  keep_vals = "Delaware")
detra <- filter_data(de,
                     filter_cat = "Site",
                     keep_vals = "Soil mesocosm")
sc <- filter_data(frol,
                  filter_cat = "Estuary",
                  keep_vals = "Waccamaw")
deinc <- filter_data(de,
                     filter_cat = "Site",
                     keep_vals = "Soil incubation")
nc <- filter_data(frol,
                  filter_cat = "Estuary",
                  keep_vals = "Alligator")

# Get environmental variables for each
env <- frol$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, Salinity_ppt_all, pH)
env_nona <- na.omit(env) # n = 110. 102 with pH.

env_sf <- sf$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, 
                Salinity_ppt_all, SO4_mgL, DOC_mgL, Fe_mgL, Mn_mgL, Cu_mgL, Zn_mgL,
                sed_pH, sed_NH4_mgL, sed_NO3_mgL, sed_PO4_mgL, sed_Cl_mgL, sed_SO4_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_Bulk_dens,
                sed_Fe_mgL, sed_Mn_mgL, sed_Cu_mgL, sed_Zn_mgL)
env_nona_sf <- na.omit(env_sf) # n = 30 (just D2)

env_detra <- detra$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, N2O_ug_m2_h, 
                Salinity_ppt_all, 
                NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL,
                Fe_mgL, Porosity, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d)
env_nona_detra <- na.omit(env_detra) # n = 8

env_deinc <- deinc$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, 
                #Salinity_ppt_all, # No salinity because in mmol/m2 and can't convert
                pH)
env_nona_deinc <- na.omit(env_deinc) # n = 8

env_sc <- sc$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, CH4_pw_air_ppmv,
                CH4_pot_umol_gdw_h, CO2_pot_umol_gdw_h,
                N2_umol_m2_h,	SOD_umol_m2_h, NO3_umol_m2_h,	NH4_umol_m2_h, SRP_umol_m2_h, DON_umol_m2_h,
                Conductivity_uS_cm, Salinity_ppt_all, DIC_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_per_org, sed_per_inorg, pH)
env_nona_sc <- na.omit(env_sc) # n = 15

env_nc <- nc$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, N2O_ug_m2_h, Salinity_ppt_all, 
                TOC_mgL, TN_mgL, NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL, Br_mgL, NO3_mgL,
                DIN_mgL, DON_mgL, pH, sed_per_C, sed_per_N, sed_CN)
env_nona_nc <- na.omit(env_nc) # n = 15



#### _Alpha ####
leveneTest(frol$map_loaded$rich ~ frol$map_loaded$Estuary) # Homogeneous
leveneTest(frol$map_loaded$rich ~ frol$map_loaded$Salt) # Homogeneous
leveneTest(frol$map_loaded$rich ~ frol$map_loaded$Depth) # Homogeneous
m <- aov(rich ~ Estuary + Salt + Depth, data = frol$map_loaded)
shapiro.test(m$residuals) # Not normally distributed
Anova(m, type = "II") # All sig.
m <- aov(rich ~ ExpEstSalt, data = frol$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "ExpEstSalt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(frol$map_loaded$rich)+(max(frol$map_loaded$rich)-min(frol$map_loaded$rich))/20)
all <- aov(rich ~ Estuary + Salt + Depth, data = frol$map_loaded)
noEst <- aov(rich ~ Salt + Depth, data = frol$map_loaded)
noSal <- aov(rich ~ Estuary + Depth, data = frol$map_loaded)
noDep <- aov(rich ~ Estuary + Salt, data = frol$map_loaded)
partial.R2(noEst, all)
partial.R2(noSal, all)
partial.R2(noDep, all)

leveneTest(frol$map_loaded$shannon ~ frol$map_loaded$Estuary) # Almost homogeneous
leveneTest(frol$map_loaded$shannon ~ frol$map_loaded$Salt) # Homogeneous
leveneTest(frol$map_loaded$shannon ~ frol$map_loaded$Depth) # Homogeneous
m1 <- aov(shannon ~ Estuary + Salt + Depth, data = frol$map_loaded)
shapiro.test(m1$residuals) # Not normally distributed but using anyway
Anova(m1, type = "II") # All sig
m1 <- aov(shannon ~ ExpEstSalt, data = frol$map_loaded)
shapiro.test(m1$residuals) # Not normal
summary(m1)
t1 <- emmeans(object = m1, specs = "ExpEstSalt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(frol$map_loaded$shannon)+(max(frol$map_loaded$shannon)-min(frol$map_loaded$shannon))/20)
all <- aov(shannon ~ Estuary + Salt + Depth, data = frol$map_loaded)
noEst <- aov(shannon ~ Salt + Depth, data = frol$map_loaded)
noSal <- aov(shannon ~ Estuary + Depth, data = frol$map_loaded)
noDep <- aov(shannon ~ Estuary + Salt, data = frol$map_loaded)
partial.R2(noEst, all)
partial.R2(noSal, all)
partial.R2(noDep, all)

# Combined Plot
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- frol$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
png("FinalFigs/Figure2.png", width = 7, height = 4.5, units = "in", res = 300)
ggplot(alpha_long, aes(ExpEstSalt, value)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Salt)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, 
              aes(fill = Depth, shape = Estuary, colour = Salt)) +
  geom_text(data = label_df, aes(ExpEstSalt, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site",
       colour = "Salinity",
       fill = "Depth (cm)",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(breaks = c( "SF", "Waccamaw", "Delaware", "Alligator"), 
                     values = c(24, 21, 23, 22)) +
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
        axis.text.x = element_text(size = 9, angle = 315, hjust = 0),
        strip.text = element_text(size = 10),
        plot.margin = margin(5, 10, 5, 5, "pt"))
dev.off()

# Remake without barplot and extract legend for beta diversity plot
good_leg_plot <- ggplot(alpha_long, aes(ExpEstSalt, value)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, 
              aes(fill = Depth, shape = Estuary, colour = Salt)) +
  geom_text(data = label_df, aes(ExpEstSalt, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site",
       colour = "Salinity",
       fill = "Depth (cm)",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(breaks = c( "SF", "Waccamaw", "Delaware", "Alligator"), 
                     values = c(24, 21, 23, 22)) +
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
        plot.margin = margin(5, 5, 5, 30, "pt"))
good_leg <- get_legend(good_leg_plot)



#### _Beta  ####
# Bray-Curtis dissimilarity
bc <- calc_dm(frol$data_loaded)

# Stats
set.seed(100)
adonis2(bc ~ frol$map_loaded$Estuary + frol$map_loaded$Salt + frol$map_loaded$Depth) # Salt and Depth sig
anova(betadisper(bc, frol$map_loaded$Estuary)) # Dispersion not homogeneous
anova(betadisper(bc, frol$map_loaded$Salt)) # Dispersion homogeneous
anova(betadisper(bc, frol$map_loaded$Depth)) # Dispersion homogeneous

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
  mutate(shortnames = c("CH4", "Salinity", "pH"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
frol$map_loaded$Axis01 <- scores(pcoa)[,1]
frol$map_loaded$Axis02 <- scores(pcoa)[,2]

# Plot
micro.hulls <- ddply(frol$map_loaded, c("EstSalt"), find_hull)
g <- ggplot(frol$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = EstSalt), fill = NA,
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(fill = Depth, colour = Salt, shape = Estuary),
             show.legend = T) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 2, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salt",
       fill = "Depth (cm)",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red", "blue", "red", "blue", "red",
                                 "blue", "red", "blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_y_continuous(expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  guides(shape = guide_legend(order = 1,),
         colour = guide_legend(order = 2),
         fill = guide_legend(override.aes = list(shape = c(16, 1)), order = 3)) +
  theme_bw() +  
  ggtitle("(a) Combined") +
  theme(legend.position = "right",
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0.25, 0, 0, -0.1, unit = "cm"),
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
# leg <- get_legend(g) # use good_leg from alpha diversity plot
g <- g + theme(legend.position = "none")
g

# Ordistep
comm_nona <- as.data.frame(t(frol$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # Salinity, CH4

# Make other PCoAs individually
bc_sf <- calc_dm(sf$data_loaded)
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
  geom_text_repel(data = vec.df_sf,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 2, color = "black", box.padding = 0.05) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  theme_bw() +
  ggtitle("(b) SF (field obs.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g2
comm_nona_sf <- as.data.frame(t(sf$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_sf))
mod0 <- rda(comm_nona_sf ~ 1, env_nona_sf)  # Model with intercept only
mod1 <- rda(comm_nona_sf ~ ., env_nona_sf)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # BD, sed_Zn, sed_pH, sed_Cl, sed_CN, sed_Mn, SO4

bc_detra <- calc_dm(detra$data_loaded)
pcoa_detra <- cmdscale(bc_detra, k = nrow(detra$map_loaded) - 1, eig = T)
set.seed(100)
ef_detra <- envfit(pcoa_detra, env_detra, permutations = 999, na.rm = TRUE)
ef_detra
ordiplot(pcoa_detra)
plot(ef_detra, p.max = 0.05, cex = 0.5)
manual_factor_detra <- 0.17
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
  geom_text_repel(data = vec.df_detra,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 2, color = "black", box.padding = 0.05) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  theme_bw() +
  ggtitle("(d) Delaware (field exp.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g3
comm_nona_detra <- as.data.frame(t(detra$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_detra))
mod0 <- rda(comm_nona_detra ~ 1, env_nona_detra)  # Model with intercept only
mod1 <- rda(comm_nona_detra ~ ., env_nona_detra)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # NH4

bc_sc <- calc_dm(sc$data_loaded)
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
  geom_text_repel(data = vec.df_sc,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 2, color = "black", box.padding = 0.05) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  theme_bw() +
  ggtitle("(c) Waccamaw (field exp.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g4
comm_nona_sc <- as.data.frame(t(sc$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_sc))
mod0 <- rda(comm_nona_sc ~ 1, env_nona_sc)  # Model with intercept only
mod1 <- rda(comm_nona_sc ~ ., env_nona_sc)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # SOD, CO2

# Note: in initial analysis, thought there was error, but looks like maybe there isn't
# Just not clear separation, no hard evidence of any sample label error
# Note no env. data here
# Update - now have pH data. Check.
bc_deinc <- calc_dm(deinc$data_loaded)
pcoa_deinc <- cmdscale(bc_deinc, k = nrow(deinc$map_loaded) - 1, eig = T)
set.seed(100)
ef_deinc <- envfit(pcoa_deinc, env_deinc, permutations = 999, na.rm = TRUE)
ef_deinc
# CH4 and pH not significant, don't plot
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
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  theme_bw() +
  ggtitle("(e) Delaware (lab exp.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g5
comm_nona_deinc <- as.data.frame(t(deinc$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_deinc))
mod0 <- rda(comm_nona_deinc ~ 1, env_nona_deinc)  # Model with intercept only
mod1 <- rda(comm_nona_deinc ~ ., env_nona_deinc)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # Can't do yet, no data

bc_nc <- calc_dm(nc$data_loaded)
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
  geom_text_repel(data = vec.df_nc,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 2, color = "black", box.padding = 0.05) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  theme_bw() +
  ggtitle("(f) Alligator (lab exp.)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -0.5, size = 12),
        plot.margin = margin(0,0,0,10, "pt"))
g6
comm_nona_nc <- as.data.frame(t(nc$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_nc))
mod0 <- rda(comm_nona_nc ~ 1, env_nona_nc)  # Model with intercept only
mod1 <- rda(comm_nona_nc ~ ., env_nona_nc)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # SO4, CH4, PO4

pcoa_plots <- plot_grid(g, g2, g4, g3, g5, g6, ncol = 2)
fig3 <- plot_grid(pcoa_plots, good_leg, rel_widths = c(0.85, 0.15))
fig3
png("FinalFigs/Figure3.png", width = 8, height = 6, units = "in", res = 300)
fig3
dev.off()



#### _Taxa ####
cliffplot_taxa_bars(input = frol, level = 1, variable = "Est3Salt")
cliffplot_taxa_bars(input = frol, level = 2, variable = "Est3Salt")
cliffplot_taxa_bars(input = frol, level = 3, variable = "Est3Salt")
cliffplot_taxa_bars(input = frol, level = 4, variable = "Est3Salt")
cliffplot_taxa_bars(input = frol, level = 5, variable = "Est3Salt")
cliffplot_taxa_bars(input = frol, level = 6, variable = "Est3Salt")
cliffplot_taxa_bars(input = frol, level = 9, variable = "Est3Salt")

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
topphy <- barsP %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean)
facet_names <- c("Obs" = "Observational", 
                 "Field" = "Field",
                 "Lab" = "Lab",
                 "Waccamaw" = "Waccamaw",
                 "Alligator" = "Alligator",
                 "Delaware" = "Delaware",
                 "SF" = "SF",
                 "Freshwater" = "Fresh",
                 "Oligohaline" = "Oligo")
phy <- ggplot(barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_nested(~ Exp + Estuary2 + Salt, space = "free", scales = "free_x", 
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 25, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 25, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
phy_leg <- get_legend(phy)
phy <- phy + theme(legend.position = "none")
phy

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
topgui <- barsG %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean)
topgui_t <- topgui %>%
  column_to_rownames(var = "taxon") %>%
  t() %>%
  as.data.frame() %>%
  mutate(MG = CH4_H2 + CH4_me + CH4_ac + CH4_mix) %>%
  mutate(MT = MOB_I + MOB_II + MOB_IIa + ANME)
gui <- ggplot(barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Exp + Estuary2 + Salt, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 10, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 10, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
gui_leg <- get_legend(gui)
gui <- gui + theme(legend.position = "none")
gui

pg <- plot_grid(phy, gui, ncol = 1)
ll <- plot_grid(phy_leg, gui_leg, ncol = 1)

plot_grid(pg, ll, rel_widths = c(0.75, 0.25))

#png("Figure3_forPPT.png", width = 6.5, height = 6, units = "in", res = 300)
pg
#dev.off()

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

# DE combined
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

# DE Field
tax_sum_phyla_detra <- summarize_taxonomy(input = detra, level = 2, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsP$taxon)
res_phy_detra <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla_detra,
                                             metadata_map = detra$map_loaded, 
                                             type_header = 'Salt', 
                                             test_type = 'MW') %>%
  mutate(Estuary = "Delaware Field",
         Level = "Phylum",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_guilds_detra <- summarize_taxonomy(input = detra, level = 9, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsG$taxon)
res_gui_detra <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_guilds_detra,
                                             metadata_map = de$map_loaded, 
                                             type_header = 'Salt', 
                                             test_type = 'MW') %>%
  mutate(Estuary = "Delaware Field",
         Level = "Guild",
         Taxon = rownames(.)) %>%
  arrange(match(Taxon, Guild_cols$Guild)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# DE Lab
tax_sum_phyla_deinc <- summarize_taxonomy(input = deinc, level = 2, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsP$taxon)
res_phy_deinc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla_deinc,
                                             metadata_map = de$map_loaded, 
                                             type_header = 'Salt', 
                                             test_type = 'MW') %>%
  mutate(Estuary = "Delaware Lab",
         Level = "Phylum",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_guilds_deinc <- summarize_taxonomy(input = deinc, level = 9, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsG$taxon)
res_gui_deinc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_guilds_deinc,
                                             metadata_map = de$map_loaded, 
                                             type_header = 'Salt', 
                                             test_type = 'MW') %>%
  mutate(Estuary = "Delaware Lab",
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

res_phy <- res_phy_sc %>%
  mutate(WA = SaltEffect) %>%
  mutate(AL = res_phy_nc$SaltEffect,
         DEF = res_phy_detra$SaltEffect,
         DEL = res_phy_deinc$SaltEffect,
         SF = res_phy_sf$SaltEffect) %>%
  dplyr::select(Level, Taxon, SF, WA, DEF, DEL, AL)


res_gui <- res_gui_sc %>%
  mutate(WA = SaltEffect) %>%
  mutate(AL = res_gui_nc$SaltEffect,
         DEF = res_gui_detra$SaltEffect,
         DEL = res_gui_deinc$SaltEffect,
         SF = res_gui_sf$SaltEffect) %>%
  dplyr::select(Level, Taxon, SF, WA, DEF, DEL, AL)

res <- rbind(res_phy, res_gui)
#write_xlsx(res, "taxa_results_unformatted.xlsx", format_headers = F)
brewer.pal(12, "Paired")



#### __AO:NOB and MG:MT ####
frol_guilds <- summarize_taxonomy(input = frol, level = 9, report_higher_tax = F) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sampleID = rownames(.)) %>%
  mutate(AO_NOB = (AOA + AOB) / NOB) %>%
  mutate(MG = CH4_ac + CH4_H2 + CH4_me + CH4_mix) %>%
  mutate(MT = MOB_I + MOB_II + MOB_IIa + ANME) %>%
  mutate(MG_MT = (CH4_ac + CH4_H2 + CH4_me + CH4_mix)/(ANME + MOB_I + MOB_II + MOB_IIa)) %>%
  left_join(., frol$map_loaded, by = "sampleID")

# Check sig., make sig. solid and not sig. dotted
summary(lm(log(MG_MT) ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "SF")))
summary(lm(log(MG_MT) ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "Waccamaw")))
summary(lm(log(MG_MT) ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "DE Field")))
summary(lm(log(MG_MT) ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "DE Lab")))
summary(lm(log(MG_MT) ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "Alligator")))

a <- ggplot(frol_guilds, aes(AO_NOB, MG_MT, colour = Estuary3, shape = Estuary3)) +
  geom_point(size = 2) +
  geom_smooth(data = subset(frol_guilds, Estuary3 == "SF" | Estuary3 == "Waccamaw"),
              method = "lm", size = 0.25, alpha = 0.1) +
  geom_smooth(data = subset(frol_guilds, Estuary3 == "DE Field" | 
                              Estuary3 == "DE Lab" | Estuary3 == "Alligator"),
              method = "lm", size = 0.5, alpha = 0.1, linetype = "dotted", se = F) +
  labs(x = "Ammonia oxidizers: Nitrite oxidizing bacteria",
       y = "Methanogens: Methanotrophs",
       colour = "Site/Experiment",
       shape = "Site/Experiment") +
  scale_shape_manual(values = c(17, 16, 18, 18, 15)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c(0.01, 0.1, 1, 10, 100)) +
  guides(colour = guide_legend(shape = c(17, 16, 18, 18, 15))) +
  theme_bw()
# This is the opposite of what we should see!
l <- get_legend(a)
a <- a + theme(legend.position = "none",
               axis.title = element_text(size = 10))

# Check sig., make sig. solid and not sig. dotted
summary(lm(MT ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "SF")))
summary(lm(MT ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "Waccamaw")))
summary(lm(MT ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "DE Field")))
summary(lm(MT ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "DE Lab")))
summary(lm(MT ~ log(AO_NOB), data = subset(frol_guilds, Estuary3 == "Alligator")))

b <- ggplot(frol_guilds, aes(AO_NOB, MT*100, colour = Estuary3, shape = Estuary3)) +
  geom_point(size = 2) +
  geom_smooth(data = subset(frol_guilds, Estuary3 == "SF"),
              method = "lm", size = 0.25, alpha = 0.1) +
  geom_smooth(data = subset(frol_guilds, Estuary3 == "DE Field" | 
                              Estuary3 == "DE Lab" | Estuary3 == "Alligator" 
                            | Estuary3 == "Waccamaw"),
              method = "lm", size = 0.5, alpha = 0.1, linetype = "dotted", se = F) +
  labs(x = "Ammonia oxidizers: Nitrite oxidizing bacteria",
       y = "Methanotroph % abundance",
       colour = "Site/Experiment",
       shape = "Site/Experiment") +
  scale_shape_manual(values = c(17, 16, 18, 18, 15)) +
  scale_x_continuous(trans = 'log10') +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 10))
# Variable!

png("FinalFigs/FigureS5.png", width = 8, height = 4, units = "in", res = 300)
plot_grid(a, b, l, ncol = 3, rel_widths = c(0.43, 0.43, 0.14), labels = c("a", "b", ""))
dev.off()



### __Methanogens ####
frol_mg <- filter_taxa_from_input(frol,
                                  taxa_to_keep = c("CH4_ac", "CH4_H2", "CH4_me", "CH4_mix"),
                                  at_spec_level = 9)
tax_sum_mg <- summarize_taxonomy(input = frol_mg, 
                                 level = 5, 
                                 report_higher_tax = F, 
                                 relative = F) %>%
  # filter(rownames(.) != "NA") %>%
  mutate_all(funs((./31264)*100))
barsMG <- plot_taxa_bars(tax_sum_mg,
                         frol_mg$map_loaded,
                         "sampleID",
                         num_taxa = 20,
                         data_only = TRUE) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  #mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., frol_mg$map_loaded, by = c("group_by" = "sampleID"))
topmg <- barsMG %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean)
barsMG <- barsMG %>%
  mutate(taxon = factor(taxon, levels = topmg$taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
facet_names <- c("Obs" = "Observational", 
                 "Field" = "Field",
                 "Lab" = "Lab",
                 "Waccamaw" = "Waccamaw",
                 "Alligator" = "Alligator",
                 "Delaware" = "Delaware",
                 "SF" = "SF",
                 "Freshwater" = "Fresh",
                 "Oligohaline" = "Oligo")
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
png("FigureS3_forPPT.png", width = 6.5, height = 6, units = "in", res = 300)
ggplot(barsMG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "% Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey90", mycolors)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_nested(~ Exp + Estuary2 + Salt, space = "free", scales = "free_x", 
               labeller = as_labeller(facet_names)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(5, 5, 5, 5, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "none")
dev.off()

# Venn
frol_mg$map_loaded$Site <- as.factor(frol_mg$map_loaded$Site)
plot_venn_diagram(frol_mg, 
                  "Estuary", 
                  0.00000000000000000000000000000000001)

# Subsets
sf_mg <- filter_data(frol_mg,
                     filter_cat = "Estuary",
                     keep_vals = "SF")
de_mg <- filter_data(frol_mg,
                     filter_cat = "Estuary",
                     keep_vals = "Delaware")
detra_mg <- filter_data(de_mg,
                        filter_cat = "Site",
                        keep_vals = "Soil mesocosm")
sc_mg <- filter_data(frol_mg,
                     filter_cat = "Estuary",
                     keep_vals = "Waccamaw")
deinc_mg <- filter_data(de_mg,
                        filter_cat = "Site",
                        keep_vals = "Soil incubation")
nc_mg <- filter_data(frol_mg,
                     filter_cat = "Estuary",
                     keep_vals = "Alligator")

# Stats
tax_sum_mg_sc <- summarize_taxonomy(input = sc_mg, level = 5, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMG$taxon)
res_mg_sc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mg_sc,
                                         metadata_map = sc_mg$map_loaded, 
                                         type_header = 'Salt', 
                                         test_type = 'MW') %>%
  mutate(Estuary = "Waccamaw",
         Level = "Family",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_mg_nc <- summarize_taxonomy(input = nc_mg, level = 5, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMG$taxon)
res_mg_nc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mg_nc,
                                         metadata_map = nc_mg$map_loaded, 
                                         type_header = 'Salt',
                                         test_type = 'MW') %>%
  mutate(Estuary = "Alligator",
         Level = "Family",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# DE Field
tax_sum_mg_detra <- summarize_taxonomy(input = detra_mg, level = 5, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMG$taxon)
res_mg_detra <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mg_detra,
                                            metadata_map = detra_mg$map_loaded, 
                                            type_header = 'Salt', 
                                            test_type = 'MW') %>%
  mutate(Estuary = "Delaware Field",
         Level = "Family",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# DE Lab
tax_sum_mg_deinc <- summarize_taxonomy(input = deinc_mg, level = 5, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMG$taxon)
res_mg_deinc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mg_deinc,
                                            metadata_map = deinc_mg$map_loaded, 
                                            type_header = 'Salt', 
                                            test_type = 'MW') %>%
  mutate(Estuary = "Delaware Lab",
         Level = "Family",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# SF
tax_sum_mg_sf <- summarize_taxonomy(input = sf_mg, level = 5, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMG$taxon)
res_mg_sf <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mg_sf,
                                         metadata_map = sf_mg$map_loaded, 
                                         type_header = 'Salt', 
                                         test_type = 'MW') %>%
  mutate(Estuary = "SF",
         Level = "Family",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# Full join, takes into account different taxa among the sites
# Order should be SF, WA, DEF, DEL, AL
topmg_noUC <- subset(topmg, taxon != "Unclassified") %>%
  mutate(ord = rev(seq(1:13)))
res_mg <- full_join(res_mg_sf, res_mg_sc, by = "Taxon", keep = F) %>%
  full_join(., res_mg_detra, by = "Taxon", keep = F) %>%
  full_join(., res_mg_deinc, by = "Taxon", keep = F) %>%
  full_join(., res_mg_nc, by = "Taxon", keep = F) %>%
  dplyr::select(Taxon, SaltEffect.x, SaltEffect.y, SaltEffect.x.x, SaltEffect.y.y, SaltEffect) %>%
  set_names(c("Taxon", "SF", "WA", "DEF", "DEL", "AL")) %>%
  replace_na(., replace = list(SF = "NP", WA = "NP", DEF = "NP", DEL = "NP", AL = "NP")) %>%
  left_join(., topmg_noUC, by = c("Taxon" = "taxon")) %>%
  arrange(ord) %>%
  dplyr::select(-mean, -ord)
#write_xlsx(res_mg, "mg_results_unformatted.xlsx", format_headers = F)
nrow(frol_mg$taxonomy_loaded)
table(frol_mg$taxonomy_loaded$taxonomy9)


### __Methanotrophs ####
frol_mt <- filter_taxa_from_input(frol,
                                  taxa_to_keep = c("ANME", "MOB_I", "MOB_II", "MOB_IIa"),
                                  at_spec_level = 9)
# Remove Methylotenera. Methylotroph but not CH4 oxidizer!
frol_mt <- filter_taxa_from_input(frol_mt,
                                  taxa_to_remove = "Methylotenera",
                                  at_spec_level = 6)
tax_sum_mt <- summarize_taxonomy(input = frol_mt, 
                                 level = 6, 
                                 report_higher_tax = F, 
                                 relative = F) %>%
  # filter(rownames(.) != "NA") %>%
  mutate_all(funs((./31264)*100))
barsMT <- plot_taxa_bars(tax_sum_mt,
                         frol_mt$map_loaded,
                         "sampleID",
                         num_taxa = 20,
                         data_only = TRUE) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., frol_mt$map_loaded, by = c("group_by" = "sampleID"))
topmt <- barsMT %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean)
barsMT <- barsMT %>%
  mutate(taxon = factor(taxon, levels = topmt$taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
facet_names <- c("Obs" = "Observational", 
                 "Field" = "Field",
                 "Lab" = "Lab",
                 "Waccamaw" = "Waccamaw",
                 "Alligator" = "Alligator",
                 "Delaware" = "Delaware",
                 "SF" = "SF",
                 "Freshwater" = "Fresh",
                 "Oligohaline" = "Oligo")
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
png("FigureS4_forPPT.png", width = 6.5, height = 6, units = "in", res = 300)
ggplot(barsMT, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "% Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0,10)) +  
  facet_nested(~ Exp + Estuary2 + Salt, space = "free", scales = "free_x", 
               labeller = as_labeller(facet_names)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(5, 5, 5, 5, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "none")
dev.off()

# Subsets
sf_mt <- filter_data(frol_mt,
                     filter_cat = "Estuary",
                     keep_vals = "SF")
de_mt <- filter_data(frol_mt,
                     filter_cat = "Estuary",
                     keep_vals = "Delaware")
detra_mt <- filter_data(de_mt,
                        filter_cat = "Site",
                        keep_vals = "Soil mesocosm")
sc_mt <- filter_data(frol_mt,
                     filter_cat = "Estuary",
                     keep_vals = "Waccamaw")
deinc_mt <- filter_data(de_mt,
                        filter_cat = "Site",
                        keep_vals = "Soil incubation")
nc_mt <- filter_data(frol_mt,
                     filter_cat = "Estuary",
                     keep_vals = "Alligator")

# Stats
tax_sum_mt_sc <- summarize_taxonomy(input = sc_mt, level = 6, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMT$taxon)
res_mt_sc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mt_sc,
                                          metadata_map = sc_mt$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "Waccamaw",
         Level = "Genus",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

tax_sum_mt_nc <- summarize_taxonomy(input = nc_mt, level = 6, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMT$taxon)
res_mt_nc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mt_nc,
                                          metadata_map = nc_mt$map_loaded, 
                                          type_header = 'Salt',
                                         test_type = 'MW') %>%
  mutate(Estuary = "Alligator",
         Level = "Genus",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# DE Field
tax_sum_mt_detra <- summarize_taxonomy(input = detra_mt, level = 6, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMT$taxon)
res_mt_detra <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mt_detra,
                                             metadata_map = detra_mt$map_loaded, 
                                             type_header = 'Salt', 
                                             test_type = 'MW') %>%
  mutate(Estuary = "Delaware Field",
         Level = "Genus",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# DE Lab
tax_sum_mt_deinc <- summarize_taxonomy(input = deinc_mt, level = 6, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMT$taxon)
res_mt_deinc <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mt_deinc,
                                             metadata_map = deinc_mt$map_loaded, 
                                             type_header = 'Salt', 
                                             test_type = 'MW') %>%
  mutate(Estuary = "Delaware Lab",
         Level = "Genus",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# SF
tax_sum_mt_sf <- summarize_taxonomy(input = sf_mt, level = 6, report_higher_tax = F) %>%
  filter(rownames(.) %in% barsMT$taxon)
res_mt_sf <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_mt_sf,
                                          metadata_map = sf_mt$map_loaded, 
                                          type_header = 'Salt', 
                                          test_type = 'MW') %>%
  mutate(Estuary = "SF",
         Level = "Genus",
         Taxon = rownames(.)) %>%
  arrange(desc(Taxon)) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Yes", "No")) %>%
  mutate(SaltEffect = ifelse(Sig == "Yes",
                             ifelse(Freshwater > Oligohaline, "-", "+"),
                             ""))

# Full join, takes into account different taxa among the sites
# Order should be SF, WA, DEF, DEL, AL
topmt_noUC <- subset(topmt, taxon != "Unclassified" & taxon != "Other") %>%
  mutate(ord = rev(seq(1:19)))
res_mt <- full_join(res_mt_sf, res_mt_sc, by = "Taxon", keep = F) %>%
  full_join(., res_mt_detra, by = "Taxon", keep = F) %>%
  full_join(., res_mt_deinc, by = "Taxon", keep = F) %>%
  full_join(., res_mt_nc, by = "Taxon", keep = F) %>%
  dplyr::select(Taxon, SaltEffect.x, SaltEffect.y, SaltEffect.x.x, SaltEffect.y.y, SaltEffect) %>%
  set_names(c("Taxon", "SF", "WA", "DEF", "DEL", "AL")) %>%
  replace_na(., replace = list(SF = "NP", WA = "NP", DEF = "NP", DEL = "NP", AL = "NP")) %>%
  left_join(., frol_mt$taxonomy_loaded, by = c("Taxon" = "taxonomy6")) %>%
  group_by(Taxon) %>%
  slice_head(n = 1) %>%
  left_join(., topmt_noUC, by = c("Taxon" = "taxon")) %>%
  arrange(ord) %>%
  dplyr::select(-mean, -ord, -taxonomy1, -taxonomy2, -taxonomy3, -taxonomy4, 
                -taxonomy5, -taxonomy7, -taxonomy8)
write_xlsx(res_mt, "mt_results_unformatted.xlsx", format_headers = F)



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
png("FinalFigs/Figure5.png", width = 5, height = 5, units = "in", res = 300)
plot_grid(top, frol.l, nrow = 2, rel_heights = c(0.87, 0.13))
dev.off()



#### _CH4 ####
# CH4 by salt and experiment
# Stats
leveneTest(CH4_ug_m2_h ~ Salt, data = sf$map_loaded)
shapiro.test(sf$map_loaded$CH4_ug_m2_h)
t.test(CH4_ug_m2_h ~ Salt, data = sf$map_loaded) # Sig.
wilcox.test(CH4_ug_m2_h ~ Salt, data = sf$map_loaded) # Sig.

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
t.test(CH4_ug_m2_h ~ Salt, data = deinc$map_loaded) # No data but sig.
wilcox.test(CH4_ug_m2_h ~ Salt, data = deinc$map_loaded) # No data but sig.

leveneTest(CH4_ug_m2_h ~ Salt, data = nc$map_loaded)
shapiro.test(nc$map_loaded$CH4_ug_m2_h)
t.test(CH4_ug_m2_h ~ Salt, data = nc$map_loaded) # Marg.
wilcox.test(CH4_ug_m2_h ~ Salt, data = nc$map_loaded) # Marg.

label_df2 <- data.frame("Exp" = c("Obs", "Obs",
                                  "Field", "Field",
                                  "Field", "Field",
                                  "Lab", "Lab", 
                                  "Lab", "Lab"),
                        "Estuary2" = c("SF", "SF",
                                       "Waccamaw", "Waccamaw",
                                       "Delaware", "Delaware",
                                       "Delaware", "Delaware",
                                       "Alligator", "Alligator"),
                        "y" = c(max(sf$map_loaded$CH4_ug_m2_h) + max(sf$map_loaded$CH4_ug_m2_h)/10,
                                max(sf$map_loaded$CH4_ug_m2_h) + max(sf$map_loaded$CH4_ug_m2_h)/10,
                                max(sc$map_loaded$CH4_ug_m2_h) + max(sc$map_loaded$CH4_ug_m2_h)/10,
                                max(sc$map_loaded$CH4_ug_m2_h) + max(sc$map_loaded$CH4_ug_m2_h)/10,
                                max(detra$map_loaded$CH4_ug_m2_h) + max(detra$map_loaded$CH4_ug_m2_h)/10,
                                max(detra$map_loaded$CH4_ug_m2_h) + max(detra$map_loaded$CH4_ug_m2_h)/10,
                                6.93 + 6.93/10,
                                6.93 + 6.93/10,
                                max(nc$map_loaded$CH4_ug_m2_h) + max(nc$map_loaded$CH4_ug_m2_h)/10,
                                max(nc$map_loaded$CH4_ug_m2_h) + max(nc$map_loaded$CH4_ug_m2_h)/10),
                        "x" = c("Freshwater", "Oligohaline", "Freshwater", "Oligohaline",
                                "Freshwater", "Oligohaline", "Freshwater", "Oligohaline",
                                "Freshwater", "Oligohaline"),
                        "label" = c("a", "b", "a", "b", "", "", "", "", "a", "b")) %>%
  mutate(Exp = factor(Exp, levels = c("Obs", "Field", "Lab"))) %>%
  mutate(Estuary2 = factor(Estuary2, levels = c("SF", "Waccamaw", "Delaware", "Alligator")))
           
# Only use top depth (data is repeated)
# Remove depth fill
#png("FinalFigs/Figure6.png", width = 8, height = 4, units = "in", res = 300)
ggplot(subset(frol$map_loaded, Depth == "0-5"), aes(Salt, CH4_ug_m2_h)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Salt)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, 
              aes(shape = Estuary2, colour = Salt)) +
  geom_text(data = label_df2, aes(x, y, label = label), 
            size = 4, color = "black", inherit.aes = F) +
  labs(y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)"),
       colour = "Salinity",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red")) +
  scale_shape_manual(breaks = c( "SF", "Waccamaw", "Delaware", "Alligator"), 
                     values = c(24, 21, 23, 22)) +
  guides(shape = guide_legend(order = 1,),
         colour = guide_legend(order = 2)) +
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
#dev.off()
  

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
sf_corMGMT <- cor.test(sf$map_loaded$MG_MT, sf$map_loaded$CH4_ug_m2_h, method = "spearman")
sf_res <- rbind(sf_cor, sf_corG, sf_corP) %>%
  add_row("Variable" = "MG:MT", "rho" = round(sf_corMGMT$estimate, digits = 2), "SpearmanPcut" = "Pfdr < 0.05") %>%
  rename(SF_rho = rho,
         SF_sig = SpearmanPcut)
sf_res$Variable <- replace(sf_res$Variable, sf_res$Variable == "sed_pH" , "pH")

# Note, DE is just field, because no lab methane data yet
detra_cor <- multiple_correlations(env_nona = env_nona_detra, var = "CH4_ug_m2_h") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
detra_corG <- meth_corr_by_taxonomy(input = detra, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
detra_corP <- meth_corr_by_taxonomy(input = detra, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
detra_corMGMT <- cor.test(detra$map_loaded$MG_MT, detra$map_loaded$CH4_ug_m2_h, method = "spearman")
de_res <- rbind(detra_cor, detra_corG, detra_corP) %>%
  add_row("Variable" = "MG:MT", "rho" = round(detra_corMGMT$estimate, digits = 2), "SpearmanPcut" = "Pfdr > 0.05") %>%
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
nc_corMGMT <- cor.test(nc$map_loaded$MG_MT, nc$map_loaded$CH4_ug_m2_h, method = "spearman")
nc_res <- rbind(nc_cor, nc_corG, nc_corP) %>%
  add_row("Variable" = "MG:MT", "rho" = round(nc_corMGMT$estimate, digits = 2), "SpearmanPcut" = "Pfdr > 0.05") %>%
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
sc_corMGMT <- cor.test(sc$map_loaded$MG_MT, sc$map_loaded$CH4_ug_m2_h, method = "spearman")
sc_res <- rbind(sc_cor, sc_corG, sc_corP) %>%
  add_row("Variable" = "MG:MT", "rho" = round(sc_corMGMT$estimate, digits = 2), "SpearmanPcut" = "Pfdr > 0.05") %>%
  rename(WA_rho = rho,
         WA_sig = SpearmanPcut)

# Make df
# Prep. data frame
# Update, need to merge WA SRP_umol_m2_h and DE SR_umol_cm3_d as SR column
CH4_res <- data.frame("Variable" = c("CH4_H2", "CH4_ac", "CH4_me", "CH4_mix", "MG:MT",  
                                     "ANME", "MOB_I", "MOB_II", "MOB_IIa",
                                     "Salinity_ppt_all", "NO3_mgL", "sed_NO3_mgL", 
                                     "N2O_ug_m2_h", "SRP_umol_m2_h", "SR_umol_cm3_d", "SO4_mgL", 
                                     "sed_SO4_mgL", "SRB", "Fe_mgL", "sed_Fe_mgL", 
                                     "FeRB", "Mn_mgL", "sed_Mn_mgL",  
                                     "CO2_ug_m2_h", "sed_per_org", "DOC_mgL", 
                                     "TOC_mgL", "sed_per_C", "sed_CN", 
                                     "sed_per_N", "DON_mgL", "DIN_mgL", "TN_mgL", 
                                     "NH4_mgL", "sed_NH4_mgL", "PO4_mgL", "sed_PO4_mgL",
                                     "pH", "Firmicutes", "Actinobacteriota"),
                      "Shortname" = c("CH4_H2", "CH4_ac", "CH4_me", "CH4_mix", "MG:MT",
                                      "ANME", "MOB_I", "MOB_II", "MOB_IIa",
                                      "Salinity", "NO3", "NO3_s", 
                                      "N2O Flux", "SR", "SR", "SO4", "SO4_s", "SRB",
                                      "Fe", "Fe_s", "FeRB", "Mn", "Mn_s", 
                                      "CO2 Flux", "SOM", "DOC", "TOC", "C", "C:N", 
                                      "N", "DON", "DIN", "TN", "NH4", "NH4_s", "PO4",
                                      "PO4_s", "pH", "Firmicutes", "Actinobacteriota"),
                      "Type" = c("Sediment", "Sediment", "Sediment", "Sediment", "Sediment",
                                 "Sediment", "Sediment", "Sediment", "Sediment",
                                 "Porewater", "Porewater", "Sediment", "Flux", 
                                 "Sediment", "Sediment", "Porewater", "Sediment", "Sediment",
                                 "Porewater", "Sediment", "Sediment", "Porewater", 
                                 "Sediment", 
                                 "Flux", "Sediment", "Porewater", "Porewater", 
                                 "Sediment", "Sediment", "Sediment", "Porewater", 
                                 "Porewater", "Porewater", "Porewater", "Sediment", 
                                 "Porewater", "Sediment", "Porewater", "Sediment", "Sediment"),
                      "Prediction" = c("Positive", "Positive", "Positive", "Positive", "Positive",
                                       "Negative", "Negative", "Negative", "Negative",
                                       "Negative", "Negative", "Negative", "Negative", 
                                       "Negative", "Negative", "Negative", "Negative", "Negative", 
                                       "Negative", "Negative", "Negative", "Negative",
                                       "Negative", 
                                       "Positive", "Positive", "Positive", "Positive",
                                       "Positive", "Negative", "Positive", "Positive", 
                                       "Positive", "Positive", "Positive", "Positive", 
                                       "Positive", "Positive", "Positive", "Positive", "Positive"),
                      "Variable Type" = c("Microbial", "Microbial", "Microbial", "Microbial", "Microbial",
                                          "Microbial", "Microbial", "Microbial", "Microbial",
                                          "Chemical", "Chemical", "Chemical", "Chemical",
                                     "Microbial", "Microbial", "Chemical", "Chemical", "Microbial",
                                     "Chemical", "Chemical", "Microbial", "Chemical",
                                     "Chemical",
                                     "Chemical", "Chemical", "Chemical", "Chemical",
                                     "Chemical", "Chemical", "Chemical", "Chemical",
                                     "Chemical", "Chemical", "Chemical", "Chemical",
                                     "Chemical", "Chemical", "Chemical", "Microbial", "Microbial"),
                      "Hypothesis" = c("Methanogens", "Methanogens", "Methanogens", "Methanogens", "Methanogens",
                                       "Methanotrophs", "Methanotrophs", "Methanotrophs", "Methanotrophs",
                                       "Alternate e-", "Alternate e-", "Alternate e-", "Alternate e-", 
                                       "Alternate e-", "Alternate e-", "Alternate e-", "Alternate e-", 
                                       "Alternate e-", "Alternate e-", "Alternate e-", "Alternate e-",
                                       "Alternate e-", "Alternate e-",
                                       "Decomposition", "Decomposition", "Decomposition", "Decomposition",
                                       "Decomposition", "Decomposition", "Decomposition", "Decomposition", 
                                       "Decomposition", "Decomposition", "Decomposition", "Decomposition", 
                                       "Decomposition", "Decomposition", "Decomposition", "Decomposition", 
                                       "Decomposition")) %>%
  left_join(., sf_res, by = "Variable") %>%
  left_join(., de_res, by = "Variable") %>%
  left_join(., nc_res, by = "Variable") %>%
  left_join(., sc_res, by = "Variable")
for (i in 1:nrow(CH4_res)) {
  if(CH4_res$Variable[i] == "SRP_umol_m2_h") {
         CH4_res$DE_rho[i] <- 0.10
  }
  if(CH4_res$Variable[i] == "SRP_umol_m2_h") {
    CH4_res$DE_sig[i] <- "Pfdr > 0.05"
  }
}
CH4_res <- subset(CH4_res, Variable != "SR_umol_cm3_d")

# Pretty heatmap
CH4_res_meta <- CH4_res %>%
  dplyr::select(Shortname, Type, Variable.Type, Hypothesis, Prediction, SF_sig, WA_sig, DE_sig, AL_sig) %>%
  mutate_if(is.character, as.factor)
CH4_res_mat <- CH4_res %>%
  rownames_to_column(var = "Num") %>%
  dplyr::select(-Num) %>%
  column_to_rownames(var = "Shortname") %>%
  dplyr::select(SF_rho, WA_rho, DE_rho, AL_rho) %>%
  rename("Wacc_rho" = WA_rho,
         "Alli_rho" = AL_rho) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(CH4_res_mat), 
                       "Alli_sig" = CH4_res_meta$`AL_sig`,
                       "DE_sig" = CH4_res_meta$`DE_sig`,
                       "Wacc_sig" = CH4_res_meta$`WA_sig`,
                       "SF_sig" = CH4_res_meta$`SF_sig`,
                       Prediction = CH4_res_meta$Prediction,
                       Hypothesis = CH4_res_meta$Hypothesis,
                       Variable = CH4_res_meta$Variable.Type,
                       Type = CH4_res_meta$Type)
ann_colors <- list(Type = c(Flux = "#FFFF99", 
                            Porewater = "#A6CEE3", 
                            Sediment = "#B15928"),
                   Variable = c(Chemical = "#440154FF",
                                Microbial = "#FDE725FF"),
                   Hypothesis = c(Methanogens = "#F8766D",
                                  Methanotrophs = "#7CAE00",
                                  `Alternate e-` = "#00BFC4",
                                  Decomposition = "#C77CFF"),
                   Prediction = c(Positive = "red", 
                                  Negative = "blue"),
                   "SF_sig" = c(`Pfdr < 0.05` = "black", 
                              `Pfdr > 0.05` = "white"),
                   "DE_sig" = c(`Pfdr < 0.05` = "black", 
                              `Pfdr > 0.05` = "white"),
                   "Alli_sig" = c(`Pfdr < 0.05` = "black", 
                              `Pfdr > 0.05` = "white"),
                   "Wacc_sig" = c(`Pfdr < 0.05` = "black", 
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
         fontsize = 6,
         fontsize_row = 9,
         fontsize_col = 9,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = c(5, 9, 22),
         filename = "FinalFigs/Figure5.png",
         width = 5,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### _Sal ####
# Correlations between salinity and BGC and Guilds by each site
# Use the custom functions written to do this

# Get all corrs
sf_cor <- multiple_correlations(env_nona = env_nona_sf, var = "Salinity_ppt_all") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sf_corG <- meth_corr_by_taxonomy(input = sf, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sf_corP <- meth_corr_by_taxonomy(input = sf, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sf_corMGMT <- cor.test(sf$map_loaded$MG_MT, sf$map_loaded$Salinity_ppt_all, method = "spearman")
sf_res <- rbind(sf_cor, sf_corG, sf_corP) %>%
  add_row("Variable" = "MG:MT", "rho" = round(sf_corMGMT$estimate, digits = 2), "SpearmanPcut" = "Pfdr > 0.05") %>%
  rename(SF_rho = rho,
         SF_sig = SpearmanPcut)
sf_res$Variable <- replace(sf_res$Variable, sf_res$Variable == "sed_pH" , "pH")

# Note, DE is just field, because no lab methane data yet
detra_cor <- multiple_correlations(env_nona = env_nona_detra, var = "Salinity_ppt_all") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
detra_corG <- meth_corr_by_taxonomy(input = detra, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
detra_corP <- meth_corr_by_taxonomy(input = detra, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
detra_corMGMT <- cor.test(detra$map_loaded$MG_MT, detra$map_loaded$Salinity_ppt_all, method = "spearman")
de_res <- rbind(detra_cor, detra_corG, detra_corP) %>%
  add_row("Variable" = "MG:MT", "rho" = round(detra_corMGMT$estimate, digits = 2), "SpearmanPcut" = "Pfdr > 0.05") %>%
  rename(DE_rho = rho,
         DE_sig = SpearmanPcut)

nc_cor <- multiple_correlations(env_nona = env_nona_nc, var = "Salinity_ppt_all") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
nc_corG <- meth_corr_by_taxonomy(input = nc, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
nc_corP <- meth_corr_by_taxonomy(input = nc, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
nc_corMGMT <- cor.test(nc$map_loaded$MG_MT, nc$map_loaded$Salinity_ppt_all, method = "spearman")
nc_res <- rbind(nc_cor, nc_corG, nc_corP) %>%
  add_row("Variable" = "MG:MT", "rho" = round(nc_corMGMT$estimate, digits = 2), "SpearmanPcut" = "Pfdr < 0.05") %>%
  rename(AL_rho = rho,
         AL_sig = SpearmanPcut)

sc_cor <- multiple_correlations(env_nona = env_nona_sc, var = "Salinity_ppt_all") %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sc_corG <- meth_corr_by_taxonomy(input = sc, level = 9, threshold = 0, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sc_corP <- meth_corr_by_taxonomy(input = sc, level = 2, threshold = 0.5, data = "Yes") %>%
  rename(Variable = Taxon) %>%
  dplyr::select(Variable, rho, SpearmanPcut)
sc_corMGMT <- cor.test(sc$map_loaded$MG_MT, sc$map_loaded$Salinity_ppt_all, method = "spearman")
sc_res <- rbind(sc_cor, sc_corG, sc_corP) %>%
  add_row("Variable" = "MG:MT", "rho" = round(sc_corMGMT$estimate, digits = 2), "SpearmanPcut" = "Pfdr > 0.05") %>%
  rename(WA_rho = rho,
         WA_sig = SpearmanPcut)

# Make df
# Prep. data frame
# Update, need to merge WA SRP_umol_m2_h and DE SR_umol_cm3_d as SR column
Sal_res <- data.frame("Variable" = c("CH4_H2", "CH4_ac", "CH4_me", "CH4_mix", "MG:MT",  
                                     "ANME", "MOB_I", "MOB_II", "MOB_IIa",
                                     "Salinity_ppt_all", "NO3_mgL", "sed_NO3_mgL", 
                                     "N2O_ug_m2_h", "SRP_umol_m2_h", "SR_umol_cm3_d", "SO4_mgL", 
                                     "sed_SO4_mgL", "SRB", "Fe_mgL", "sed_Fe_mgL", 
                                     "FeRB", "Mn_mgL", "sed_Mn_mgL",  
                                     "CO2_ug_m2_h", "sed_per_org", "DOC_mgL", 
                                     "TOC_mgL", "sed_per_C", "sed_CN", 
                                     "sed_per_N", "DON_mgL", "DIN_mgL", "TN_mgL", 
                                     "NH4_mgL", "sed_NH4_mgL", "PO4_mgL", "sed_PO4_mgL",
                                     "pH", "Firmicutes", "Actinobacteriota"),
                      "Shortname" = c("CH4_H2", "CH4_ac", "CH4_me", "CH4_mix", "MG:MT",
                                      "ANME", "MOB_I", "MOB_II", "MOB_IIa",
                                      "Salinity", "NO3", "NO3_s", 
                                      "N2O Flux", "SR", "SR", "SO4", "SO4_s", "SRB",
                                      "Fe", "Fe_s", "FeRB", "Mn", "Mn_s", 
                                      "CO2 Flux", "SOM", "DOC", "TOC", "C", "C:N", 
                                      "N", "DON", "DIN", "TN", "NH4", "NH4_s", "PO4",
                                      "PO4_s", "pH", "Firmicutes", "Actinobacteriota"),
                      "Type" = c("Sediment", "Sediment", "Sediment", "Sediment", "Sediment",
                                 "Sediment", "Sediment", "Sediment", "Sediment",
                                 "Porewater", "Porewater", "Sediment", "Flux", 
                                 "Sediment", "Sediment", "Porewater", "Sediment", "Sediment",
                                 "Porewater", "Sediment", "Sediment", "Porewater", 
                                 "Sediment", 
                                 "Flux", "Sediment", "Porewater", "Porewater", 
                                 "Sediment", "Sediment", "Sediment", "Porewater", 
                                 "Porewater", "Porewater", "Porewater", "Sediment", 
                                 "Porewater", "Sediment", "Porewater", "Sediment", "Sediment"),
                      "Prediction" = c("Positive", "Positive", "Positive", "Positive", "Positive",
                                       "Negative", "Negative", "Negative", "Negative",
                                       "Negative", "Negative", "Negative", "Negative", 
                                       "Negative", "Negative", "Negative", "Negative", "Negative", 
                                       "Negative", "Negative", "Negative", "Negative",
                                       "Negative", 
                                       "Positive", "Positive", "Positive", "Positive",
                                       "Positive", "Negative", "Positive", "Positive", 
                                       "Positive", "Positive", "Positive", "Positive", 
                                       "Positive", "Positive", "Positive", "Positive", "Positive"),
                      "Variable Type" = c("Microbial", "Microbial", "Microbial", "Microbial", "Microbial",
                                          "Microbial", "Microbial", "Microbial", "Microbial",
                                          "Chemical", "Chemical", "Chemical", "Chemical",
                                          "Microbial", "Microbial", "Chemical", "Chemical", "Microbial",
                                          "Chemical", "Chemical", "Microbial", "Chemical",
                                          "Chemical",
                                          "Chemical", "Chemical", "Chemical", "Chemical",
                                          "Chemical", "Chemical", "Chemical", "Chemical",
                                          "Chemical", "Chemical", "Chemical", "Chemical",
                                          "Chemical", "Chemical", "Chemical", "Microbial", "Microbial"),
                      "Hypothesis" = c("Methanogens", "Methanogens", "Methanogens", "Methanogens", "Methanogens",
                                       "Methanotrophs", "Methanotrophs", "Methanotrophs", "Methanotrophs",
                                       "Alternate e-", "Alternate e-", "Alternate e-", "Alternate e-", 
                                       "Alternate e-", "Alternate e-", "Alternate e-", "Alternate e-", 
                                       "Alternate e-", "Alternate e-", "Alternate e-", "Alternate e-",
                                       "Alternate e-", "Alternate e-",
                                       "Decomposition", "Decomposition", "Decomposition", "Decomposition",
                                       "Decomposition", "Decomposition", "Decomposition", "Decomposition", 
                                       "Decomposition", "Decomposition", "Decomposition", "Decomposition", 
                                       "Decomposition", "Decomposition", "Decomposition", "Decomposition", 
                                       "Decomposition")) %>%
  left_join(., sf_res, by = "Variable") %>%
  left_join(., de_res, by = "Variable") %>%
  left_join(., nc_res, by = "Variable") %>%
  left_join(., sc_res, by = "Variable")
for (i in 1:nrow(Sal_res)) {
  if(Sal_res$Variable[i] == "SRP_umol_m2_h") {
    Sal_res$DE_rho[i] <- 0.79
  }
  if(Sal_res$Variable[i] == "SRP_umol_m2_h") {
    Sal_res$DE_sig[i] <- "Pfdr > 0.05"
  }
}
Sal_res <- subset(Sal_res, Variable != "SR_umol_cm3_d")

# Pretty heatmap
Sal_res_meta <- Sal_res %>%
  dplyr::select(Shortname, Type, Variable.Type, Hypothesis, Prediction, SF_sig, WA_sig, DE_sig, AL_sig) %>%
  mutate_if(is.character, as.factor)
Sal_res_mat <- Sal_res %>%
  rownames_to_column(var = "Num") %>%
  dplyr::select(-Num) %>%
  column_to_rownames(var = "Shortname") %>%
  dplyr::select(SF_rho, WA_rho, DE_rho, AL_rho) %>%
  rename("Wacc_rho" = WA_rho,
         "Alli_rho" = AL_rho) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(Sal_res_mat), 
                       "Alli_sig" = Sal_res_meta$`AL_sig`,
                       "DE_sig" = Sal_res_meta$`DE_sig`,
                       "Wacc_sig" = Sal_res_meta$`WA_sig`,
                       "SF_sig" = Sal_res_meta$`SF_sig`,
                       Prediction = Sal_res_meta$Prediction,
                       Hypothesis = Sal_res_meta$Hypothesis,
                       Variable = Sal_res_meta$Variable.Type,
                       Type = Sal_res_meta$Type)
ann_colors <- list(Type = c(Flux = "#FFFF99", 
                            Porewater = "#A6CEE3", 
                            Sediment = "#B15928"),
                   Variable = c(Chemical = "#440154FF",
                                Microbial = "#FDE725FF"),
                   Hypothesis = c(Methanogens = "#F8766D",
                                  Methanotrophs = "#7CAE00",
                                  `Alternate e-` = "#00BFC4",
                                  Decomposition = "#C77CFF"),
                   Prediction = c(Positive = "red", 
                                  Negative = "blue"),
                   "SF_sig" = c(`Pfdr < 0.05` = "black", 
                                `Pfdr > 0.05` = "white"),
                   "DE_sig" = c(`Pfdr < 0.05` = "black", 
                                `Pfdr > 0.05` = "white"),
                   "Alli_sig" = c(`Pfdr < 0.05` = "black", 
                                `Pfdr > 0.05` = "white"),
                   "Wacc_sig" = c(`Pfdr < 0.05` = "black", 
                                `Pfdr > 0.05` = "white"))
pheatmap(Sal_res_mat,
         legend = T,
         legend_breaks = c(-0.5, 0, 0.5, max(na.omit(Sal_res_mat))),
         legend_labels = c("-0.5", "0", "0.5", "rho\n"),
         main = "",
         #color = bluered(100),
         #border_color = NA,
         #na_col = "gray",
         scale = "none",
         angle_col = 315,
         fontsize = 6,
         fontsize_row = 9,
         fontsize_col = 9,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = c(5, 9, 21),
         filename = "FinalFigs/FigureS2.png",
         width = 5,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



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
  scale_y_continuous(expand = c(max(tallest_barG$sum)/100, max(tallest_barG$sum)/100)) + 
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
  scale_y_continuous(expand = c(max(tallest_barM$sum)/100, max(tallest_barM$sum)/100)) + 
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
#png("Figure7_3panel.png", width = 5, height = 7, units = "in", res = 300)
plot_grid(p, l, ncol = 2, align = "h", rel_widths = c(0.57, 0.43))
#dev.off()

# Updated Figure - just show Methanobacteriaceae, and for all samples
# Hypothesis is specifically about OTUs in same family/genus having different salt responses
lab_barsMethanobacteriaceae <- plot_taxa_bars(Methanobacteriaceae_wTax, 
                                              lab$map_loaded, 
                                              "sampleID", 
                                              num_taxa = nrow(Methanobacteriaceae_wTax), 
                                              data_only = T) %>%
  left_join(., lab$map_loaded, by = c("group_by" = "sampleID"), keep = TRUE) %>%
  mutate(mean_value = mean_value/26429*100)
tallest_barM <- lab_barsMethanobacteriaceae %>%
  group_by(sampleID) %>%
  summarise(sum = sum(mean_value))
png("FinalFigs/Figure7.png", width = 7, height = 5, units = "in", res = 300)
ggplot(lab_barsMethanobacteriaceae, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Methanobacteriaceae\nGenus; Species; OTU ID") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_barM$sum)/1000, max(tallest_barM$sum)/1000)) + 
  facet_nested(~ Estuary + Salt, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.line.y = element_blank(),
        strip.background = element_rect(size = 0.2),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10))
dev.off()



#### 4. Ecology ####
# Need to delve into some ecological dynamics to explain discrepancies between lab and field
# Make Venn, Jaccard vs. Bray, NTI

#### _Venn ####
source("~/Documents/GitHub/EastCoast/plot_venn_diagram2.R")

plot_venn_diagram2_mirror(sf, "Salt", 0.000000000000000000000000000001)
plot_venn_diagram2(sc, "Salt", 0.000000000000000000000000000001)
plot_venn_diagram2(detra, "Salt", 0.000000000000000000000000000001)
plot_venn_diagram2_mirror(deinc, "Salt", 0.000000000000000000000000000001)
plot_venn_diagram2(nc, "Salt", 0.000000000000000000000000000001)

png("FinalFigs/FigureS6.png", width = 9, height = 6, units = "in", res = 300)
plot_grid(plot_venn_diagram2_mirror(sf, "Salt", 0.000000000000000000000000000001),
          plot_venn_diagram2(sc, "Salt", 0.000000000000000000000000000001),
          plot_venn_diagram2(detra, "Salt", 0.000000000000000000000000000001),
          plot_venn_diagram2_mirror(deinc, "Salt", 0.000000000000000000000000000001),
          plot_venn_diagram2(nc, "Salt", 0.000000000000000000000000000001),
          labels = c("(a) SF (Field obs.)", "(b) Waccamaw (Field exp.)", "(c) Delaware (Field exp.)", 
                     "(d) Delaware (Lab exp.)", "(e) Alligator (Lab exp.)"),
          label_size = 10,
          label_x = -0.15,
          label_y = 0.9)
dev.off()

# Also get site Venn at all levels - Figure S11
# Need short site names
frol$map_loaded$EstuaryShort <- recode_factor(frol$map_loaded$Estuary2,
                                            "SF" = "SF",
                                            "Waccamaw" = "Wacc.",
                                            "Delaware" = "DE",
                                            "Alligator" = "Alli.")
frol$map_loaded$EstuaryShort <- factor(frol$map_loaded$EstuaryShort,
                                       levels = c("SF", "Wacc.", "DE", "Alli."))
phy <- summarize_taxonomy(frol, level = 2, report_higher_tax = F)
cla <- summarize_taxonomy(frol, level = 3, report_higher_tax = F)
ord <- summarize_taxonomy(frol, level = 4, report_higher_tax = F)
fam <- summarize_taxonomy(frol, level = 5, report_higher_tax = F)
gen <- summarize_taxonomy(frol, level = 6, report_higher_tax = F)

input_phylum <- frol
input_phylum$data_loaded <- phy
input_class <- frol
input_class$data_loaded <- cla
input_order <- frol
input_order$data_loaded <- ord
input_family <- frol
input_family$data_loaded <- fam
input_genus <- frol
input_genus$data_loaded <- gen

png("FinalFigs/FigureS11.png", width = 9, height = 6, units = "in", res = 300)
plot_grid(plot_venn_diagram(input_phylum, "EstuaryShort", 0.00000000000000001),
          plot_venn_diagram(input_class, "EstuaryShort", 0.00000000000000001),
          plot_venn_diagram(input_order, "EstuaryShort", 0.00000000000000001),
          plot_venn_diagram(input_family, "EstuaryShort", 0.00000000000000001),
          plot_venn_diagram(input_genus, "EstuaryShort", 0.00000000000000001),
          plot_venn_diagram(frol, "EstuaryShort", 0.00000000000000001),
          labels = c("(a) Phylum", "(b) Class", "(c) Order", 
                     "(d) Family", "(e) Genus", "(f) OTU"))
dev.off()



#### _Jaccard/Bray ####
# Bray
bray.distance <- calc_dm(frol$data_loaded, method = "bray_sq_trans")
bac_bray_mat <- as.matrix(bray.distance)
bac_bray_mat[upper.tri(bac_bray_mat, diag = TRUE)] <- NA
bac_bray_df <- as.data.frame(bac_bray_mat)
bac_bray_df$sampleID <- rownames(bac_bray_df)
bac_bray_df_long <- melt(bac_bray_df, id.vars = "sampleID")
bac_bray_df_long <- na.omit(bac_bray_df_long)
nrow(bac_bray_df_long) == (nrow(frol$map_loaded)*(nrow(frol$map_loaded)-1))/2
bac_bray_df_long$sampleID <- as.factor(bac_bray_df_long$sampleID)
site_sal <- dplyr::select(frol$map_loaded, sampleID, Estuary3, Salt, Estuary2, Exp)
bac_bray_df_long <- inner_join(bac_bray_df_long, site_sal, 
                               by = c("sampleID" = "sampleID"))
bac_bray_df_long <- inner_join(bac_bray_df_long, site_sal, 
                               by = c("variable" = "sampleID"))
bac_bray_df_long <- subset(bac_bray_df_long, Estuary3.x == Estuary3.y)
for (i in 1:nrow(bac_bray_df_long)) {
  ifelse(bac_bray_df_long$Salt.x[i] == bac_bray_df_long$Salt.y[i],
         bac_bray_df_long$comparison[i] <- "within",
         bac_bray_df_long$comparison[i] <- "between")
}
bac_bray_df_long$comparison <- as.factor(bac_bray_df_long$comparison)
bac_bray_df_long_sf <- subset(bac_bray_df_long, Estuary3.x == "SF")
bac_bray_df_long_wa <- subset(bac_bray_df_long, Estuary3.x == "Waccamaw")
bac_bray_df_long_def <- subset(bac_bray_df_long, Estuary3.x == "DE Field")
bac_bray_df_long_del <- subset(bac_bray_df_long, Estuary3.x == "DE Lab")
bac_bray_df_long_al <- subset(bac_bray_df_long, Estuary3.x == "Alligator")
table(bac_bray_df_long_sf$comparison)
table(bac_bray_df_long_wa$comparison)
table(bac_bray_df_long_def$comparison)
table(bac_bray_df_long_del$comparison)
table(bac_bray_df_long_al$comparison)
t.test(value ~ comparison, data = bac_bray_df_long_sf) # Sig
t.test(value ~ comparison, data = bac_bray_df_long_wa) # Sig
t.test(value ~ comparison, data = bac_bray_df_long_def) # NSD
t.test(value ~ comparison, data = bac_bray_df_long_del) # NSD
t.test(value ~ comparison, data = bac_bray_df_long_al) # Sig
cohensD(value ~ comparison, data = bac_bray_df_long_sf) # 1.23
cohensD(value ~ comparison, data = bac_bray_df_long_wa) # 0.97
cohensD(value ~ comparison, data = bac_bray_df_long_def) # 0.14
cohensD(value ~ comparison, data = bac_bray_df_long_del) # 0.03
cohensD(value ~ comparison, data = bac_bray_df_long_al) # 1.02
bac_bray_df_long_multi <- rbind(bac_bray_df_long_sf,
                                bac_bray_df_long_wa,
                                bac_bray_df_long_def,
                                bac_bray_df_long_del,
                                bac_bray_df_long_al)

# Jaccard
jac.distance <- calc_dm(frol$data_loaded, method = "jaccard")
bac_jac_mat <- as.matrix(jac.distance)
bac_jac_mat[upper.tri(bac_jac_mat, diag = TRUE)] <- NA
bac_jac_df <- as.data.frame(bac_jac_mat)
bac_jac_df$sampleID <- rownames(bac_jac_df)
bac_jac_df_long <- melt(bac_jac_df, id.vars = "sampleID")
bac_jac_df_long <- na.omit(bac_jac_df_long)
nrow(bac_jac_df_long) == (nrow(frol$map_loaded)*(nrow(frol$map_loaded)-1))/2
bac_jac_df_long$sampleID <- as.factor(bac_jac_df_long$sampleID)
site_sal <- dplyr::select(frol$map_loaded, sampleID, Estuary3, Salt, Estuary2, Exp)
bac_jac_df_long <- inner_join(bac_jac_df_long, site_sal, 
                               by = c("sampleID" = "sampleID"))
bac_jac_df_long <- inner_join(bac_jac_df_long, site_sal, 
                               by = c("variable" = "sampleID"))
bac_jac_df_long <- subset(bac_jac_df_long, Estuary3.x == Estuary3.y)
for (i in 1:nrow(bac_jac_df_long)) {
  ifelse(bac_jac_df_long$Salt.x[i] == bac_jac_df_long$Salt.y[i],
         bac_jac_df_long$comparison[i] <- "within",
         bac_jac_df_long$comparison[i] <- "between")
}
bac_jac_df_long$comparison <- as.factor(bac_jac_df_long$comparison)
bac_jac_df_long_sf <- subset(bac_jac_df_long, Estuary3.x == "SF")
bac_jac_df_long_wa <- subset(bac_jac_df_long, Estuary3.x == "Waccamaw")
bac_jac_df_long_def <- subset(bac_jac_df_long, Estuary3.x == "DE Field")
bac_jac_df_long_del <- subset(bac_jac_df_long, Estuary3.x == "DE Lab")
bac_jac_df_long_al <- subset(bac_jac_df_long, Estuary3.x == "Alligator")
table(bac_jac_df_long_sf$comparison)
table(bac_jac_df_long_wa$comparison)
table(bac_jac_df_long_def$comparison)
table(bac_jac_df_long_del$comparison)
table(bac_jac_df_long_al$comparison)
t.test(value ~ comparison, data = bac_jac_df_long_sf) # Sig
t.test(value ~ comparison, data = bac_jac_df_long_wa) # Sig
t.test(value ~ comparison, data = bac_jac_df_long_def) # NSD
t.test(value ~ comparison, data = bac_jac_df_long_del) # NSD
t.test(value ~ comparison, data = bac_jac_df_long_al) # Sig
cohensD(value ~ comparison, data = bac_jac_df_long_sf) # 1.31
cohensD(value ~ comparison, data = bac_jac_df_long_wa) # 1.13
cohensD(value ~ comparison, data = bac_jac_df_long_def) # 0.25
cohensD(value ~ comparison, data = bac_jac_df_long_del) # 0.04
cohensD(value ~ comparison, data = bac_jac_df_long_al) # 1.15
bac_jac_df_long_multi <- rbind(bac_jac_df_long_sf,
                                bac_jac_df_long_wa,
                                bac_jac_df_long_def,
                                bac_jac_df_long_del,
                                bac_jac_df_long_al)

# Graph
bac_bray_df_long_multi$Dissim <- "Bray-Curtis"
bac_jac_df_long_multi$Dissim <- "Jaccard"
comb <- rbind(bac_bray_df_long_multi, bac_jac_df_long_multi) %>%
  mutate(Exp.x = factor(Exp.x, levels = c("Obs", "Field", "Lab")),
         Estuary2.x = factor(Estuary2.x, levels = c("SF", "Waccamaw", "Delaware", "Alligator")))
label.df <- data.frame(Dissim = c("Bray-Curtis", "Bray-Curtis", 
                                  "Bray-Curtis", "Bray-Curtis", 
                                  "Bray-Curtis", "Bray-Curtis",
                                  "Bray-Curtis", "Bray-Curtis",
                                  "Bray-Curtis", "Bray-Curtis",
                                  "Jaccard", "Jaccard", 
                                  "Jaccard", "Jaccard", 
                                  "Jaccard", "Jaccard",
                                  "Jaccard", "Jaccard",
                                  "Jaccard", "Jaccard"),
                       Exp.x = c("Obs", "Obs", "Field", "Field", "Field", "Field",
                               "Lab", "Lab", "Lab", "Lab",
                               "Obs", "Obs", "Field", "Field", "Field", "Field",
                               "Lab", "Lab", "Lab", "Lab"),
                       Estuary2.x = c("SF", "SF",
                                    "Waccamaw", "Waccamaw",
                                    "Delaware", "Delaware",
                                    "Delaware", "Delaware",
                                    "Alligator", "Alligator",
                                    "SF", "SF",
                                    "Waccamaw", "Waccamaw",
                                    "Delaware", "Delaware",
                                    "Delaware", "Delaware",
                                    "Alligator", "Alligator"),
                       comparison = c("between","within","between","within",
                                      "between","within","between","within",
                                      "between","within","between","within",
                                      "between","within","between","within",
                                      "between","within","between","within"),
                       Value = c(1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,
                                 1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05),
                       Sig = c("a","b","a","b","","","","","a","b",
                               "a","b","a","b","","","","","a","b")) %>%
  mutate(Exp.x = factor(Exp.x, levels = c("Obs", "Field", "Lab")),
         Estuary2.x = factor(Estuary2.x, levels = c("SF", "Waccamaw", "Delaware", "Alligator")))
label.df2 <- data.frame(Dissim = c("Bray-Curtis", "Bray-Curtis", 
                                   "Bray-Curtis", "Bray-Curtis", 
                                   "Bray-Curtis", "Jaccard",
                                   "Jaccard", "Jaccard",
                                   "Jaccard", "Jaccard"),
                       Exp.x = c("Obs", "Field", "Field", "Lab", "Lab", 
                                 "Obs", "Field", "Field", "Lab", "Lab"),
                       Estuary2.x = c("SF",
                                      "Waccamaw",
                                      "Delaware",
                                      "Delaware",
                                      "Alligator",
                                      "SF",
                                      "Waccamaw",
                                      "Delaware",
                                      "Delaware",
                                      "Alligator"),
                       x = c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5),
                       Value = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                       es = c("d = 1.23","d = 0.97","d = 0.14","d = 0.03","d = 1.02",
                              "d = 1.31","d = 1.13","d = 0.25","d = 0.04","d = 1.23")) %>%
  mutate(Exp.x = factor(Exp.x, levels = c("Obs", "Field", "Lab")),
         Estuary2.x = factor(Estuary2.x, levels = c("SF", "Waccamaw", "Delaware", "Alligator")))
#png("FinalFigs/FigureS8.png", width = 8, height = 5, units = "in", res = 300)
ggplot(data = comb, aes(comparison, value)) +
  geom_jitter(data = subset(comb, Estuary2.x == "SF"),
              size = 1, alpha = 0.2) +
  geom_jitter(data = subset(comb, Estuary2.x == "Waccamaw"),
              size = 1, alpha = 0.3) +
  geom_jitter(data = subset(comb, Estuary2.x == "Delaware"), 
              size = 1, alpha = 0.5) +
  geom_jitter(data = subset(comb, Estuary2.x == "Alligator"), 
              size = 1, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, color = "blue", fill = NA) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group=NULL), size = 4) +
  geom_text(data = label.df2, aes(x = x, y = Value, label = es, group=NULL), size = 3.5) +
  labs(x = "Salinity Comparison",
       y = "Dissimilarity") +
  facet_nested(Dissim ~ Exp.x + Estuary2.x, scales = "free_y") +
  ylim(0,1.05) +
  theme_bw() +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 11))
#dev.off()



#### _NTI ####
# The NTI values between −2 and 2 indicate stochastic community assembly, whereas NTI values less than −2 or higher than 2 indicate that deterministic processes play a more important role in structuring the community (deterministic environmental filtering)

# Align seqs with MUSCLE and Make tree with QIIME
# Remove rare taxa or tree is too big
frol_abund <- filter_taxa_from_input(frol, filter_thresh = 0.5)
bac_asv <- as.data.frame(t(frol_abund$data_loaded))
tree <- read.tree("rep_phylo.tre")
colnames(bac_asv)
tree$tip.label
tree <- prune.sample(bac_asv, tree)
bac_asv <- bac_asv[,tree$tip.label]
phy.dist <- cophenetic(tree)
NTI <- NTI.p(bac_asv, 
             phy.dist, 
             nworker = 4, 
             memo.size.GB = 50,
             weighted = TRUE, 
             rand = 1000,
             check.name = TRUE, 
             output.MNTD = FALSE,
             sig.index = "NTI",
             silent = FALSE)
saveRDS(NTI, "NTI.rds")

# But, what we really should do is NTI for each experiment.
# (Do this part on server for speed)
sf_asv <- as.data.frame(t(sf$data_loaded))
sf_tree <- read.tree("rep_phylo.tre")
sf_tree <- prune.sample(sf_asv, tree)
sf_asv <- sf_asv[,sf_tree$tip.label]
sf_phy.dist <- cophenetic(sf_tree)
sf_NTI <- NTI.p(sf_asv, 
             sf_phy.dist, 
             nworker = 4, 
             memo.size.GB = 50,
             weighted = TRUE, 
             rand = 1000,
             check.name = TRUE, 
             output.MNTD = FALSE,
             sig.index = "NTI",
             silent = FALSE)

sc_asv <- as.data.frame(t(sc$data_loaded))
sc_tree <- read.tree("rep_phylo.tre")
sc_tree <- prune.sample(sc_asv, tree)
sc_asv <- sc_asv[,sc_tree$tip.label]
sc_phy.dist <- cophenetic(sc_tree)
sc_NTI <- NTI.p(sc_asv, 
                sc_phy.dist, 
                nworker = 4, 
                memo.size.GB = 50,
                weighted = TRUE, 
                rand = 1000,
                check.name = TRUE, 
                output.MNTD = FALSE,
                sig.index = "NTI",
                silent = FALSE)

detra_asv <- as.data.frame(t(detra$data_loaded))
detra_tree <- read.tree("rep_phylo.tre")
detra_tree <- prune.sample(detra_asv, tree)
detra_asv <- detra_asv[,detra_tree$tip.label]
detra_phy.dist <- cophenetic(detra_tree)
detra_NTI <- NTI.p(detra_asv, 
                detra_phy.dist, 
                nworker = 4, 
                memo.size.GB = 50,
                weighted = TRUE, 
                rand = 1000,
                check.name = TRUE, 
                output.MNTD = FALSE,
                sig.index = "NTI",
                silent = FALSE)

deinc_asv <- as.data.frame(t(deinc$data_loaded))
deinc_tree <- read.tree("rep_phylo.tre")
deinc_tree <- prune.sample(deinc_asv, tree)
deinc_asv <- deinc_asv[,deinc_tree$tip.label]
deinc_phy.dist <- cophenetic(deinc_tree)
deinc_NTI <- NTI.p(deinc_asv, 
                deinc_phy.dist, 
                nworker = 4, 
                memo.size.GB = 50,
                weighted = TRUE, 
                rand = 1000,
                check.name = TRUE, 
                output.MNTD = FALSE,
                sig.index = "NTI",
                silent = FALSE)

nc_asv <- as.data.frame(t(nc$data_loaded))
nc_tree <- read.tree("rep_phylo.tre")
nc_tree <- prune.sample(nc_asv, tree)
nc_asv <- nc_asv[,nc_tree$tip.label]
nc_phy.dist <- cophenetic(nc_tree)
nc_NTI <- NTI.p(nc_asv, 
                nc_phy.dist, 
                nworker = 4, 
                memo.size.GB = 50,
                weighted = TRUE, 
                rand = 1000,
                check.name = TRUE, 
                output.MNTD = FALSE,
                sig.index = "NTI",
                silent = FALSE)

NTI <- readRDS("NTI_concat.rds") %>%
  mutate(sampleID = rownames(.))

# Merge, test, and plot
frol$map_loaded <- frol$map_loaded %>%
  left_join(., NTI, by = "sampleID")

m2 <- aov(NTI ~ Estuary + Salt + Depth, data = frol$map_loaded)
summary(m2)
Anova(m2, type = "II")

m2 <- aov(NTI ~ EstSalt, data = frol$map_loaded)
shapiro.test(m2$residuals) # Normal
summary(m2)
t2 <- emmeans(object = m2, specs = "EstSalt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(frol$map_loaded$NTI)+(max(frol$map_loaded$NTI)-min(frol$map_loaded$NTI))/20)

t.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "SF")) # Sig.
wilcox.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "SF")) # Sig.

t.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "Waccamaw")) # Sig.
wilcox.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "Waccamaw")) # Sig.

t.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "DE Field")) # NSD
wilcox.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "DE Field")) # NSD

t.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "DE Lab")) # NSD
wilcox.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "DE Lab")) # NSD

t.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "Alligator")) # Sig.
wilcox.test(NTI ~ Salt, data = subset(frol$map_loaded, Estuary3 == "Alligator")) # Marg.

label_df3 <- data.frame("Exp" = c("Obs", "Obs",
                                  "Field", "Field",
                                  "Field", "Field",
                                  "Lab", "Lab", 
                                  "Lab", "Lab"),
                        "Estuary2" = c("SF", "SF",
                                       "Waccamaw", "Waccamaw",
                                       "Delaware", "Delaware",
                                       "Delaware", "Delaware",
                                       "Alligator", "Alligator"),
                        "y" = c(10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1, 10.1),
                        "x" = c("Freshwater", "Oligohaline", "Freshwater", "Oligohaline",
                                "Freshwater", "Oligohaline", "Freshwater", "Oligohaline",
                                "Freshwater", "Oligohaline"),
                        "label" = c("a", "b", "a", "b", "", "", "", "", "a", "b")) %>%
  mutate(Exp = factor(Exp, levels = c("Obs", "Field", "Lab"))) %>%
  mutate(Estuary2 = factor(Estuary2, levels = c("SF", "Waccamaw", "Delaware", "Alligator")))

#png("FinalFigs/FigureS7.png", width = 8, height = 4, units = "in", res = 300)
ggplot(frol$map_loaded, aes(Salt, NTI)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, aes(colour = Salt)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, 
              aes(fill = Depth, shape = Estuary2, colour = Salt)) +
  geom_text(data = label_df3, aes(x, y, label = label), 
            size = 4, color = "black", inherit.aes = F) +
  labs(y = "NTI",
       colour = "Salinity",
       shape = "Site",
       fill = "Depth (cm)") +
  scale_colour_manual(values = c("blue", "red")) +
  scale_shape_manual(breaks = c( "SF", "Waccamaw", "Delaware", "Alligator"), 
                     values = c(24, 21, 23, 22)) +
  scale_fill_manual(values = c("black", "white")) +
  guides(shape = guide_legend(order = 1,),
         colour = guide_legend(order = 2),
         fill = guide_legend(override.aes = list(shape = c(16, 1)), order = 3)) +
  facet_nested_wrap(~ Exp + Estuary2, nrow = 1) +
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
#dev.off()



#### 5. Map ####
# Need map as Figure S1
# Just basic US map with points on the 4 sites
coords <- data.frame(Estuary2 = c("SF", "Waccamaw", "Delaware", "Alligator"),
                     Longitude = c(-121.624, -79.0919, -75.17314722, -76.1569),
                     Latitude = c(38.001, 33.5250, 39.85918056, 35.9061)) %>%
  mutate(Estuary2 = factor(Estuary2,
                           levels = c("SF", "Waccamaw", "Delaware", "Alligator")))
test_data <- data.frame(lon = coords$Longitude, lat = coords$Latitude)
transformed_data <- usmap_transform(test_data)
coords$x <- transformed_data$x
coords$y <- transformed_data$y

#png("FinalFigs/FigureS1.png", width = 6, height = 4, units = "in", res = 300)
plot_usmap(exclude = c("AK", "HI"),
           color = "white",
           fill = "grey80",
           size = 0.3) +
  geom_point(data = coords, 
             aes(x = x, y = y, shape = Estuary2),
             fill = "red",
             color = "red",
             size = 4) +
  geom_text(data = coords,
            aes(x = x, y = y + 130000, label = Estuary2),
            size = 4) +
  scale_shape_manual(breaks = c( "SF", "Waccamaw", "Delaware", "Alligator"), 
                     values = c(24, 21, 23, 22)) +
  theme(legend.position = "none")
#dev.off()



#### 6. Biogeochem. ####
# It would probably be helpful for some readers to see how sites compare in the raw BGC data, not just the correlations with CH4 or salinity.
# But, what to show? Most useful would be to show variables that were measured in at least 2 sites, that way at least some comparison can be done.

metaComb_long <- melt(frol$map_loaded,
                      id.vars = c("Estuary2", "Salt", "Depth", "sampleID"),
                      measure.vars = names(frol$map_loaded)[10:68]) %>%
  filter(variable == "CH4_ug_m2_h" |
           variable == "CO2_ug_m2_h" |
           variable == "Salinity_ppt_all" |
           variable == "pH" |
           variable == "sed_per_C" |
           variable == "sed_per_N" |
           variable == "sed_CN" |
           variable == "Cl_mgL" |
           variable == "SO4_mgL" |
           variable == "TN_mgL" |
           variable == "Fe_mgL" |
           variable == "N2O_ug_m2_h" |
           variable == "NH4_mgL" |
           variable == "PO4_mgL")

# Sort by sample size
metaComb_long_n <- metaComb_long %>%
  group_by(variable) %>%
  dplyr::summarize(n = sum(!is.na(value))) %>%
  arrange(desc(n))

metaComb_long$variable <- factor(metaComb_long$variable,
                                 levels = metaComb_long_n$variable)
metaComb_long <- subset(metaComb_long, variable != "CH4_ug_m2_h")
metaComb_long <- subset(metaComb_long, sampleID != "WestPond_CattailA_D1")
metaComb_long <- subset(metaComb_long, sampleID != "WestPond_CattailA_D2")

# SF Cl is in wrong units - divide by 35.453
metaComb_long <- metaComb_long %>%
  mutate(value2 = value/35.453)
SF_Cl <- subset(metaComb_long, Estuary2 == "SF" & variable == "Cl_mgL")
SF_Cl$value <- SF_Cl$value2
metaComb_long <- rbind(metaComb_long, SF_Cl)

facet_names <- c("CO2_ug_m2_h" = "CO2 Flux (ug/m2/h)", 
                 "Salinity_ppt_all" = "Salinity (ppt)",
                 "pH" = "pH",
                 "sed_per_C" = "% C",
                 "sed_per_N" = "% N",
                 "sed_CN" = "C:N",
                 "Cl_mgL" = "Cl (mg/L)",
                 "SO4_mgL" = "SO4 (mg/L)",
                 "TN_mgL" = "Total N (mg/L)",
                 "Fe_mgL" = "Fe (mg/L)",
                 "N2O_ug_m2_h" = "N2O Flux (ug/m2/h)",
                 "NH4_mgL" = "NH4 (mg/L)",
                 "PO4_mgL" = "PO4 (mg/L)")
strip <- strip_themed(background_x = elem_list_rect(fill = c("#FFFF99", "#B15928",
                                                             "#B15928", "#B15928",
                                                             "#A6CEE3", "#A6CEE3", 
                                                             "#A6CEE3", "#A6CEE3",
                                                             "#A6CEE3", "#A6CEE3", 
                                                             "#FFFF99", "#A6CEE3",
                                                             "#A6CEE3")))

meanY <- metaComb_long %>%
  group_by(variable) %>%
  summarise(varmax = max(value, na.rm = T)) %>%
  mutate(y = varmax + varmax/1.9)

meanlab <- metaComb_long %>%
  group_by(Estuary2, Salt, variable) %>%
  summarise(mean = round(mean(value, na.rm = T), 1),
            max = max(value, na.rm = T),
            min = min(value, na.rm = T)) %>%
  left_join(., meanY, by = "variable")

png("FinalFigs/FigureS7.png", width = 8, height = 6, units = "in", res = 300)
ggplot(metaComb_long, aes(Estuary2, value, color = Estuary2)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Salt)) +
  geom_point(size = 2, alpha = 0.75, position = position_jitterdodge(),
              aes(shape = Estuary2, colour = Salt)) +
  geom_text_repel(data = meanlab,
            aes(Estuary2, y, group = Salt, label = mean),
            position = position_dodge(width = 0.9),
            size = 2, color = "black", inherit.aes = F,
            min.segment.length = 2,
            direction = "y",
            box.padding = 0.1) +
  labs(x = "Site",
       colour = "Salinity",
       fill = "Depth (cm)",
       shape = "Site") +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(breaks = c( "SF", "Waccamaw", "Delaware", "Alligator"), 
                     values = c(24, 21, 23, 22)) +
  guides(shape = guide_legend(order = 1,),
         colour = guide_legend(order = 2),
         fill = guide_legend(override.aes = list(shape = c(16, 1)), order = 3)) +
  facet_wrap2(~ variable, scales = "free_y", labeller = as_labeller(facet_names), 
              strip = strip, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank()) +
  theme(strip.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        legend.position = "right")
dev.off()



#### 7. NCBI ####
# For NCBI submission
# Filter repset to OTU's actually analyzed
library(microseq)
f <- microseq::readFasta("data/repset.comb.fasta") %>%
  separate(Header, into = c("Header", "Junk"), sep = " ") %>%
  dplyr::select(-Junk) %>%
  filter(Header %in% rownames(frol$data_loaded))
microseq::writeFasta(f, "data/repset_used.fasta")

# Map sequences to samples - Info is in data_loaded
# This takes a while!
info <- frol$data_loaded
for (i in 1:ncol(info)) {
  for (j in 1:nrow(info)) {
    ifelse(info[j, i] > 0, info[j, i] <- names(info)[i], info[j, i] <- "")
  }
}

# Merge columns
info_cat <- info
info_cat <- info_cat %>%
  mutate_all(na_if, "") %>%
  mutate(unite(., "sample_name", c(names(info)), sep = ", ")) %>%
  mutate(sample_name = gsub("NA, ", "", sample_name)) %>%
  mutate(sample_name = gsub(", NA", "", sample_name)) %>%
  rownames_to_column(var = "Sequence_ID") %>%
  dplyr::select(Sequence_ID, sample_name)

# They don't accept multiple samples, so make another one with just the first sample
info_first <- info_cat %>%
  separate(sample_name, into = c("sample_name", "Junk"), sep = ", ") %>%
  dplyr::select(Sequence_ID, sample_name)

# Save
write_tsv(info_cat, file = "data/biosample_assignment.tsv")
write_tsv(info_first, file = "data/biosample_assignment_first.tsv")

# Filter out NCBI flagged sequences
flagged <- read.csv("data/NCBI_filter.csv")
f <- readFasta("data/repset_used.fasta") %>%
  filter(Header %notin% flagged$ASV_ID)
microseq::writeFasta(f, "data/repset_used_filtered.fasta")
sum(f$Header %in% flagged$ASV_ID)

info_cat <- read_tsv("data/biosample_assignment.tsv") %>%
  filter(Sequence_ID %notin% flagged$ASV_ID)
info_first <- read_tsv("data/biosample_assignment_first.tsv") %>%
  filter(Sequence_ID %notin% flagged$ASV_ID)
write_tsv(info_cat, file = "data/biosample_assignment_filt.tsv")
write_tsv(info_first, file = "data/biosample_assignment_first_filt.tsv")
