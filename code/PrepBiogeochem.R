# Prepare Biogeochemistry data for SF/East Coast Synthesis Project
# by Cliff Bueno de Mesquita, JGI, October 2022
# Input data files from the 4 projects, homogenize (get all same columns), and combine
# If not present, columns need to be added to each to facilitate rbind (row bind)
# Calculate approximate salinity (ppt) using 0.0018066 * Cl-
# Don't do main analysis here - merge with the input$map_loaded file in the main microbial analysis
# Variables present in all are Salinity, CH4, Cl, SO4, NH4, PO4

#### Setup ####
# Libraries
library(plyr)
library(tidyverse)
library(readxl)
library(reshape2)
library(cowplot)

# Working directory (repo)
setwd("~/Documents/GitHub/EastCoast/")

# Functions
plot_panel <- function(data, x, y) {
  ggplot(data, aes_string(x, y)) +
    geom_boxplot() +
    scale_x_discrete(labels = c("SC", "NC", "DE", "SF")) +
    labs(x = NULL,
         y = vars[i]) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8),
          legend.position = "none")
}


#### DE Biogeochem ####
# Dealt with in Excel, see Biogeochem_DE.xlsx
# Assigned sampleIDs, matched to sequence data
# Note: all porewater data
metaDE <- read_excel("~/Documents/GitHub/EastCoast/Biogeochem_DE.xlsx",
                     sheet = 2) %>%
  filter(sampleID != "NA") %>%
  mutate_at(14:33, as.numeric) %>%
  mutate(Salinity_ppt_all = Cl_mgL/1000,
         N2O = `N2O Flux umol m-2 hr-1`,
         N2O_ug_m2_h = (44.013 * N2O),
         TotalVFA_uM = `Total VFA (uM)`,
         SR_umol_cm3_d = `SR (umol cm-3 d-1)`,
         AMG_umol_cm3_d = `AMG (umol cm-3 d-1)`,
         CO2_ug_m2_h = NA,
         TOC_mgL = NA, 
         TN_mgL = NA, 
         Br_mgL = NA,
         NO3_mgL = NA, 
         DIN_mgL = NA, 
         DON_mgL = NA,
         pH = NA,
         DOC_mgL = NA,
         Na_mgL = NA,
         K_mgL = NA,
         Ca_mgL = NA,
         Mn_mgL = NA,
         Mg_mgL = NA,
         Cu_mgL = NA,
         Zn_mgL = NA,
         Estuary = "Delaware",
         Salinity = NA,
         sed_per_C = NA, 
         sed_per_N = NA, 
         sed_CN = NA, 
         sed_per_org = NA, 
         sed_per_inorg = NA,
         Conductivity_uS_cm = NA, 
         CH4_pw_air_ppmv = NA,
         N2_umol_m2_h = NA,	
         SOD_umol_m2_h = NA,	
         NO3_umol_m2_h = NA,	
         NH4_umol_m2_h = NA,	
         SRP_umol_m2_h = NA,	
         DON_umol_m2_h = NA,
         NEE_mgC_m2_m = NA,
         GEP_mgC_m2_m = NA,
         PAR_uE_m2_s = NA,
         CH4_pot_umol_gdw_h = NA, 
         CO2_pot_umol_gdw_h = NA,
         sed_pH = NA, 
         sed_NH4_mgL = NA, 
         sed_NO3_mgL = NA, 
         sed_PO4_mgL = NA, 
         sed_Cl_mgL = NA, 
         sed_SO4_mgL = NA, 
         sed_Bulk_dens = NA, 
         sed_Fe_mgL = NA, 
         sed_Mn_mgL = NA, 
         sed_Cu_mgL = NA, 
         sed_Zn_mgL = NA,
         DIC_mgL = NA) %>%
  dplyr::select(Estuary, sampleID, Salinity_calcd_ppt, Salinity, Salinity_ppt_all,
         Conductivity_uS_cm, CH4_pw_air_ppmv,
         CH4_ug_m2_h, N2O_ug_m2_h, CO2_ug_m2_h,
         NEE_mgC_m2_m, GEP_mgC_m2_m, PAR_uE_m2_s,
         CH4_pot_umol_gdw_h, CO2_pot_umol_gdw_h,
         Cl_mgL, SO4_mgL, NH4_mgL, NO3_mgL, pH, PO4_mgL, 
         Fe_mgL, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d,
         TOC_mgL, TN_mgL, Br_mgL, DIN_mgL, DON_mgL, DOC_mgL, DIC_mgL, Na_mgL,
         K_mgL, Ca_mgL, Mn_mgL, Mg_mgL, Cu_mgL, Zn_mgL,
         sed_pH, sed_NH4_mgL, sed_NO3_mgL, sed_PO4_mgL, sed_Cl_mgL, sed_SO4_mgL,
         sed_per_C, sed_per_N, sed_CN, sed_per_org, sed_per_inorg,
         sed_Bulk_dens, sed_Fe_mgL, sed_Mn_mgL, sed_Cu_mgL, sed_Zn_mgL,
         N2_umol_m2_h,	SOD_umol_m2_h,	NO3_umol_m2_h,	NH4_umol_m2_h,	
         SRP_umol_m2_h,	DON_umol_m2_h, Porosity)



#### NC Biogeochem ####
# Marcelo/Emily sent three files, water chem, GHG, and pH
# Note: all porewater
# Data just for incubation, not field samples
# Take just flooded samples,  final timepoint, 5 and 15 cm depth
# Note: _AF1, _AF3, _AF4, BF3, BF4, BF5 sequence IDs not here and prob should be removed
# Note: SO4 trt has 2, 3, 4, 5a, 5b instead of 1,2,3,4,5. can only match 2,3,4
NC_chem <- read_excel("~/Documents/GitHub/EastCoast/Biogeochem_NC.xlsx",
                      sheet = 1)
NC_GHG <- read_excel("~/Documents/GitHub/EastCoast/Biogeochem_NC.xlsx",
                     sheet = 2)
NC_pH <- read_excel("~/Documents/GitHub/EastCoast/Soil pH-Jessie.xls",
                    sheet = 2) %>%
  group_by(sampleID) %>%
  dplyr::summarize(pH = mean(pH)) %>%
  filter(sampleID != "Not sequenced")
metaNC <- left_join(NC_chem, NC_GHG, by = "bgcID") %>%
  filter(`Hydrology 2` == "Flooded",
         Week.x == 14,
         Depth == "5cm" | Depth == "15cm") %>%
  mutate(Treatment = recode_factor(`Water Chemistry.x`,
                                   "DI" = "DI_ctrl",
                                   "ASW-SO4" = "ASW_noS")) %>%
  mutate(Prefix = "TL_inc",
         Depth2 = ifelse(Depth == "5cm", "d1", "d2")) %>%
  mutate(sampleID = paste(Prefix, Depth2, Treatment, Rep.x, sep = "_")) %>%
  left_join(., NC_pH, by = "sampleID") %>%
  mutate(TOC_mgL = `NPOCmg/L`,
         TN_mgL = TNmgL,
         NH4_mgL = NH4mgL,
         PO4_mgL = PO4mgL,
         Cl_mgL = ClppmmgL,
         SO4_mgL = SO4mgL,
         Br_mgL = BrmgL,
         NO3_mgL = NO3NmgL,
         DIN_mgL = `DIN (mg/L)`,
         DON_mgL = `DON (mg/L)`,
         CH4_ug_m2_h = `CH4 (mg CH4/m2/hr)` * 1000,
         CO2_ug_m2_h = `CO2 flux (mg/m2/hr)` * 1000,
         N2O_ug_m2_h = `N2O flux (ug/m2/hr)`,
         Salinity_calcd_ppt = Cl_mgL * 0.0018066,
         Salinity_ppt_all = Cl_mgL/1000) %>%
  mutate(Fe_mgL = NA, 
         Acetate_mgL = NA, 
         TotalVFA_uM = NA,
         SR_umol_cm3_d = NA,
         AMG_umol_cm3_d = NA,
         DOC_mgL = NA,
         Na_mgL = NA,
         K_mgL = NA,
         Ca_mgL = NA,
         Mn_mgL = NA,
         Mg_mgL = NA,
         Cu_mgL = NA,
         Zn_mgL = NA,
         Estuary = "Alligator",
         Salinity = NA,
         sed_per_C = NA, 
         sed_per_N = NA, 
         sed_CN = NA, 
         sed_per_org = NA, 
         sed_per_inorg = NA,
         Conductivity_uS_cm = NA, 
         CH4_pw_air_ppmv = NA,
         N2_umol_m2_h = NA,	
         SOD_umol_m2_h = NA,	
         NO3_umol_m2_h = NA,	
         NH4_umol_m2_h = NA,	
         SRP_umol_m2_h = NA,	
         DON_umol_m2_h = NA,
         NEE_mgC_m2_m = NA,
         GEP_mgC_m2_m = NA,
         PAR_uE_m2_s = NA,
         CH4_pot_umol_gdw_h = NA, 
         CO2_pot_umol_gdw_h = NA,
         sed_pH = NA, 
         sed_NH4_mgL = NA, 
         sed_NO3_mgL = NA, 
         sed_PO4_mgL = NA, 
         sed_Cl_mgL = NA, 
         sed_SO4_mgL = NA, 
         sed_Bulk_dens = NA, 
         sed_Fe_mgL = NA, 
         sed_Mn_mgL = NA, 
         sed_Cu_mgL = NA, 
         sed_Zn_mgL = NA,
         DIC_mgL = NA,
         Porosity = NA) %>%
  dplyr::select(Estuary, sampleID, Salinity_calcd_ppt, Salinity, Salinity_ppt_all,
                Conductivity_uS_cm, CH4_pw_air_ppmv,
                CH4_ug_m2_h, N2O_ug_m2_h, CO2_ug_m2_h,
                NEE_mgC_m2_m, GEP_mgC_m2_m, PAR_uE_m2_s,
                CH4_pot_umol_gdw_h, CO2_pot_umol_gdw_h,
                Cl_mgL, SO4_mgL, NH4_mgL, NO3_mgL, pH, PO4_mgL, 
                Fe_mgL, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d,
                TOC_mgL, TN_mgL, Br_mgL, DIN_mgL, DON_mgL, DOC_mgL, DIC_mgL, Na_mgL,
                K_mgL, Ca_mgL, Mn_mgL, Mg_mgL, Cu_mgL, Zn_mgL,
                sed_pH, sed_NH4_mgL, sed_NO3_mgL, sed_PO4_mgL, sed_Cl_mgL, sed_SO4_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_per_org, sed_per_inorg,
                sed_Bulk_dens, sed_Fe_mgL, sed_Mn_mgL, sed_Cu_mgL, sed_Zn_mgL,
                N2_umol_m2_h,	SOD_umol_m2_h,	NO3_umol_m2_h,	NH4_umol_m2_h,	
                SRP_umol_m2_h,	DON_umol_m2_h, Porosity)



#### SC Biogeochem ####
# Note: flux, porewater, and soil data, so need to specify
# Also 2 types of flux: field measured emissions and potential production
# Other fluxes: N2_umol_m2_h,	SOD_umol_m2_h,	NO3_umol_m2_h,	NH4_umol_m2_h,	SRP_umol_m2_h,	DON_umol_m2_h
# Porewater pH, conductivity, salinity, CH4, DIC
# Soil %C, %N, CN, % organic, % inorganic
metaSC <- read_excel("~/Documents/GitHub/EastCoast/Nov 2011 data from Brookgreen.xlsx",
                     sheet = 5) %>%
  mutate(Conductivity_uS_cm = as.numeric(Conductivity_uS_cm),
         Salinity_calcd_ppt = NA,
         Salinity_ppt_all = as.numeric(Salinity),
         CO2_ug_m2_h = CO2_mg_m2_m * 1000 * 60,
         CH4_ug_m2_h = CH4_mg_m2_m * 1000 * 60,
         CH4_pw_air_ppmv = as.numeric(CH4_pw_air_ppmv),
         DIC_mmolL = as.numeric(DIC_mmolL),
         CH4_pot_umol_gdw_h = CH4_pot_nmol_gdw_h/1000) %>%
  mutate(Cl_mgL = NA,
         SO4_mgL = NA,
         NH4_mgL = NA,
         PO4_mgL = NA,
         TN_mgL = NA,
         Fe_mgL = NA, 
         Acetate_mgL = NA, 
         TotalVFA_uM = NA,
         SR_umol_cm3_d = NA,
         AMG_umol_cm3_d = NA,
         NO3_mgL = NA,
         DOC_mgL = NA,
         DIC_mgL = DIC_mmolL * 12.01, 
         Na_mgL = NA,
         K_mgL = NA,
         Ca_mgL = NA,
         Mn_mgL = NA,
         Mg_mgL = NA,
         Cu_mgL = NA,
         Zn_mgL = NA,
         N2O_ug_m2_h = NA,
         Acetate_mgL = NA, 
         TotalVFA_uM = NA,
         SR_umol_cm3_d = NA,
         AMG_umol_cm3_d = NA,
         TOC_mgL = NA, 
         Br_mgL = NA,
         DIN_mgL = NA, 
         DON_mgL = NA,
         Estuary = "Waccamaw",
         sed_pH = NA, 
         sed_NH4_mgL = NA, 
         sed_NO3_mgL = NA, 
         sed_PO4_mgL = NA, 
         sed_Cl_mgL = NA, 
         sed_SO4_mgL = NA, 
         sed_Bulk_dens = NA, 
         sed_Fe_mgL = NA, 
         sed_Mn_mgL = NA, 
         sed_Cu_mgL = NA, 
         sed_Zn_mgL = NA,
         Porosity = NA) %>%
  dplyr::select(Estuary, sampleID, Salinity_calcd_ppt, Salinity, Salinity_ppt_all,
                Conductivity_uS_cm, CH4_pw_air_ppmv,
                CH4_ug_m2_h, N2O_ug_m2_h, CO2_ug_m2_h,
                NEE_mgC_m2_m, GEP_mgC_m2_m, PAR_uE_m2_s,
                CH4_pot_umol_gdw_h, CO2_pot_umol_gdw_h,
                Cl_mgL, SO4_mgL, NH4_mgL, NO3_mgL, pH, PO4_mgL, 
                Fe_mgL, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d,
                TOC_mgL, TN_mgL, Br_mgL, DIN_mgL, DON_mgL, DOC_mgL, DIC_mgL, Na_mgL,
                K_mgL, Ca_mgL, Mn_mgL, Mg_mgL, Cu_mgL, Zn_mgL,
                sed_pH, sed_NH4_mgL, sed_NO3_mgL, sed_PO4_mgL, sed_Cl_mgL, sed_SO4_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_per_org, sed_per_inorg,
                sed_Bulk_dens, sed_Fe_mgL, sed_Mn_mgL, sed_Cu_mgL, sed_Zn_mgL,
                N2_umol_m2_h,	SOD_umol_m2_h,	NO3_umol_m2_h,	NH4_umol_m2_h,	
                SRP_umol_m2_h,	DON_umol_m2_h, Porosity)


#### SF Biogeochem ####
# Get meta file from other repo
# Get porewater data (but NO3 and NH4 are from sediment extraction)
# Get correct units, calc salinity
# Note Cl was in meq/L so needed to convert
metaSF <- read.delim("~/Documents/GitHub/SF_microbe_methane/data/meta/SF_sal_meta_FIX3.5.txt") %>%
  mutate(Cl_mgL = Cl_pw*35.45) %>%
  mutate(Salinity_calcd_ppt = Salinity.x,
         Salinity_ppt_all = Salinity.x,
         Salinity = Salinity.x,
         sampleID = Sample,
         CO2_ug_m2_h = CO2_mg_m2_h * 1000,
         SO4_mgL = SO4_pw,
         sed_pH = pH,
         sed_NH4_mgL = NH4_N,
         sed_NO3_mgL = NO3_N,
         sed_PO4_mgL = Olsen_P,
         sed_Cl_mgL = Cl*35.45,
         sed_SO4_mgL = SO4,
         Fe_mgL = Fe_pw,
         TN_mgL = N,
         DOC_mgL = DOC_mg_L,
         Na_mgL = Na_pw,
         K_mgL = K_pw,
         Ca_mgL = Ca_pw,
         Mn_mgL = Mn_pw,
         Mg_mgL = Mg_pw,
         Cu_mgL = Cu_pw,
         Zn_mgL = Mg_pw) %>%
  mutate(N2O_ug_m2_h = NA,
         Acetate_mgL = NA, 
         TotalVFA_uM = NA,
         SR_umol_cm3_d = NA,
         AMG_umol_cm3_d = NA,
         TOC_mgL = NA, 
         Br_mgL = NA,
         DIN_mgL = NA, 
         DON_mgL = NA,
         NH4_mgL = NA,
         NO3_mgL = NA,
         PO4_mgL = NA,
         Estuary = "SF",
         sed_per_C = C, 
         sed_per_N = N, 
         sed_CN = CN, 
         sed_P = P,
         sed_CP = CP,
         sed_NP = NP,
         sed_per_org = NA, 
         sed_per_inorg = NA,
         sed_Bulk_dens = Bulk_dens,
         sed_Fe_mgL = Fe,
         sed_Mn_mgL = Mn,
         sed_Cu_mgL = Cu,
         sed_Zn_mgL = Zn,
         Conductivity_uS_cm = NA, 
         CH4_pw_air_ppmv = NA,
         N2_umol_m2_h = NA,	
         SOD_umol_m2_h = NA,	
         NO3_umol_m2_h = NA,	
         NH4_umol_m2_h = NA,	
         SRP_umol_m2_h = NA,	
         DON_umol_m2_h = NA,
         NEE_mgC_m2_m = NA,
         GEP_mgC_m2_m = NA,
         PAR_uE_m2_s = NA,
         CH4_pot_umol_gdw_h = NA, 
         CO2_pot_umol_gdw_h = NA,
         DIC_mgL = NA,
         Porosity = NA) %>%
  dplyr::select(Estuary, sampleID, Salinity_calcd_ppt, Salinity, Salinity_ppt_all,
                Conductivity_uS_cm, CH4_pw_air_ppmv,
                CH4_ug_m2_h, N2O_ug_m2_h, CO2_ug_m2_h,
                NEE_mgC_m2_m, GEP_mgC_m2_m, PAR_uE_m2_s,
                CH4_pot_umol_gdw_h, CO2_pot_umol_gdw_h,
                Cl_mgL, SO4_mgL, NH4_mgL, NO3_mgL, pH, PO4_mgL, 
                Fe_mgL, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d,
                TOC_mgL, TN_mgL, Br_mgL, DIN_mgL, DON_mgL, DOC_mgL, DIC_mgL, Na_mgL,
                K_mgL, Ca_mgL, Mn_mgL, Mg_mgL, Cu_mgL, Zn_mgL,
                sed_pH, sed_NH4_mgL, sed_NO3_mgL, sed_PO4_mgL, sed_Cl_mgL, sed_SO4_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_per_org, sed_per_inorg,
                sed_Bulk_dens, sed_Fe_mgL, sed_Mn_mgL, sed_Cu_mgL, sed_Zn_mgL,
                N2_umol_m2_h,	SOD_umol_m2_h,	NO3_umol_m2_h,	NH4_umol_m2_h,	
                SRP_umol_m2_h,	DON_umol_m2_h, Porosity)

# Something is of about salinity and Cl
m <- lm(Cl_mgL/1000 ~ Salinity, data = metaSF)
summary(m)
ggplot(metaSF, aes(Salinity, Cl_mgL/1000)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Water salinity (YSI ppt)",
       y = "Porewater salinity (Cl- ppt)") +
  theme_classic()

m <- lm(sed_Cl_mgL/1000 ~ Salinity, data = metaSF)
summary(m)
ggplot(metaSF, aes(Salinity, sed_Cl_mgL/1000)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Water salinity (YSI ppt)",
       y = "Sediment salinity (Cl- ppt)") +
  theme_classic()

#### Merge and Save ####
sum(names(metaDE) != names(metaNC))
sum(names(metaDE) != names(metaSC))
sum(names(metaDE) != names(metaSF))
metaComb <- rbind(metaDE, metaNC, metaSC, metaSF)

write.csv(metaComb, "biogeochem_all_clean.csv")

metaComb <- read.csv("biogeochem_all_clean.csv") %>%
  mutate(Estuary = factor(Estuary,
                          levels = c("Waccamaw", "Alligator", "Delaware", "SF"))) %>%
  dplyr::select(-Salinity_calcd_ppt, -Salinity)



#### Quick Look ####
# Start to look at some biogeochem
# Plots:
# All, key vars, methane, methane/salinity

#### _All ####
metaComb_long <- melt(metaComb,
                      id.vars = c("Estuary", "sampleID"),
                      measure.vars = names(metaComb)[4:ncol(metaComb)])

# Sort by sample size
metaComb_long_n <- metaComb_long %>%
  group_by(variable) %>%
  dplyr::summarize(n = sum(!is.na(value))) %>%
  arrange(desc(n))

metaComb_long$variable <- factor(metaComb_long$variable,
                                 levels = metaComb_long_n$variable)

pdf("InitialFigs/Comb_All_Biogeochem.pdf", width = 10, height = 10)
ggplot(metaComb_long, aes(Estuary, value, color = Estuary)) +
  geom_boxplot() +
  scale_color_viridis_d() +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank()) +
  theme(strip.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        legend.position = "none")
dev.off()

#### _Key ####
# Some key partially shared variables (flux, CNPS, salinity, pH)
# CH4 and CO2 flux
# Salinity, pH, SO4, sed SO4
# TOC, DIC, DOC
# TN, DIN, DON
# NH4, sed_NH4, NO3, sed_NO3
# PO4, sed_PO4
# % C, % N, C:N
metaComb_long_shared <- metaComb_long %>%
  filter(variable == "CH4_ug_m2_h" |
           variable == "CO2_ug_m2_h" |
           variable == "Salinity_ppt_all" |
           variable == "pH" |
           variable == "SO4_mgL" |
           variable == "sed_SO4_mgL" |
           variable == "TOC_mgL" |
           variable == "DIC_mgL" |
           variable == "DOC_mgL" |
           variable == "TN_mgL" |
           variable == "DIN_mgL" |
           variable == "DON_mgL" |
           variable == "NH4_mgL" |
           variable == "sed_NH4_mgL" |
           variable == "NO3_mgL" |
           variable == "sed_NO3_mgL" |
           variable == "PO4_mgL" |
           variable == "sed_PO4_mgL" |
           variable == "sed_per_C" |
           variable == "sed_per_N" |
           variable == "sed_CN") %>%
  mutate(variable = factor(variable,
                           levels = c("CH4_ug_m2_h", "CO2_ug_m2_h", 
                                      "Salinity_ppt_all", "pH", "SO4_mgL", 
                                      "sed_SO4_mgL", "TOC_mgL", "DIC_mgL", 
                                      "DOC_mgL", "TN_mgL", "DIN_mgL", "DON_mgL", 
                                      "NH4_mgL", "sed_NH4_mgL", "NO3_mgL", 
                                      "sed_NO3_mgL", "PO4_mgL", "sed_PO4_mgL",
                                      "sed_per_C", "sed_per_N", "sed_CN")))

fw <- ggplot(metaComb_long_shared, aes(Estuary, value, color = Estuary)) +
  geom_boxplot() +
  scale_color_viridis_d() +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.position = "right")
leg <- get_legend(fw)

vars <- levels(metaComb_long_shared$variable)
p <- list()
for (i in 1:length(vars)) {
  p[[i]] <- plot_panel(data = metaComb, 
                       x = "Estuary", 
                       y = vars[i])
}

pdf("InitialFigs/Comb_All_Biogeochem_CNPS.pdf", width = 8, height = 8)
plot_grid(p[[1]], p[[2]], NULL, NULL,
          p[[3]], p[[4]], p[[5]], p[[6]],
          p[[7]], p[[8]], p[[9]], NULL,
          p[[10]], p[[11]], p[[12]], NULL,
          p[[13]], p[[14]], p[[15]], p[[16]],
          p[[17]], p[[18]], NULL, NULL,
          p[[19]], p[[20]], p[[21]], NULL,
          align = "hv",
          ncol = 4)
dev.off()



#### _CH4 ####
ggplot(metaComb, aes(Salinity_ppt_all, CH4_ug_m2_h)) +
  geom_point(aes(color = Estuary)) +
  geom_smooth() +
  scale_color_viridis_d() +
  labs(x = "Salinity (ppt)",
       y = "CH4 (ug/m2/h)") +
  theme_bw()

ggplot(metaComb, aes(Salinity_ppt_all, CH4_ug_m2_h)) +
  geom_point(aes(color = Estuary)) +
  geom_smooth() +
  scale_y_log10() +
  scale_color_viridis_d() +
  labs(x = "Salinity (ppt)",
       y = "CH4 (ug/m2/h)") +
  theme_bw()

ggplot(metaComb, aes(Estuary, CH4_ug_m2_h)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank())

pdf("InitialFigs/Comb_All_CH4.pdf", width = 7, height = 5)
ggplot(metaComb, aes(Estuary, CH4_ug_m2_h, colour = Estuary)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_y_log10() +
  scale_color_viridis_d() +
  labs(x = "Estuary",
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

#### _CH4/Salinity ####
# Show trends by estuary
pdf("InitialFigs/Comb_All_CH4_Salinity.pdf", width = 8, height = 3)
ggplot(metaComb, aes(Salinity_ppt_all, CH4_ug_m2_h)) +
  geom_point(aes(color = Estuary)) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_viridis_d() +
  facet_wrap(~ Estuary, ncol = 4) +
  labs(x = "Salinity (ppt)",
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# Subset SF data to the lower salinity sites
# Try mesohaline < 18 ppt and Oligohaline < 5 ppt, and max NC or max DE)
max(metaDE$Salinity_ppt_all, na.rm = TRUE)
max(metaNC$Salinity_ppt_all, na.rm = TRUE)
max(metaSC$Salinity_ppt_all, na.rm = TRUE)
max(metaSF$Salinity_ppt_all, na.rm = TRUE)

metaComb_18ppt <- subset(metaComb, Salinity_ppt_all <= 18)
pdf("InitialFigs/Comb_All_CH4_Salinity_18ppt.pdf", width = 8, height = 3)
ggplot(metaComb_18ppt, aes(Salinity_ppt_all, CH4_ug_m2_h)) +
  geom_point(aes(color = Estuary)) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_viridis_d() +
  facet_wrap(~ Estuary, ncol = 4) +
  labs(x = "Salinity (ppt)",
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

metaComb_5ppt <- subset(metaComb, Salinity_ppt_all <= 5)
pdf("InitialFigs/Comb_All_CH4_Salinity_5ppt.pdf", width = 8, height = 3)
ggplot(metaComb_5ppt, aes(Salinity_ppt_all, CH4_ug_m2_h)) +
  geom_point(aes(color = Estuary)) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_viridis_d() +
  facet_wrap(~ Estuary, ncol = 4) +
  labs(x = "Salinity (ppt)",
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

metaComb_12.5ppt <- subset(metaComb, Salinity_ppt_all <= max(metaNC$Salinity_ppt_all, na.rm = TRUE))
pdf("InitialFigs/Comb_All_CH4_Salinity_12.5ppt.pdf", width = 8, height = 3)
ggplot(metaComb_12.5ppt, aes(Salinity_ppt_all, CH4_ug_m2_h)) +
  geom_point(aes(color = Estuary)) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_viridis_d() +
  facet_wrap(~ Estuary, ncol = 4) +
  labs(x = "Salinity (ppt)",
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

metaComb_7.5ppt <- subset(metaComb, Salinity_ppt_all <= max(metaDE$Salinity_ppt_all, na.rm = TRUE))
pdf("InitialFigs/Comb_All_CH4_Salinity_7.5ppt.pdf", width = 8, height = 3)
ggplot(metaComb_7.5ppt, aes(Salinity_ppt_all, CH4_ug_m2_h)) +
  geom_point(aes(color = Estuary)) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_viridis_d() +
  facet_wrap(~ Estuary, ncol = 4) +
  labs(x = "Salinity (ppt)",
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

#### _CH4/CO2 ####
pdf("InitialFigs/Comb_All_CH4_CO2.pdf", width = 8, height = 3)
ggplot(metaComb, aes(CO2_ug_m2_h, CH4_ug_m2_h)) +
  geom_point(aes(color = Estuary)) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_viridis_d() +
  facet_wrap(~ Estuary, ncol = 4) +
  labs(x = expression(""*CO[2]*" flux (µg/"*m^2*"/h)"),
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"))
dev.off()

