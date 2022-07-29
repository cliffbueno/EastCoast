# East Coast 16S data analysis
# Analyze microbial communities from samples from Weston, Neubauer, Bernhardt labs
# Samples are from Delaware River, Alligator River NC, and Waccamaw River SC
# Compare communities to Tringe Lab SF Bay/Delta samples
# Two plates were sequenced and processed with iTagger
# Cliff reassigned taxonomy with SILVA v 138.1
# Analysis by Cliff Bueno de Mesquita Summer 2022



# Key papers/experimental information:
# Neubauer 2013 Estuaries and Coasts
## - South Carolina, experiment with Control, +salt, +fresh plots

# Ardón et al. 2013 Global Change Biology
## - North Carolina, experiment with hydrology, salt, sulfate microcosms
## 50 samples: 10 initial, 40 experimental (2 x 2 x 2 factorial)
## Control = DI water
## ASW = artificial salt-water with 5 ppt salinity
## ASW-SO4 = artificial salt-water with 5 ppt salinity without sulfate
## DI+SO4 = DI water with same sulfate as ASW
## Flooded samples (water level maintained at surface)
## There was also a drought treatment but it wasn't sequenced
## 30˚C, 12 weeks

# Weston et al. 2014 Biogeochemistry
## - Delaware River, field sampling
## Tidal freshwater, oligohaline, mesohaline

## - Delaware River, experiment
## Freshwater or ASW, 0, 4 wk, 7 wk, 12 wk sampling

## - Delaware River, transplants
## Transplanted TFM to another TFM, Oligo, and Meso at surface and 40 cm below surface



#### Setup ####
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

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
`%notin%` <- Negate(`%in%`)
paste_ranks = function(sm_taxa){
  k = data.frame(k ="k__", sm_taxa['taxonomy1'])
  k2 <- do.call(paste, c(k, sep = ""))
  
  p = data.frame(k ="p__", sm_taxa['taxonomy2'])
  p2 <- do.call(paste, c(p, sep = ""))
  
  c = data.frame(k ="c__", sm_taxa['taxonomy3'])
  c2 <- do.call(paste, c(c, sep = ""))
  
  o = data.frame(k ="o__", sm_taxa['taxonomy4'])
  o2 <- do.call(paste, c(o, sep = ""))
  
  f = data.frame(k ="f__", sm_taxa['taxonomy5'])
  f2 <- do.call(paste, c(f, sep = ""))
  
  g = data.frame(k ="g__", sm_taxa['taxonomy6'])
  g2 <- do.call(paste, c(g, sep = ""))
  
  # NOT USING SPECIES HERE! OTU preprocessing doesn't!
  s = data.frame(k ="s__", sm_taxa['taxonomy7'])
  s2 <- do.call(paste, c(s, sep = ""))
  
  # combine all
  lineage_df = data.frame(k2, p2, c2, o2, f2, g2)
  lineage = do.call(paste, c(lineage_df, sep = ';'))
  return(lineage)
}

# Guild subsetting module from other repository
source("~/Documents/GitHub/SF_microbe_methane/modules/3_OTU_subsetting_modules_v.0.4_strip.r")

# Repository path
setwd("~/Documents/GitHub/EastCoast/")

# Wyatt Hartman's guild color palette
# Note that MeOB don't exist in this dataset, so removed
# Extra methanogen guilds added so colors added too
Guild_cols <- read.table("~/Documents/GitHub/SF_microbe_methane/data/colors/Guild_color_palette.txt",
                         sep='\t') %>%
  dplyr::select(Guild, G_index, color) %>%
  set_names(c("Guild", "Index", "color")) %>%
  mutate(Index = rev(Index)) %>%
  filter(Guild != "MeOB")
Guild_cols[15,] = c("CH4_me", 16, "#FDC086")
Guild_cols[16,] = c("CH4_mix", 17, "#FFFF99")
Guild_cols <- Guild_cols %>%
  mutate(Index = as.integer(Index)) %>%
  arrange(Index)



#### _Delaware ####
# Make mapping file using the 2 Excel files that Wyatt sent Tijana for iTag sequencing
# Note, already fixed "Freswater" to "Freshwater" typo in Excel
p1 <- read_excel("SPITS Wyatt 1520 itags4.xlsx")
p2 <- read_excel("SPITS Wyatt 1520 itags pl2v2.xlsx")
metadata <- rbind(p1, p2) %>%
  dplyr::select(`Sample Name*`, `Collection Year*`, `Collection Month*`, `Collection Day*`,
         `Sample Isolated From*`, `Collection Site or Growth Conditions`, 
         `Latitude*`, `Longitude*`, `Altitude or Depth*`) %>%
  set_names(c("sampleID", "Year", "Month", "Day", "Experiment", "Treatment", 
              "Latitude", "Longitude", "Depth")) %>%
  mutate(Estuary = "NA")
for (i in 1:nrow(metadata)) {
  if (metadata$Latitude[i] == 33.5250) {
    metadata$Estuary[i] <- "Waccamaw"
  }
  if (metadata$Latitude[i] == 35.9061) {
    metadata$Estuary[i] <- "Alligator"
  }
  if (metadata$Latitude[i] > 36) {
    metadata$Estuary[i] <- "Delaware"
  }
}
write.table(metadata, "delaware_metadata.txt", sep = "\t", row.names = F)

# Metadata has 184 samples
# ASV table has 177 samples
# 7 samples lost/not sequenced - check which ones
otu_table <- read.table("seqtab_wTax_mctoolsr.txt", header = 2)
missing <- metadata %>%
  filter(sampleID %notin% names(otu_table))
write.table(missing, "not_sequenced.txt", sep = "\t", row.names = F)

# Import Data (n = 177), Filter, Rarefy, Calc richness
tax_table_fp <- "seqtab_wTax_mctoolsr.txt"
map_fp <- "delaware_metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp)

# Filter chloroplast, mitochondria, eukaryotes, unassigned at domain
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 292 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 787 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Eukarya") # none
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 34 removed

# Guilds
# Use updated Wyatt guild calling script
# Version that pulls out 4 different methanogen guilds
# Need to first reorganize OTU table
# Then place guilds as a column in input$taxonomy_loaded to analyze as any other taxonomic group
# Follow Wyatt pipeline to prep an mctoolsr object for guild analysis
taxa = input_filt$taxonomy_loaded
OTUs = input_filt$data_loaded
raw_OTU_table = data.frame(OTUs, taxa)
taxa[taxa=='NA'] <- ""
Consensus.lineage = paste_ranks(taxa)
reformed_OTU_table = data.frame(OTUs, Consensus.lineage) %>%
  mutate_if(is.integer, as.numeric) %>%
  mutate(OTU = rownames(.)) %>%
  dplyr::select(OTU, everything()) %>%
  left_join(., input_filt$taxonomy_loaded, by = c("OTU" = "taxonomy8")) %>%
  dplyr::rename(Kingdom = taxonomy1,
                Phylum= taxonomy2,
                Class = taxonomy3,
                Order = taxonomy4,
                Family = taxonomy5,
                Genus = taxonomy6) %>%
  dplyr::select(-taxonomy7) %>%
  mutate(Taxonomy = Phylum)
rownames(reformed_OTU_table) <- reformed_OTU_table$OTU
Guild_OTUs <- Get_16S_Guilds_alt(reformed_OTU_table)
levels(as.factor(Guild_OTUs$Guild))

# Now add as 9th column to input_filt$taxonomy_loaded
input_filt$taxonomy_loaded <- input_filt$taxonomy_loaded %>%
  left_join(., Guild_OTUs, by = c("taxonomy8" = "OTU")) %>%
  rename(taxonomy9 = Guild) %>%
  mutate(taxonomy9 = as.character(taxonomy9)) %>%
  mutate(taxonomy9 = replace_na(taxonomy9, "NA"))
rownames(input_filt$taxonomy_loaded) <- input_filt$taxonomy_loaded$taxonomy8
saveRDS(input_filt, "input_filt.rds")

# Rarefy at minimum (26429)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded))
se(colSums(input_filt$data_loaded))
# Depth 92024 ± 2919
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 26429)
sort(colSums(input_filt_rare$data_loaded))

# OTU Richness
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, 
                                              MARGIN = 2)

# Shannon diversity
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", 
                                                MARGIN = 2)

# Guilds
taxa = input_filt_rare$taxonomy_loaded
OTUs = input_filt_rare$data_loaded
raw_OTU_table = data.frame(OTUs, taxa)
taxa[taxa=='NA'] <- ""
Consensus.lineage = paste_ranks(taxa)
reformed_OTU_table = data.frame(OTUs, Consensus.lineage) %>%
  mutate_if(is.integer, as.numeric) %>%
  mutate(OTU = rownames(.)) %>%
  dplyr::select(OTU, everything()) %>%
  left_join(., input_filt_rare$taxonomy_loaded, by = c("OTU" = "taxonomy8")) %>%
  dplyr::rename(Kingdom = taxonomy1,
                Phylum= taxonomy2,
                Class = taxonomy3,
                Order = taxonomy4,
                Family = taxonomy5,
                Genus = taxonomy6) %>%
  dplyr::select(-taxonomy7) %>%
  mutate(Taxonomy = Phylum)
rownames(reformed_OTU_table) <- reformed_OTU_table$OTU
Guild_OTUs <- Get_16S_Guilds_alt(reformed_OTU_table)
levels(as.factor(Guild_OTUs$Guild))

# Now add as 9th column to input_filt_rare$taxonomy_loaded
input_filt_rare$taxonomy_loaded <- input_filt_rare$taxonomy_loaded %>%
  left_join(., Guild_OTUs, by = c("taxonomy8" = "OTU")) %>%
  rename(taxonomy9 = Guild) %>%
  mutate(taxonomy9 = as.character(taxonomy9)) %>%
  mutate(taxonomy9 = replace_na(taxonomy9, "NA"))
rownames(input_filt_rare$taxonomy_loaded) <- input_filt_rare$taxonomy_loaded$taxonomy8

# Save
saveRDS(input_filt_rare, "input_filt_rare.rds")



#### _Combined ####
# Make combined metadata table
de <- read.delim("delaware_metadata.txt") %>%
  dplyr::select(sampleID, Experiment, Treatment, Depth, Estuary) %>%
  set_names(c("sampleID", "Site", "Detail", "Depth", "Estuary")) %>%
  mutate(Info = Detail)
sf <- read.delim("~/Desktop/Wyatt Manuscript/metadata_mctoolsr.txt") %>%
  dplyr::select(SampleID, Location, Vegetation, Depth) %>%
  set_names(c("sampleID", "Site", "Detail", "Depth")) %>%
  mutate(Estuary = "SF") %>%
  mutate(Info = Site)
comb <- rbind(de, sf)

# Metadata has 352 samples
# ASV table has 345 samples
# 7 lost from Delaware

write.table(comb, "combined_metadata.txt", sep = "\t", row.names = F)

# Import Data (n = 345), Filter, Rarefy, Calc richness
tax_table_fp <- "seqtab_wTax_mctoolsr_comb.txt"
map_fp <- "combined_metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp)

# Filter chloroplast, mitochondria, eukaryotes, unassigned at domain
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 368 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 815 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Eukarya") # none
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 54 removed

# Rarefy at 26429
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded))
se(colSums(input_filt$data_loaded))
# Original depth 118864 ± 2323
# Drop Sandmound_TuleB_D1 (1296 reads) and Muzzi_PWB_D2 (5287 reads)
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 26429) # Now n = 343
sort(colSums(input_filt_rare$data_loaded))

# OTU Richness
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, 
                                              MARGIN = 2)

# Shannon diversity
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", 
                                                MARGIN = 2)

# Guilds
taxa = input_filt_rare$taxonomy_loaded
OTUs = input_filt_rare$data_loaded
raw_OTU_table = data.frame(OTUs, taxa)
taxa[taxa=='NA'] <- ""
Consensus.lineage = paste_ranks(taxa)
reformed_OTU_table = data.frame(OTUs, Consensus.lineage) %>%
  mutate_if(is.integer, as.numeric) %>%
  mutate(OTU = rownames(.)) %>%
  dplyr::select(OTU, everything()) %>%
  left_join(., input_filt_rare$taxonomy_loaded, by = c("OTU" = "taxonomy8")) %>%
  dplyr::rename(Kingdom = taxonomy1,
                Phylum= taxonomy2,
                Class = taxonomy3,
                Order = taxonomy4,
                Family = taxonomy5,
                Genus = taxonomy6) %>%
  dplyr::select(-taxonomy7) %>%
  mutate(Taxonomy = Phylum)
rownames(reformed_OTU_table) <- reformed_OTU_table$OTU
Guild_OTUs <- Get_16S_Guilds_alt(reformed_OTU_table)
levels(as.factor(Guild_OTUs$Guild))

# Now add as 9th column to input_filt_rare$taxonomy_loaded
input_filt_rare$taxonomy_loaded <- input_filt_rare$taxonomy_loaded %>%
  left_join(., Guild_OTUs, by = c("taxonomy8" = "OTU")) %>%
  rename(taxonomy9 = Guild) %>%
  mutate(taxonomy9 = as.character(taxonomy9)) %>%
  mutate(taxonomy9 = replace_na(taxonomy9, "NA"))
rownames(input_filt_rare$taxonomy_loaded) <- input_filt_rare$taxonomy_loaded$taxonomy8

# Save
saveRDS(input_filt_rare, "input_filt_rare_comb.rds")



#### ................................. ####
#### East Coast Overview ####
input_filt_rare <- readRDS("input_filt_rare.rds")



#### Map ####
# Summarize data, get unique coordinates
coords <- input_filt_rare$map_loaded %>%
  group_by(Latitude, Longitude) %>%
  slice_head()
min(coords$Latitude)
max(coords$Latitude)
min(coords$Longitude)
max(coords$Longitude)

del <- get_stamenmap(bbox = c(left = min(coords$Longitude) - 0.5, 
                              bottom = min(coords$Latitude) - 0.5, 
                              right = max(coords$Longitude) + 0.5, 
                              top = max(coords$Latitude) + 0.5),
                     zoom = 10, 
                     maptype = "terrain-background")
del_attributes <- attributes(del)
del_transparent <- matrix(adjustcolor(del, alpha.f = 0.4), nrow = nrow(del))
attributes(del_transparent) <- del_attributes

pdf("Figs/Map_EastCoast.pdf", width = 6, height = 6)
ggmap(del_transparent, extent = "device") + # the base map
  geom_point(data = coords,
             aes(x = Longitude, y = Latitude), size = 4) +
#  geom_text(aes(x = -121.7, y = 38.2, label = "Delta"), 
#            colour = "black", size = 6, fontface = "italic", check_overlap = T) +
#  geom_text(aes(x = -122.4, y = 38.05, label = "San Pablo\nBay"), 
#            colour = "white", size = 3, check_overlap = T) +
#  geom_text(aes(x = -122.44, y = 37.75, label = "San Francisco"), 
#            colour = "grey40", size = 3, check_overlap = T) +
#  geom_segment(aes(x = -121.54, xend = -121.54, y = 37.46, yend = 37.49), 
#               arrow = arrow(length = unit(0.30, "cm"))) +
#  geom_text(aes(x = -121.54, y = 37.51, label = "N"), 
#            colour = "black", size = 4, check_overlap = T) +
  xlab(NULL) + 
  ylab(NULL) +
#  scalebar(x.min = -122.6, y.min = 37.48, x.max = -121.6, y.max = 37.98, 
#           dist = 10, dist_unit = "km", height = 0.02, st.dist = 0.03, st.size = 4,
#           transform = TRUE, model = "WGS84", location = "bottomright") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,-1,0,-1), "cm"),
        axis.text = element_text(size = 8, color = "black"))
dev.off()

# There are 5 Delaware River sites, but also a North Carolina and South Carolina site!
# Zoom in on Delaware River
coords_del_only <- input_filt_rare$map_loaded %>%
  group_by(Latitude, Longitude) %>%
  slice_head() %>%
  filter(Latitude > 37)
del_only <- get_stamenmap(bbox = c(left = min(coords_del_only$Longitude) - 0.25, 
                                   bottom = min(coords_del_only$Latitude) - 0.2, 
                                   right = max(coords_del_only$Longitude) + 0.1, 
                                   top = max(coords_del_only$Latitude) + 0.1),
                     zoom = 10, 
                     maptype = "terrain-background")
del_only_attributes <- attributes(del_only)
del_only_transparent <- matrix(adjustcolor(del_only, alpha.f = 0.4), nrow = nrow(del_only))
attributes(del_only_transparent) <- del_only_attributes
pdf("Figs/Map_Delaware.pdf", width = 6, height = 6)
ggmap(del_only_transparent, extent = "device") + # the base map
  geom_point(data = coords_del_only,
             aes(x = Longitude, y = Latitude), size = 4) +
  geom_text(aes(x = -75.2, y = 39.95, label = "Philadelphia"), 
            colour = "grey40", size = 3, check_overlap = T) +
  geom_text(aes(x = -75.56, y = 39.75, label = "Wilmington"), 
            colour = "grey40", size = 3, check_overlap = T) +
  geom_text(aes(x = -75.39, y = 39.33, label = "Delaware River"), 
            colour = "white", angle = 315, size = 4, fontface = "bold", check_overlap = T) +
  geom_segment(aes(x = -74.8, xend = -74.8, y = 39.25, yend = 39.28), 
               arrow = arrow(length = unit(0.30, "cm"))) +
  geom_text(aes(x = -74.8, y = 39.3, label = "N"), 
            colour = "black", size = 4, check_overlap = T) +
  xlab(NULL) + 
  ylab(NULL) +
  scalebar(x.min = min(coords_del_only$Longitude), 
           y.min = min(coords_del_only$Latitude)-0.15, 
           x.max = max(coords_del_only$Longitude)-0.05, 
           y.max = max(coords_del_only$Latitude), 
           dist = 10, dist_unit = "km", height = 0.02, st.dist = 0.03, st.size = 4,
           transform = TRUE, model = "WGS84", location = "bottomright") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,-1,0,-1), "cm"),
        axis.text = element_text(size = 8, color = "black"))
dev.off()



#### Alpha Diversity ####
#### _OTU rich ####
leveneTest(input_filt_rare$map_loaded$rich ~ input_filt_rare$map_loaded$Treatment)
# Variance homogeneous (p > 0.05)
m <- aov(input_filt_rare$map_loaded$rich ~ input_filt_rare$map_loaded$Treatment)
shapiro.test(m$residuals)
# Residuals not normally distributed (p < 0.05)
summary(m)

pdf("Figs/EastCoast_Richness.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Treatment, rich, colour = Experiment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.25) +
  labs(x = "Site", y = "Number of OTUs", colour = "Experiment") +
  facet_wrap(~ Experiment, ncol = 3, scales = "free") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = c(0.75, 0.2),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

#### _Shannon ####
leveneTest(input_filt_rare$map_loaded$shannon ~ input_filt_rare$map_loaded$Treatment)
# Variance homogeneous (p > 0.05)
m1 <- aov(input_filt_rare$map_loaded$shannon ~ input_filt_rare$map_loaded$Treatment)
shapiro.test(m1$residuals)
# Residuals not normally distributed (p < 0.05)
summary(m1)

pdf("Figs/EastCoast_Shannon.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Treatment, shannon, colour = Experiment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.25) +
  labs(x = "Treatment", y = "Shannon diversity", colour = "Experiment") +
  facet_wrap(~ Experiment, ncol = 3, scales = "free") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = c(0.75, 0.2),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### Beta Diversity ####
bc <- calc_dm(input_filt_rare$data_loaded)
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_filt_rare$map_loaded$Axis01 <- scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Experiment"), find_hull)
pdf("Figs/EastCoast_PCoA.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02, colour = Experiment)) +
  geom_polygon(data = micro.hulls, aes(colour = Experiment, fill = Experiment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5, aes(shape = Estuary)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Experiment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
dev.off()

# Stats
set.seed(1150)
input_filt_rare$map_loaded$Estuary <- as.factor(input_filt_rare$map_loaded$Estuary)
adonis2(bc ~ Estuary + Experiment + Depth, data = input_filt_rare$map_loaded)
anova(betadisper(bc, input_filt_rare$map_loaded$Estuary)) # Dispersion not homogeneous
anova(betadisper(bc, input_filt_rare$map_loaded$Experiment)) # Dispersion not homogeneous
anova(betadisper(bc, input_filt_rare$map_loaded$Depth)) # Dispersion not homogeneous



#### Taxa ####
# Prelim exploration but don't save anything. 
# Will redo with SF Bay data included and save figures

#### _Indicators ####
sim <- simper(t(input_filt_rare$data_loaded), 
              input_filt_rare$map_loaded$Experiment)
s <- summary(sim)
head(s$`Soil incubation_Soil Incubation`, n = 10)
head(s$`Soil_Soil Field plots`, n = 10)

# MULTIPATT (list ASVs associated with each group)
set.seed(1202)
mp <- multipatt(t(input_filt_rare$data_loaded), 
                input_filt_rare$map_loaded$Experiment, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp)

#### _Domain ####
tax_sum_domain <- summarize_taxonomy(input_filt_rare, level = 1, report_higher_tax = F)
plot_ts_heatmap(tax_sum_domain, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_domain,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative Abundance", fill = "Domain") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_domain, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')



#### _Phylum ####
tax_sum_phyla <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = F)
plot_ts_heatmap(tax_sum_phyla, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_phyla,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_phyla, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Look at archaeal phyla
tax_sum_phyla_ar <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = T)
tax_sum_phyla_ar <- tax_sum_phyla_ar[grep("Archaea", rownames(tax_sum_phyla_ar)),]
bars_ar <- plot_taxa_bars(tax_sum_phyla_ar,
                          input_filt_rare$map_loaded,
                          "Experiment",
                          num_taxa = 13,
                          data_only = TRUE)
nb.cols <- nrow(tax_sum_phyla_ar)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_ar, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

#### _Class ####
# Check sulfur reducers
tax_sum_phyla_su <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = T)
tax_sum_phyla_su <- tax_sum_phyla_su[grep("Desulfo", rownames(tax_sum_phyla_su)),]
bars_su <- plot_taxa_bars(tax_sum_phyla_su,
                          input_filt_rare$map_loaded,
                          "Experiment",
                          num_taxa = 13,
                          data_only = TRUE) %>%
  mutate(taxon = substring(taxon, 11))
nb.cols <- nrow(tax_sum_phyla_su)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_su, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))



#### _Family ####
tax_sum_families <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = FALSE)
plot_ts_heatmap(tax_sum_families, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars2 <- plot_taxa_bars(tax_sum_families,
                        input_filt_rare$map_loaded,
                        "Experiment",
                        num_taxa = 10,
                        data_only = TRUE)
ggplot(bars2, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative Abundance", fill = "Genus") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_families, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Look at methanogens
tax_sum_families_meth <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = F)
tax_sum_families_meth <- tax_sum_families_meth[grep("Methano", rownames(tax_sum_families_meth)),]
bars_meth <- plot_taxa_bars(tax_sum_families_meth,
                            input_filt_rare$map_loaded,
                            "Experiment",
                            num_taxa = nrow(tax_sum_families_meth),
                            data_only = TRUE)
nb.cols <- nrow(tax_sum_families_meth)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_meth, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

#### _Genus ####
tax_sum_genera <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = TRUE)
plot_ts_heatmap(tax_sum_genera, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_genera,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative Abundance", fill = "Genus") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_genera, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Guilds ####
tax_sum_guilds <- summarize_taxonomy(input_filt_rare, level = 9, report_higher_tax = F)
plot_ts_heatmap(tax_sum_guilds, 
                input_filt_rare$map_loaded, 
                0, 
                'Experiment',
                rev_taxa = T,
                remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))

bars <- plot_taxa_bars(tax_sum_guilds,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = nrow(tax_sum_guilds),
                       data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_guilds, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')



#### _Venn ####
phy <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = F)
cla <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = F)
ord <- summarize_taxonomy(input_filt_rare, level = 4, report_higher_tax = F)
fam <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = F)
gen <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = F)

input_phylum <- input_filt_rare
input_phylum$data_loaded <- phy
input_class <- input_filt_rare
input_class$data_loaded <- cla
input_order <- input_filt_rare
input_order$data_loaded <- ord
input_family <- input_filt_rare
input_family$data_loaded <- fam
input_genus <- input_filt_rare
input_genus$data_loaded <- gen

plot_venn_diagram(input_filt_rare,
                  "Estuary",
                  0.00000000000000001)

plot_grid(plot_venn_diagram(input_phylum, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_class, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_order, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_family, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_genus, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_filt_rare, "Estuary", 0.00000000000000001),
          labels = c("(a) Phylum", "(b) Class", "(c) Order", 
                     "(d) Family", "(e) Genus", "(f) OTU"))



#### ...................................... ####
#### East Coast Experiments ####
# Analyze each of the East Coast experiments separately
# Look at microbial responses to the different sample types/manipulations/time points



#### _South Carolina ####
input_filt <- readRDS("input_filt.rds")
sc <- filter_data(input_filt,
                  filter_cat = "Estuary",
                  keep_vals = "Waccamaw")
set.seed(530)
sc <- single_rarefy(sc, min(colSums(sc$data_loaded))) # 79282
sc$map_loaded <- sc$map_loaded %>%
  mutate(rich = specnumber(sc$data_loaded, MARGIN = 2),
         shannon = diversity(sc$data_loaded, index = "shannon", MARGIN = 2),
         TrtDepth = paste(sc$map_loaded$Treatment, sc$map_loaded$Depth, sep = ""),
         Depth = as.factor(Depth)) %>%
  mutate_if(is.character, as.factor)

#### __Alpha ####
leveneTest(sc$map_loaded$rich ~ sc$map_loaded$Treatment)
m <- aov(rich ~ Treatment + Depth, data = sc$map_loaded)
Anova(m, type = "III") # Treatment, not Depth
m <- aov(rich ~ Treatment, data = sc$map_loaded)
shapiro.test(m$residuals)
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(sc$map_loaded$rich)+(max(sc$map_loaded$rich)-min(sc$map_loaded$rich))/20)
leveneTest(sc$map_loaded$shannon ~ sc$map_loaded$Treatment)
m1 <- aov(shannon ~ Treatment + Depth, data = sc$map_loaded)
Anova(m1, type = "III") # Treatment, not Depth
m1 <- aov(shannon ~ Treatment, data = sc$map_loaded)
shapiro.test(m1$residuals)
summary(m1)
t1 <- emmeans(object = m1, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(sc$map_loaded$shannon)+(max(sc$map_loaded$shannon)-min(sc$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- sc$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/SC_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                       colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_x_discrete(labels = c("+Saltwater", "+Freshwater", "Control")) +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
sc_bc <- calc_dm(sc$data_loaded)
set.seed(1150)
adonis2(sc_bc ~ sc$map_loaded$Treatment + sc$map_loaded$Depth) # Both sig
anova(betadisper(sc_bc, sc$map_loaded$Treatment)) # Dispersion not homogeneous
anova(betadisper(sc_bc, sc$map_loaded$Depth)) # Dispersion homogeneous
sc_pcoa <- cmdscale(sc_bc, k = nrow(sc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(sc_pcoa)/sum(eigenvals(sc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(sc_pcoa)/sum(eigenvals(sc_pcoa)))[2]*100, digits = 1)
sc$map_loaded$Axis01 <- scores(sc_pcoa)[,1]
sc$map_loaded$Axis02 <- scores(sc_pcoa)[,2]
micro.hulls <- ddply(sc$map_loaded, c("TrtDepth"), find_hull)
g <- ggplot(sc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = TrtDepth, fill = TrtDepth),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = TrtDepth, shape = TrtDepth),
             show.legend = F) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("#440154FF", "#440154FF", "#21908CFF", "#21908CFF",
                      "#FDE725FF", "#FDE725FF")) +
  scale_fill_manual(values = c("#440154FF", "#440154FF", "#21908CFF", "#21908CFF",
                                 "#FDE725FF", "#FDE725FF")) +
  scale_shape_manual(values = c(16,17,16,17,16,17)) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, -5, 5, 5, "pt"))
leg <- get_legend(ggplot(sc$map_loaded, aes(Axis01, Axis02)) +
                    geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
                    theme_bw() +
                    scale_colour_viridis_d(labels = c("Control", "+Freshwater", "+Saltwater"),
                                           guide = guide_legend(reverse = T,
                                                                override.aes = list(shape = 15))) +
                    labs(shape = "Depth (m)"))
pdf("Figs/SC_PCoA.pdf", width = 6, height = 4)
plot_grid(g, leg, rel_widths = c(4, 1))
dev.off()

#### __Taxa ####
bar_text <- data.frame(group_by = c("Freshwater amended0.02", "Freshwater amended0.1"),
                       y = c(1.05, 1.05),
                       label = c("|-----2 cm-----|",
                                 "|-----10 cm-----|"))

sc_phyla <- summarize_taxonomy(sc, level = 2, report_higher_tax = F)
plot_ts_heatmap(sc_phyla, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsP <- plot_taxa_bars(sc_phyla, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/SC_Phyla.pdf", width = 7, height = 5)
ggplot(sc_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(sc_phyla, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_class <- summarize_taxonomy(sc, level = 3, report_higher_tax = F)
plot_ts_heatmap(sc_class, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsC <- plot_taxa_bars(sc_class, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sc_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sc_class, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_order <- summarize_taxonomy(sc, level = 4, report_higher_tax = F)
plot_ts_heatmap(sc_order, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsO <- plot_taxa_bars(sc_order, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sc_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sc_order, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_family <- summarize_taxonomy(sc, level = 5, report_higher_tax = F)
plot_ts_heatmap(sc_family, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsF <- plot_taxa_bars(sc_family, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sc_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sc_family, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_genus <- summarize_taxonomy(sc, level = 6, report_higher_tax = F)
plot_ts_heatmap(sc_genus, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsG <- plot_taxa_bars(sc_genus, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sc_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sc_genus, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_guilds <- summarize_taxonomy(sc, level = 9, report_higher_tax = F)
plot_ts_heatmap(sc_guilds, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsGu <- plot_taxa_bars(sc_guilds,
                       sc$map_loaded,
                       "TrtDepth",
                       num_taxa = 20,
                       data_only = TRUE) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
bar_textGu <- data.frame(group_by = c("Freshwater amended0.02", "Freshwater amended0.1"),
                       y = c(0.35, 0.35),
                       label = c("|-----2 cm-----|",
                                 "|-----10 cm-----|"))
pdf("Figs/SC_Guilds.pdf", width = 7, height = 5)
ggplot(sc_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_textGu,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  ylim(0, 0.35) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(sc_guilds, 
                            sc$map_loaded, 
                            type_header = 'TrtDepth', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### __ Simper ####
sc_sim <- simper(t(sc$data_loaded), 
                 sc$map_loaded$TrtDepth)
sc_s <- summary(sc_sim)
head(sc_s$`Control0.02_Saltwater amended0.02`)
head(sc_s$`Control0.1_Saltwater amended0.1`)
sc_df1 <- head(sc_s$`Control0.02_Saltwater amended0.02`, n = 20) %>%
  mutate(Comparison = "ASW/Control 2 cm",
         ASV = rownames(.))
sc_df2 <- head(sc_s$`Control0.1_Saltwater amended0.1`, n = 20) %>%
  mutate(Comparison = "ASW/Control 10 cm",
         ASV = rownames(.))
sc_simper_results <- rbind(sc_df1, sc_df2) %>%
  left_join(., sc$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(ava > avb, "Postive", "Negative")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanSalt" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanSalt, MeanControl, CumulativeContribution)
write_xlsx(sc_simper_results, 
           "simper_results_sc.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
sc_mp <- multipatt(t(sc$data_loaded), 
                   sc$map_loaded$TrtDepth, 
                   func = "r.g", 
                   control = how(nperm=999))
sc_mp_results <- sc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.Control0.02`, `s.Control0.1`, `s.Freshwater amended0.02`,
                                   `s.Freshwater amended0.1`, `s.Saltwater amended0.02`,
                                   `s.Saltwater amended0.1`)),
         q.value = qvalue(sc_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(sc_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.85) %>%
  filter(num_sites <= 2) %>%
  left_join(., sc$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$s.Control0.02[i] == 1) {
    sc_mp_results$Group[i] <- "Control_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$s.Control0.1[i] == 1) {
    sc_mp_results$Group[i] <- "Control_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Freshwater amended0.02`[i] == 1) {
    sc_mp_results$Group[i] <- "+Fresh_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Freshwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Fresh_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Saltwater amended0.02`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Saltwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 2 & sc_mp_results$`s.Saltwater amended0.02`[i] == 1 &
      sc_mp_results$`s.Saltwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_both"
  }
}
table(sc_mp_results$Group)
sc_asv <- summarize_taxonomy(sc, level = 8, report_higher_tax = F)
sc_asv_all <- data.frame("RelAbundance" = round(rowMeans(sc_asv) * 100, digits = 4)) %>%
  mutate(ASV = rownames(.))
sc_mp_corrs <- as.data.frame(sc_mp$str) %>%
  dplyr::select(1:length(levels(sc$map_loaded$TrtDepth)), 
                `Saltwater amended0.02+Saltwater amended0.1`) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% sc_mp_results$ASV) %>%
  set_names(c("Control_0.02", "Control_0.1", "+Fresh_0.02", "+Fresh_0.1", "+Salt_0.02", 
              "+Salt_0.1", "+Salt_both", "ASV"))
# Add corrs and taxonomy
sc_mp_results <- sc_mp_results %>%
  filter(Group == "+Salt_0.02" | Group == "+Salt_0.1" | Group == "+Salt_both") %>%
  left_join(., sc_asv_all, by = "ASV") %>%
  left_join(., sc_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(sc_mp_corrs)[1:7], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(sc_mp_corrs)[1:7], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
sc.hm.melted <- sc_mp_results %>%
  dplyr::select(taxon, names(sc_mp_corrs)[1:7]) %>%
  melt(., id.vars = c("taxon"))
sc.hm <- ggplot(data = sc.hm.melted, 
             aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(sc.hm.melted$taxon), labels = unique(sc.hm.melted$taxon),
                   limits = rev(levels(sc.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
sc.l <- get_legend(sc.hm)
sc.hm.clean <- sc.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
sc.bp.y <- ggplot(data = sc_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(sc_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/SC_Multipatt.pdf", width = 8, height = 5)
plot_grid(sc.hm.clean, sc.bp.y, sc.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, OTUs are too rare, not even 1%. So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
sc_tax <- sc$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
sc_genus <- summarize_taxonomy(sc, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
sc_mp <- multipatt(t(sc_genus), 
                   sc$map_loaded$TrtDepth, 
                   func = "r.g", 
                   control = how(nperm=999))
sc_mp_results <- sc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.Control0.02`, `s.Control0.1`, `s.Freshwater amended0.02`,
                                   `s.Freshwater amended0.1`, `s.Saltwater amended0.02`,
                                   `s.Saltwater amended0.1`)),
         q.value = qvalue(sc_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(sc_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.7) %>%
  filter(num_sites <= 2) %>%
  left_join(., sc_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$s.Control0.02[i] == 1) {
    sc_mp_results$Group[i] <- "Control_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$s.Control0.1[i] == 1) {
    sc_mp_results$Group[i] <- "Control_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Freshwater amended0.02`[i] == 1) {
    sc_mp_results$Group[i] <- "+Fresh_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Freshwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Fresh_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Saltwater amended0.02`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Saltwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 2 & sc_mp_results$`s.Saltwater amended0.02`[i] == 1 &
      sc_mp_results$`s.Saltwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_both"
  }
}
table(sc_mp_results$Group)
sc_genus_all <- data.frame("RelAbundance" = round(rowMeans(sc_genus)/min(colSums(sc$data_loaded)) * 100, digits = 4)) %>%
  mutate(Genus = rownames(.))
sc_mp_corrs <- as.data.frame(sc_mp$str) %>%
  dplyr::select(1:length(levels(sc$map_loaded$TrtDepth)), 
                `Saltwater amended0.02+Saltwater amended0.1`) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% sc_mp_results$Genus) %>%
  set_names(c("Control_0.02", "Control_0.1", "+Fresh_0.02", "+Fresh_0.1", "+Salt_0.02", 
              "+Salt_0.1", "+Salt_both", "Genus"))
# Add corrs and taxonomy
sc_mp_results <- sc_mp_results %>%
  filter(Group == "+Salt_0.02" | Group == "+Salt_0.1" | Group == "+Salt_both") %>%
  left_join(., sc_genus_all, by = "Genus") %>%
  left_join(., sc_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(sc_mp_corrs)[1:7], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(sc_mp_corrs)[1:7], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
sc.hm.melted <- sc_mp_results %>%
  dplyr::select(taxon, names(sc_mp_corrs)[1:7]) %>%
  melt(., id.vars = c("taxon"))
sc.hm <- ggplot(data = sc.hm.melted, 
                aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(sc.hm.melted$taxon), labels = unique(sc.hm.melted$taxon),
                   limits = rev(levels(sc.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
sc.l <- get_legend(sc.hm)
sc.hm.clean <- sc.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
sc.bp.y <- ggplot(data = sc_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(sc_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/SC_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(sc.hm.clean, sc.bp.y, sc.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### _North Carolina ####
input_filt <- readRDS("input_filt.rds")
nc <- filter_data(input_filt,
                  filter_cat = "Estuary",
                  keep_vals = "Alligator")
nc$map_loaded$sampleID <- rownames(nc$map_loaded)
nc <- filter_data(nc,
                  filter_cat = "sampleID",
                  filter_vals = c("TL_nw_d1_DI_ctrl_AF1", "TL_nw_d1_DI_ctrl_AF3", 
                                  "TL_nw_d1_DI_ctrl_AF4", "TL_nw_d1_ASW_noS_BF3",
                                  "TL_nw_d1_ASW_noS_BF4", "TL_nw_d1_ASW_noS_BF5"))
set.seed(530)
nc <- single_rarefy(nc, min(colSums(nc$data_loaded))) # 82350
nc$map_loaded <- nc$map_loaded %>%
  mutate(rich = specnumber(nc$data_loaded, MARGIN = 2),
         shannon = diversity(nc$data_loaded, index = "shannon", MARGIN = 2),
         TrtDepth = paste(nc$map_loaded$Treatment, nc$map_loaded$Depth, sep = ""),
         Depth = as.factor(Depth)) %>%
  mutate_if(is.character, as.factor)

#### __Alpha ####
leveneTest(nc$map_loaded$rich ~ nc$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment + Depth, data = nc$map_loaded)
Anova(m, type = "III") # Treatment and Depth
m <- aov(rich ~ Treatment, data = nc$map_loaded)
shapiro.test(m$residuals) # Not normal, but not too bad
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(nc$map_loaded$rich)+(max(nc$map_loaded$rich)-min(nc$map_loaded$rich))/20)
leveneTest(nc$map_loaded$shannon ~ nc$map_loaded$Treatment) # Homogeneous
m1 <- aov(shannon ~ Treatment + Depth, data = nc$map_loaded)
Anova(m1, type = "III") # Treatment and Depth
m1 <- aov(shannon ~ Treatment, data = nc$map_loaded)
shapiro.test(m1$residuals) # Normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(nc$map_loaded$shannon)+(max(nc$map_loaded$shannon)-min(nc$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- nc$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/NC_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                       colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_x_discrete(labels = c("+Salt +SO4", "+Salt", "+SO4", "Control", "Field")) +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
nc_bc <- calc_dm(nc$data_loaded)
set.seed(1150)
adonis2(nc_bc ~ nc$map_loaded$Treatment + nc$map_loaded$Depth) # Both sig
anova(betadisper(nc_bc, nc$map_loaded$Treatment)) # Dispersion not homogeneous
anova(betadisper(nc_bc, nc$map_loaded$Depth)) # Dispersion homogeneous
nc_pcoa <- cmdscale(nc_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
pdf("Figs/NC_PCoA.pdf", width = 7, height = 5)
ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (m)") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme_bw() +  
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))
dev.off()

#### __Taxa ####
nc_phyla <- summarize_taxonomy(nc, level = 2, report_higher_tax = F)
plot_ts_heatmap(nc_phyla, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsP <- plot_taxa_bars(nc_phyla, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/NC_Phyla.pdf", width = 7, height = 5)
ggplot(nc_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+Salt", "+Salt +SO4")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(nc_phyla, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_class <- summarize_taxonomy(nc, level = 3, report_higher_tax = F)
plot_ts_heatmap(nc_class, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsC <- plot_taxa_bars(nc_class, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(nc_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+Salt", "+Salt +SO4")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(nc_class, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_order <- summarize_taxonomy(nc, level = 4, report_higher_tax = F)
plot_ts_heatmap(nc_order, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsO <- plot_taxa_bars(nc_order, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(nc_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+Salt", "+Salt +SO4")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(nc_order, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_family <- summarize_taxonomy(nc, level = 5, report_higher_tax = F)
plot_ts_heatmap(nc_family, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsF <- plot_taxa_bars(nc_family, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(nc_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+Salt", "+Salt +SO4")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(nc_family, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_genus<- summarize_taxonomy(nc, level = 6, report_higher_tax = F)
plot_ts_heatmap(nc_genus, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsG <- plot_taxa_bars(nc_genus, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(nc_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+Salt", "+Salt +SO4")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(nc_genus, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_guilds <- summarize_taxonomy(nc, level = 9, report_higher_tax = F)
plot_ts_heatmap(nc_guilds, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsGu <- plot_taxa_bars(nc_guilds, nc$map_loaded, "Treatment", num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW")))
pdf("Figs/NC_Guilds.pdf", width = 7, height = 5)
ggplot(nc_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Treatment", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+Salt", "+Salt +SO4")) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
dev.off()
taxa_summary_by_sample_type(nc_guilds, 
                            nc$map_loaded, 
                            type_header = 'Treatment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### __Simper ####
nc_sim <- simper(t(nc$data_loaded), 
                 nc$map_loaded$TrtDepth)
nc_s <- summary(nc_sim)
nc_df1 <- head(nc_s$`5ppt ASW0.025_DI_ctrl0.025`, n = 20) %>%
  mutate(Comparison = "ASW/Control 2.5 cm",
         ASV = rownames(.))
nc_df2 <- head(nc_s$`5ppt ASW0.125_DI_ctrl0.125`, n = 20) %>%
  mutate(Comparison = "ASW/Control 12.5 cm",
         ASV = rownames(.))
nc_simper_results <- rbind(nc_df1, nc_df2) %>%
  left_join(., nc$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(ava > avb, "Postive", "Negative")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanSalt" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
         Species, OTU, MeanSalt, MeanControl, CumulativeContribution)
write_xlsx(nc_simper_results, 
           "simper_results_nc.xlsx",
           format_headers = F)


#### __Multipatt ####
set.seed(1202)
nc_mp <- multipatt(t(nc$data_loaded), 
                   nc$map_loaded$Treatment, 
                   func = "r.g", 
                   control = how(nperm=999))
multipatt_results <- nc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.5ppt ASW`, `s.DI_ctrl`, `s.Field Reference`,
                                   `s.SO4 amended`, `s.SW_noSO4`)),
         q.value = qvalue(nc_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(nc_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.65) %>%
  filter(num_sites <= 2) %>%
  left_join(., nc$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$`s.5ppt ASW`[i] == 1) {
    multipatt_results$Group[i] <- "+Salt +SO4"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$s.DI_ctrl[i] == 1) {
    multipatt_results$Group[i] <- "Control"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$`s.Field Reference`[i] == 1) {
    multipatt_results$Group[i] <- "Field"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$`s.SO4 amended`[i] == 1) {
    multipatt_results$Group[i] <- "+SO4"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$s.SW_noSO4[i] == 1) {
    multipatt_results$Group[i] <- "+Salt"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 2 & multipatt_results$`s.5ppt ASW`[i] == 1 &
      multipatt_results$s.SW_noSO4[i] == 1) {
    multipatt_results$Group[i] <- "+Salt, +Salt+SO4"
  }
}
table(multipatt_results$Group)
nc_asv <- summarize_taxonomy(nc, level = 8, report_higher_tax = F)
nc_asv_all <- data.frame("RelAbundance" = round(rowMeans(nc_asv) * 100, digits = 4)) %>%
  mutate(ASV = rownames(.))
nc_mp_corrs <- as.data.frame(nc_mp$str) %>%
  dplyr::select(1:5, 9) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% multipatt_results$ASV) %>%
  set_names(c("+Salt +SO4", "Control", "Field", "+SO4", "+Salt", "+Salt or +Salt+SO4", "ASV"))
# Add corrs and taxonomy
multipatt_results <- multipatt_results %>%
  filter(Group == "+Salt" | Group == "+Salt, +Salt+SO4" | Group == "+Salt+SO4" |
           Group == "+SO4") %>%
  left_join(., nc_asv_all, by = "ASV") %>%
  left_join(., nc_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, "+Salt +SO4", "Control", "Field", "+SO4", "+Salt",
                "+Salt or +Salt+SO4", "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", "+Salt +SO4", "Control", "Field", "+SO4", "+Salt",
              "+Salt+SO4/+Salt", "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
hm.melted <- multipatt_results %>%
  dplyr::select(taxon, Field, Control, "+SO4", "+Salt", "+Salt +SO4", "+Salt+SO4/+Salt") %>%
  melt(., id.vars = c("taxon"))
hm <- ggplot(data = hm.melted, 
             aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-0.8, 0.8)) +
  scale_x_discrete(breaks = unique(hm.melted$taxon), labels = unique(hm.melted$taxon),
                   limits = rev(levels(hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
l <- get_legend(hm)
hm.clean <- hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-3,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
bp.y <- ggplot(data = multipatt_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(multipatt_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-3,0))), 
        plot.margin = margin(c(0,-2,0,-5))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/NC_Multipatt.pdf", width = 8, height = 5)
plot_grid(hm.clean, bp.y, l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Do at genus level
nc_tax <- nc$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
nc_genus <- summarize_taxonomy(nc, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
nc_mp <- multipatt(t(nc_genus), 
                   nc$map_loaded$Treatment, 
                   func = "r.g", 
                   control = how(nperm=999))
nc_mp_results <- nc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.5ppt ASW`, `s.DI_ctrl`, `s.Field Reference`,
                                   `s.SO4 amended`, `s.SW_noSO4`)),
         q.value = qvalue(nc_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(nc_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.5) %>%
  filter(num_sites <= 2) %>%
  left_join(., nc_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$`s.5ppt ASW`[i] == 1) {
    nc_mp_results$Group[i] <- "+Salt +SO4"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$s.DI_ctrl[i] == 1) {
    nc_mp_results$Group[i] <- "Control"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$`s.Field Reference`[i] == 1) {
    nc_mp_results$Group[i] <- "Field"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$`s.SO4 amended`[i] == 1) {
    nc_mp_results$Group[i] <- "+SO4"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$s.SW_noSO4[i] == 1) {
    nc_mp_results$Group[i] <- "+Salt"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 2 & nc_mp_results$`s.5ppt ASW`[i] == 1 &
      nc_mp_results$s.SW_noSO4[i] == 1) {
    nc_mp_results$Group[i] <- "+Salt, +Salt+SO4"
  }
}
table(nc_mp_results$Group)
nc_genus_all <- data.frame("RelAbundance" = round(rowMeans(nc_genus)/min(colSums(nc$data_loaded)) * 100, digits = 4)) %>%
  mutate(Genus = rownames(.))
nc_mp_corrs <- as.data.frame(nc_mp$str) %>%
  dplyr::select(1:5, 9) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% nc_mp_results$Genus) %>%
  set_names(c("+Salt +SO4", "Control", "Field", "+SO4", "+Salt", "+Salt or +Salt+SO4", "Genus"))
# Add corrs and taxonomy
nc_mp_results <- nc_mp_results %>%
  filter(Group == "+Salt" | Group == "+Salt, +Salt+SO4" | Group == "+Salt +SO4" |
           Group == "+SO4") %>%
  left_join(., nc_genus_all, by = "Genus") %>%
  left_join(., nc_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, "+Salt +SO4", "Control", "Field", "+SO4", "+Salt",
                "+Salt or +Salt+SO4", "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", "+Salt +SO4", "Control", "Field", "+SO4", "+Salt",
              "+Salt+SO4/+Salt", "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
nc.hm.melted <- nc_mp_results %>%
  dplyr::select(taxon, Field, Control, "+SO4", "+Salt", "+Salt +SO4", "+Salt+SO4/+Salt") %>%
  melt(., id.vars = c("taxon"))
nc.hm <- ggplot(data = nc.hm.melted, 
                aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-0.8, 0.8)) +
  scale_x_discrete(breaks = unique(nc.hm.melted$taxon), labels = unique(nc.hm.melted$taxon),
                   limits = rev(levels(nc.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
nc.l <- get_legend(nc.hm)
nc.hm.clean <- nc.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
nc.bp.y <- ggplot(data = nc_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(nc_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/nc_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(nc.hm.clean, nc.bp.y, nc.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### _Delaware field ####
# Delaware field (3 sites, Weston et al. 2014)
# Tidal freshwater, oligohaline, mesohaline
input_filt <- readRDS("input_filt.rds")
defie <- filter_data(input_filt,
                     filter_cat = "Estuary",
                     keep_vals = "Delaware")
defie <- filter_data(defie,
                     filter_cat = "Experiment",
                     keep_vals = "Soil")
set.seed(530)
defie <- single_rarefy(defie, min(colSums(defie$data_loaded))) # 52468
defie$map_loaded <- defie$map_loaded %>%
  mutate(rich = specnumber(defie$data_loaded, MARGIN = 2),
         shannon = diversity(defie$data_loaded, index = "shannon", MARGIN = 2),
         Treatment = as.factor(Treatment),
         Salt = recode_factor(Treatment, "TFM1_source" = "Freshwater",
                              "TFM2_source" = "Freshwater",
                              "MesoHal_source" = "Mesohaline",
                              "OligoHal_source" = "Oligohaline")) %>%
  mutate(Salt = factor(Salt, levels = c("Freshwater", "Oligohaline", "Mesohaline"))) %>%
  unite("TrtDepth", c("Salt", "Depth"), sep = "", remove = F) %>%
  mutate(Depth = as.factor(Depth)) %>%
  mutate_if(is.character, as.factor)

#### __Alpha ####
leveneTest(defie$map_loaded$rich ~ defie$map_loaded$Salt) # Homogeneous
m <- aov(rich ~ Salt + Depth, data = defie$map_loaded)
Anova(m, type = "III") # Salt, not Depth
m <- aov(rich ~ Salt, data = defie$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(defie$map_loaded$rich)+(max(defie$map_loaded$rich)-min(defie$map_loaded$rich))/20)
leveneTest(defie$map_loaded$shannon ~ defie$map_loaded$Salt) # Homogeneous
m1 <- aov(shannon ~ Salt + Depth, data = defie$map_loaded)
Anova(m1, type = "III") # Salt, not Depth
m1 <- aov(shannon ~ Salt, data = defie$map_loaded)
shapiro.test(m1$residuals) # Normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(defie$map_loaded$shannon)+(max(defie$map_loaded$shannon)-min(defie$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- defie$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/DEfie_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value, 
                       colour = Salt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Salt, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
defie_bc <- calc_dm(defie$data_loaded)
set.seed(1150)
adonis2(defie_bc ~ defie$map_loaded$Salt + defie$map_loaded$Depth) # Salt sig
anova(betadisper(defie_bc, defie$map_loaded$Salt)) # Dispersion homogeneous
anova(betadisper(defie_bc, defie$map_loaded$Depth)) # Dispersion homogeneous
defie_pcoa <- cmdscale(defie_bc, k = nrow(defie$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(defie_pcoa)/sum(eigenvals(defie_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(defie_pcoa)/sum(eigenvals(defie_pcoa)))[2]*100, digits = 1)
defie$map_loaded$Axis01 <- scores(defie_pcoa)[,1]
defie$map_loaded$Axis02 <- scores(defie_pcoa)[,2]
micro.hulls <- ddply(defie$map_loaded, c("Salt"), find_hull)
pdf("Figs/DEfie_PCoA.pdf", width = 7, height = 5)
ggplot(defie$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

#### __Taxa ####
defie_phyla <- summarize_taxonomy(defie, level = 2, report_higher_tax = F)
plot_ts_heatmap(defie_phyla, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsP <- plot_taxa_bars(defie_phyla, defie$map_loaded, "Salt", num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/DEfie_Phyla.pdf", width = 7, height = 5)
ggplot(defie_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(defie_phyla, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_class <- summarize_taxonomy(defie, level = 3, report_higher_tax = F)
plot_ts_heatmap(defie_class, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsC <- plot_taxa_bars(defie_class, defie$map_loaded, "Salt", 
                                 num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(defie_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(defie_class, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_order <- summarize_taxonomy(defie, level = 4, report_higher_tax = F)
plot_ts_heatmap(defie_order, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsO <- plot_taxa_bars(defie_order, defie$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(defie_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(defie_order, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_family <- summarize_taxonomy(defie, level = 5, report_higher_tax = F)
plot_ts_heatmap(defie_family, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsF <- plot_taxa_bars(defie_family, defie$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(defie_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(defie_family, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_genus <- summarize_taxonomy(defie, level = 6, report_higher_tax = F)
plot_ts_heatmap(defie_genus, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsG <- plot_taxa_bars(defie_genus, defie$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(defie_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(defie_genus, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_guilds <- summarize_taxonomy(defie, level = 9, report_higher_tax = F)
plot_ts_heatmap(defie_guilds, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsGu <- plot_taxa_bars(defie_guilds, defie$map_loaded, "Salt", num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
pdf("Figs/DEfie_Guilds.pdf", width = 7, height = 5)
ggplot(defie_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Salt", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        axis.title.x = element_blank())
dev.off()
taxa_summary_by_sample_type(defie_guilds, 
                            defie$map_loaded, 
                            type_header = 'Salt', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### __Simper ####
defie_sim <- simper(t(defie$data_loaded), 
                 defie$map_loaded$Salt)
defie_s <- summary(defie_sim)
head(defie_s$Freshwater_Mesohaline)
defie_simper_results <- head(defie_s$Freshwater_Mesohaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline") %>%
  left_join(., defie$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Postive", "Negative")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanSalt" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanSalt, MeanControl, CumulativeContribution)
write_xlsx(defie_simper_results, 
           "simper_results_DEfie.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
defie_mp <- multipatt(t(defie$data_loaded), 
                   defie$map_loaded$Salt, 
                   func = "r.g", 
                   control = how(nperm=999))
# None with Q, use P
defie_mp_results <- defie_mp$sign %>%
  mutate(num_sites = rowSums(cbind(s.Freshwater, s.Oligohaline, s.Mesohaline)),
         q.value = qvalue(defie_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(defie_mp$sign)) %>%
  filter(p.value == 0.001) %>%
  filter(stat >= 0.9) %>%
  filter(num_sites <= 2) %>%
  left_join(., defie$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Freshwater[i] == 1) {
    defie_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Oligohaline[i] == 1) {
    defie_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Mesohaline[i] == 1) {
    defie_mp_results$Group[i] <- "Mesohaline"
  }
}
table(defie_mp_results$Group)
defie_asv <- summarize_taxonomy(defie, level = 8, report_higher_tax = F)
defie_asv_all <- data.frame("RelAbundance" = round(rowMeans(defie_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
defie_mp_corrs <- as.data.frame(defie_mp$str) %>%
  dplyr::select(1:length(levels(defie$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% defie_mp_results$ASV) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesolhaline", "ASV"))
# Add corrs and taxonomy
defie_mp_results <- defie_mp_results %>%
  filter(Group == "Mesohaline") %>%
  left_join(., defie_asv_all, by = "ASV") %>%
  left_join(., defie_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(defie_mp_corrs)[1:3], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(defie_mp_corrs)[1:3], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
defie.hm.melted <- defie_mp_results %>%
  dplyr::select(taxon, names(defie_mp_corrs)[1:3]) %>%
  melt(., id.vars = c("taxon"))
defie.hm <- ggplot(data = defie.hm.melted, 
                aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(defie.hm.melted$taxon), labels = unique(defie.hm.melted$taxon),
                   limits = rev(levels(defie.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
defie.l <- get_legend(defie.hm)
defie.hm.clean <- defie.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
defie.bp.y <- ggplot(data = defie_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(defie_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/DEfie_Multipatt.pdf", width = 8, height = 5)
plot_grid(defie.hm.clean, defie.bp.y, defie.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, many OTUs are too rare, not even 1%. So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
defie_tax <- defie$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
defie_genus <- summarize_taxonomy(defie, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
defie_mp <- multipatt(t(defie_genus), 
                   defie$map_loaded$Salt, 
                   func = "r.g", 
                   control = how(nperm=999))
defie_mp_results <- defie_mp$sign %>%
  mutate(num_sites = rowSums(cbind(s.Freshwater, s.Oligohaline, s.Mesohaline)),
         q.value = qvalue(defie_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(defie_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.7) %>%
  filter(num_sites <= 2) %>%
  left_join(., defie_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Freshwater[i] == 1) {
    defie_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Oligohaline[i] == 1) {
    defie_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Mesohaline[i] == 1) {
    defie_mp_results$Group[i] <- "Mesohaline"
  }
}
table(defie_mp_results$Group)
defie_genus_all <- data.frame("RelAbundance" = round(rowMeans(defie_genus)/min(colSums(defie$data_loaded)) * 100, digits = 4)) %>%
  mutate(Genus = rownames(.))
defie_mp_corrs <- as.data.frame(defie_mp$str) %>%
  dplyr::select(1:length(levels(defie$map_loaded$Salt))) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% defie_mp_results$Genus) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesohaline", "Genus"))
# Add corrs and taxonomy
defie_mp_results <- defie_mp_results %>%
  filter(Group == "Mesohaline") %>%
  left_join(., defie_genus_all, by = "Genus") %>%
  left_join(., defie_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(defie_mp_corrs)[1:3], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(defie_mp_corrs)[1:3], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
defie.hm.melted <- defie_mp_results %>%
  dplyr::select(taxon, names(defie_mp_corrs)[1:3]) %>%
  melt(., id.vars = c("taxon"))
defie.hm <- ggplot(data = defie.hm.melted, 
                aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(defie.hm.melted$taxon), labels = unique(defie.hm.melted$taxon),
                   limits = rev(levels(defie.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
defie.l <- get_legend(defie.hm)
defie.hm.clean <- defie.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
defie.bp.y <- ggplot(data = defie_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(defie_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/DEfie_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(defie.hm.clean, defie.bp.y, defie.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### _Delaware inc ####
input_filt_rare <- readRDS("input_filt.rds")
deinc <- filter_data(input_filt_rare,
                    filter_cat = "Estuary",
                    keep_vals = "Delaware")
deinc <- filter_data(deinc,
                    filter_cat = "Experiment",
                    keep_vals = "Soil incubation") # 28 samples

# Diagnose outliers and errors
deinc_bc <- calc_dm(deinc$data_loaded)
deinc_pcoa <- cmdscale(deinc_bc, k = nrow(deinc$map_loaded) - 1, eig = T)
eigenvals(deinc_pcoa)/sum(eigenvals(deinc_pcoa)) # 30.3, 17.1 % variation explained
deinc$map_loaded$Axis01 <- scores(deinc_pcoa)[,1]
deinc$map_loaded$Axis02 <- scores(deinc_pcoa)[,2]
micro.hulls <- ddply(deinc$map_loaded, c("Treatment"), find_hull)
ggplotly(ggplot(deinc$map_loaded, aes(Axis01, Axis02)) +
           geom_polygon(data = micro.hulls, 
                        aes(colour = Treatment, fill = Treatment),
                        alpha = 0.1, show.legend = F) +
           geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
           labs(x = "PC1: 30.3%", 
                y = "PC2: 17.1%",
                shape = "Depth (m)") +
           scale_colour_viridis_d() +
           scale_fill_viridis_d() +
           theme_bw() +  
           theme(legend.position = c(0,0),
                 legend.justification = c(0,0),
                 legend.background = element_blank(),
                 axis.title = element_text(face = "bold", size = 12), 
                 axis.text = element_text(size = 10),
                 plot.margin = margin(5, 5, 5, 5, "pt")))

# Looks like there are two outliers and potentially 2 with depths mixed up!
# Restart
input_filt_rare <- readRDS("input_filt.rds")
deinc <- filter_data(input_filt_rare,
                     filter_cat = "Estuary",
                     keep_vals = "Delaware")
deinc <- filter_data(deinc,
                     filter_cat = "Experiment",
                     keep_vals = "Soil incubation") # 28 samples
deinc$map_loaded$sampleID <- rownames(deinc$map_loaded)
deinc <- filter_data(deinc, 
                     filter_cat = "sampleID", 
                     filter_vals = c("TS_FW_d1_12_2", "TS_FW_d2_12_2"))
for (i in 1:nrow(deinc$map_loaded)) {
  if (deinc$map_loaded$sampleID[i] == "TS_FW_d2_04_1") {
    deinc$map_loaded$Depth[i] <- "0.02"
  }
}
for (i in 1:nrow(deinc$map_loaded)) {
  if (deinc$map_loaded$sampleID[i] == "TS_FW_d1_04_2") {
    deinc$map_loaded$Depth[i] <- "0.12"
  }
}

set.seed(530)
deinc <- single_rarefy(deinc, min(colSums(deinc$data_loaded))) # 31264
deinc$map_loaded <- deinc$map_loaded %>%
  mutate(rich = specnumber(deinc$data_loaded, MARGIN = 2),
         shannon = diversity(deinc$data_loaded, index = "shannon", MARGIN = 2),
         Depth = as.factor(Depth),
         sampleID = rownames(.)) %>%
  mutate(Treatment = gsub("5 ppt ASW", "ASW", Treatment)) %>%
  mutate(Treatment = gsub("wk", "", Treatment)) %>%
  mutate(Treatment = gsub("Freshwater", "Fresh", Treatment)) %>%
  separate(Treatment, into = c("Salt", "Time"), sep = " ", remove = F) %>%
  mutate(Time = replace_na(Time, 0)) %>%
  mutate(Time = as.integer(Time))
for (i in 1:nrow(deinc$map_loaded)) {
  if (deinc$map_loaded$Time[i] == 0) {
    deinc$map_loaded$Treatment[i] <- "Initial"
    deinc$map_loaded$Salt[i] <- "Initial"
  }
}
deinc$map_loaded <- deinc$map_loaded %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Treatment = factor(Treatment,
                          levels = c("Initial", "Fresh 4", "Fresh 7",
                                     "Fresh 12", "ASW 4", "ASW 7", "ASW 12")),
         TrtDepth = paste(deinc$map_loaded$Treatment, deinc$map_loaded$Depth, sep = ""))

#### __Alpha ####
leveneTest(deinc$map_loaded$rich ~ deinc$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment + Depth, data = deinc$map_loaded)
Anova(m, type = "III") # Neither
m <- aov(rich ~ Treatment, data = deinc$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(deinc$map_loaded$rich)+(max(deinc$map_loaded$rich)-min(deinc$map_loaded$rich))/20)
leveneTest(deinc$map_loaded$shannon ~ deinc$map_loaded$Treatment) # Not homogeneous
m1 <- aov(shannon ~ Treatment + Depth, data = deinc$map_loaded)
Anova(m1, type = "III") # Neither
m1 <- aov(shannon ~ Treatment, data = deinc$map_loaded)
shapiro.test(m1$residuals) # Normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(deinc$map_loaded$shannon)+(max(deinc$map_loaded$shannon)-min(deinc$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- deinc$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/DEinc_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(Treatment, value, 
                       colour = Salt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
# From here on just look at Initial and week 12
deinc <- filter_data(deinc, 
                     filter_cat = "Time",
                     keep_vals = c("0", "12"))
deinc_bc <- calc_dm(deinc$data_loaded)
deinc_pcoa <- cmdscale(deinc_bc, k = nrow(deinc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(deinc_pcoa)/sum(eigenvals(deinc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(deinc_pcoa)/sum(eigenvals(deinc_pcoa)))[2]*100, digits = 1)
deinc$map_loaded$Axis01 <- scores(deinc_pcoa)[,1]
deinc$map_loaded$Axis02 <- scores(deinc_pcoa)[,2]
micro.hulls <- ddply(deinc$map_loaded, c("Salt"), find_hull)
pdf("Figs/DEinc_PCoA.pdf", width = 7, height = 5)
ggplot(deinc$map_loaded, aes(Axis01, Axis02)) +
           geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Depth)) +
           labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
                y = paste("PC2: ", pcoaA2, "%", sep = ""),
                shape = "Depth (m)") +
           scale_colour_viridis_d() +
           scale_fill_viridis_d() +
           theme_bw() +  
           theme(legend.position = "right",
                 axis.title = element_text(face = "bold", size = 12), 
                 axis.text = element_text(size = 10),
                 plot.margin = margin(5, 5, 5, 5, "pt"))
dev.off()

set.seed(1150)
adonis2(deinc_bc ~ Salt + Depth, data = deinc$map_loaded) # Depth
anova(betadisper(deinc_bc, deinc$map_loaded$Treatment)) # Dispersion homogeneous
anova(betadisper(deinc_bc, deinc$map_loaded$Depth)) # Dispersion homogeneous
anova(betadisper(deinc_bc, deinc$map_loaded$Time)) # Dispersion homogeneous

#### __Taxa ####
bar_text <- data.frame(group_by = c("Fresh 120.02", "Fresh 120.12"),
                       y = c(1.05, 1.05),
                       label = c("|-----2 cm-----|",
                                 "|-----12 cm-----|"))
deinc_phyla <- summarize_taxonomy(deinc, level = 2, report_higher_tax = F)
plot_ts_heatmap(deinc_phyla, deinc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
plot_taxa_bars(deinc_phyla, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = F) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
deinc_barsP <- plot_taxa_bars(deinc_phyla, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/DEinc_Phyla.pdf", width = 7, height = 5)
ggplot(deinc_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(deinc_phyla, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_class <- summarize_taxonomy(deinc, level = 3, report_higher_tax = F)
plot_ts_heatmap(deinc_class, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsC <- plot_taxa_bars(deinc_class, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(deinc_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(deinc_class, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_order <- summarize_taxonomy(deinc, level = 4, report_higher_tax = F)
plot_ts_heatmap(deinc_order, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsO <- plot_taxa_bars(deinc_order, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(deinc_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(deinc_order, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_family <- summarize_taxonomy(deinc, level = 5, report_higher_tax = F)
plot_ts_heatmap(deinc_family, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsF <- plot_taxa_bars(deinc_family, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(deinc_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(deinc_family, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_genus<- summarize_taxonomy(deinc, level = 6, report_higher_tax = F)
plot_ts_heatmap(deinc_genus, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsG <- plot_taxa_bars(deinc_genus, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(deinc_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(deinc_genus, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_guilds <- summarize_taxonomy(deinc, level = 9, report_higher_tax = F)
plot_ts_heatmap(deinc_guilds, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsGu <- plot_taxa_bars(deinc_guilds, deinc$map_loaded, "TrtDepth", num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12")))
bar_textGu <- data.frame(group_by = c("Fresh 120.02", "Fresh 120.12"),
                       y = c(0.25, 0.25),
                       label = c("|-----2 cm-----|",
                                 "|-----12 cm-----|"))
pdf("Figs/DEinc_Guilds.pdf", width = 7, height = 5)
ggplot(deinc_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_textGu,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "TrtDepth", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(deinc_guilds, 
                            deinc$map_loaded, 
                            type_header = 'Treatment', 
                            filter_level = 0.01, 
                            test_type = 'KW')


#### __Simper ####
deinc_sim <- simper(t(deinc$data_loaded), 
                    deinc$map_loaded$TrtDepth)
deinc_s <- summary(deinc_sim)
head(deinc_s$`Fresh 120.02_ASW 120.02`)
head(deinc_s$`Fresh 120.12_ASW 120.12`)
deinc_df1 <- head(deinc_s$`Fresh 120.02_ASW 120.02`, n = 20) %>%
  mutate(Comparison = "ASW/Control 2 cm",
         ASV = rownames(.))
deinc_df2 <- head(deinc_s$`Fresh 120.12_ASW 120.12`, n = 20) %>%
  mutate(Comparison = "ASW/Control 12 cm",
         ASV = rownames(.))
deinc_simper_results <- rbind(deinc_df1, deinc_df2) %>%
  left_join(., deinc$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(ava > avb, "Postive", "Negative")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanSalt" = "ava",
           "Meadeincontrol" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanSalt, Meadeincontrol, CumulativeContribution)
write_xlsx(deinc_simper_results, 
           "simper_results_DEinc.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
deinc_mp <- multipatt(t(deinc$data_loaded), 
                   deinc$map_loaded$TrtDepth, 
                   func = "r.g", 
                   control = how(nperm=999))
# Note nothing with significant q value, no strong correlations
deinc_mp_results <- deinc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.ASW 120.02`, `s.ASW 120.12`, `s.Fresh 120.02`,
                                   `s.Fresh 120.12`, `s.Initial0.02`, `s.Initial0.12`)),
         q.value = qvalue(deinc_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(deinc_mp$sign)) %>%
  filter(p.value < 0.05) %>%
  filter(stat >= 0.3) %>%
  filter(num_sites < 2) %>%
  left_join(., deinc$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$`s.ASW 120.02`[i] == 1) {
    deinc_mp_results$Group[i] <- "ASW 2 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$`s.ASW 120.12`[i] == 1) {
    deinc_mp_results$Group[i] <- "ASW 12 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$`s.Fresh 120.02`[i] == 1) {
    deinc_mp_results$Group[i] <- "Fresh 2 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$`s.Fresh 120.12`[i] == 1) {
    deinc_mp_results$Group[i] <- "Fresh 12 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$s.Initial0.02[i] == 1) {
    deinc_mp_results$Group[i] <- "Initial 2 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$s.Initial0.12[i] == 1) {
    deinc_mp_results$Group[i] <- "Initial 12 cm"
  }
}
table(deinc_mp_results$Group)
deinc_asv <- summarize_taxonomy(deinc, level = 8, report_higher_tax = F)
deinc_asv_all <- data.frame("RelAbundance" = round(rowMeans(deinc_asv) * 100, digits = 4)) %>%
  mutate(ASV = rownames(.))
deinc_mp_corrs <- as.data.frame(deinc_mp$str) %>%
  dplyr::select(1:6) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% deinc_mp_results$ASV) %>%
  set_names(c("ASW 2 cm", "ASW 12 cm", "Fresh 2 cm", "Fresh 12 cm", "Initial 2 cm", 
              "Initial 12 cm", "ASV"))
# Add corrs and taxonomy
deinc_mp_results <- deinc_mp_results %>%
  filter(Group == "ASW 2 cm" | Group == "ASW 12 cm") %>%
  left_join(., deinc_asv_all, by = "ASV") %>%
  left_join(., deinc_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, "ASW 2 cm", "ASW 12 cm", 
                "Fresh 2 cm", "Fresh 12 cm", "Initial 2 cm", "Initial 12 cm", "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", "ASW 2 cm", "ASW 12 cm", 
              "Fresh 2 cm", "Fresh 12 cm", "Initial 2 cm", "Initial 12 cm", "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
hm.melted <- deinc_mp_results %>%
  dplyr::select(taxon, "ASW 2 cm", "ASW 12 cm", "Fresh 2 cm", "Fresh 12 cm", "Initial 2 cm", 
                "Initial 12 cm") %>%
  melt(., id.vars = c("taxon"))
hm <- ggplot(data = hm.melted, 
             aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-0.8, 0.8)) +
  scale_x_discrete(breaks = unique(hm.melted$taxon), labels = unique(hm.melted$taxon),
                   limits = rev(levels(hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
l <- get_legend(hm)
hm.clean <- hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-3,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
bp.y <- ggplot(data = deinc_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(deinc_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-3,0))), 
        plot.margin = margin(c(0,-2,0,-5))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/DEinc_Multipatt.pdf", width = 8, height = 5)
plot_grid(hm.clean, bp.y, l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### _Delaware transplant ####
input_filt <- readRDS("input_filt.rds")
detra <- filter_data(input_filt,
                     filter_cat = "Estuary",
                     keep_vals = "Delaware")
detra <- filter_data(detra,
                     filter_cat = "Experiment",
                     keep_vals = c("Soil mesocosm"))
set.seed(530)
detra <- single_rarefy(detra, min(colSums(detra$data_loaded))) # 26450
detra$map_loaded <- detra$map_loaded %>%
  mutate(rich = specnumber(detra$data_loaded, MARGIN = 2),
         shannon = diversity(detra$data_loaded, index = "shannon", MARGIN = 2),
         Treatment = as.factor(Treatment)) %>%
  unite("TrtDepth", c("Treatment", "Depth"), sep = "", remove = F) %>%
  mutate(Depth = as.factor(Depth)) %>%
  mutate_if(is.character, as.factor)

#### __Alpha ####
leveneTest(detra$map_loaded$rich ~ detra$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment + Depth, data = detra$map_loaded)
Anova(m, type = "III") # Treatment and Depth
m <- aov(rich ~ Treatment, data = detra$map_loaded)
shapiro.test(m$residuals) # Not normal
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(detra$map_loaded$rich)+(max(detra$map_loaded$rich)-min(detra$map_loaded$rich))/20)
leveneTest(detra$map_loaded$shannon ~ detra$map_loaded$Treatment) # Homogeneous
m1 <- aov(shannon ~ Treatment + Depth, data = detra$map_loaded)
Anova(m1, type = "III") # Depth
m1 <- aov(shannon ~ Treatment, data = detra$map_loaded)
shapiro.test(m1$residuals)
summary(m1)
t1 <- emmeans(object = m1, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(detra$map_loaded$shannon)+(max(detra$map_loaded$shannon)-min(detra$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- detra$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/DEtra_Alpha.pdf", width = 7, height = 4)
ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                       colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
detra_bc <- calc_dm(detra$data_loaded)
set.seed(1150)
adonis2(detra_bc ~ detra$map_loaded$Treatment + detra$map_loaded$Depth) # Both sig
anova(betadisper(detra_bc, detra$map_loaded$Treatment)) # Dispersion homogeneous
anova(betadisper(detra_bc, detra$map_loaded$Depth)) # Dispersion homogeneous
detra_pcoa <- cmdscale(detra_bc, k = nrow(detra$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(detra_pcoa)/sum(eigenvals(detra_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(detra_pcoa)/sum(eigenvals(detra_pcoa)))[2]*100, digits = 1)
detra$map_loaded$Axis01 <- scores(detra_pcoa)[,1]
detra$map_loaded$Axis02 <- scores(detra_pcoa)[,2]
micro.hulls <- ddply(detra$map_loaded, c("Treatment"), find_hull)
pdf("Figs/DEtra_PCoA.pdf", width = 7, height = 5)
ggplot(detra$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

#### __Taxa ####
detra_phyla <- summarize_taxonomy(detra, level = 2, report_higher_tax = F)
plot_ts_heatmap(detra_phyla, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsP <- plot_taxa_bars(detra_phyla, detra$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/DEtra_Phyla.pdf", width = 7, height = 5)
ggplot(detra_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(detra_phyla, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_class <- summarize_taxonomy(detra, level = 3, report_higher_tax = F)
plot_ts_heatmap(detra_class, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsC <- plot_taxa_bars(detra_class, detra$map_loaded, "Treatment", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(detra_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(detra_class, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_order <- summarize_taxonomy(detra, level = 4, report_higher_tax = F)
plot_ts_heatmap(detra_order, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsO <- plot_taxa_bars(detra_order, detra$map_loaded, "Treatment", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(detra_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(detra_order, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_family <- summarize_taxonomy(detra, level = 5, report_higher_tax = F)
plot_ts_heatmap(detra_family, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsF <- plot_taxa_bars(detra_family, detra$map_loaded, "Treatment", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(detra_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(detra_family, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_genus <- summarize_taxonomy(detra, level = 6, report_higher_tax = F)
plot_ts_heatmap(detra_genus, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsG <- plot_taxa_bars(detra_genus, detra$map_loaded, "Treatment", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(detra_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(detra_genus, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_guilds <- summarize_taxonomy(detra, level = 9, report_higher_tax = F)
plot_ts_heatmap(detra_guilds, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsGu <- plot_taxa_bars(detra_guilds, detra$map_loaded, "Treatment", num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
pdf("Figs/DEtra_Guilds.pdf", width = 7, height = 6)
ggplot(detra_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Treatment", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(detra_guilds, 
                            detra$map_loaded, 
                            type_header = 'Treatment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### __ Simper ####
detra_sim <- simper(t(detra$data_loaded), 
                    detra$map_loaded$Treatment)
detra_s <- summary(detra_sim)
detra_df1 <- head(detra_s$`TFM1@MesoHal_TFM1`, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline")
detra_df2 <- head(detra_s$`TFM1@MesoHal-40cm_TFM1-40cm`, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline 40 cm")
detra_simper_results <- rbind(detra_df1, detra_df2) %>%
  left_join(., detra$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(TreatmentResponse = ifelse(avb > ava, "Postive", "Negative")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanTreatment" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, TreatmentResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanTreatment, MeanControl, CumulativeContribution)
write_xlsx(detra_simper_results, 
           "simper_results_DEtra.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
detra_mp <- multipatt(t(detra$data_loaded), 
                      detra$map_loaded$Treatment, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 2)
# None with Q, use P
detra_mp_results <- detra_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.TFM1@MesoHal`, `s.TFM1@MesoHal-40cm`)),
         q.value = qvalue(detra_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(detra_mp$sign)) %>%
  filter(p.value == 0.001) %>%
  filter(stat >= 0.4) %>%
  # filter(num_sites <= 2) %>%
  left_join(., detra$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(detra_mp_results)) {
  if (detra_mp_results$num_sites[i] == 1 & detra_mp_results$`s.TFM1@MesoHal`[i] == 1) {
    detra_mp_results$Group[i] <- "TFM1@MesoHal"
  }
}
for (i in 1:nrow(detra_mp_results)) {
  if (detra_mp_results$num_sites[i] == 1 & detra_mp_results$`s.TFM1@MesoHal-40cm`[i] == 1) {
    detra_mp_results$Group[i] <- "TFM1@MesoHal-40cm"
  }
}
table(detra_mp_results$Group)
detra_asv <- summarize_taxonomy(detra, level = 8, report_higher_tax = F)
detra_asv_all <- data.frame("RelAbundance" = round(rowMeans(detra_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
detra_mp_corrs <- as.data.frame(detra_mp$str) %>%
  dplyr::select(`TFM1@MesoHal`, `TFM1@MesoHal-40cm`) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% detra_mp_results$ASV) %>%
  set_names(c("TFM1@MesoHal", "TFM1@MesoHal-40cm", "ASV"))
# Add corrs and taxonomy
detra_mp_results <- detra_mp_results %>%
  filter(Group == "TFM1@MesoHal" | Group == "TFM1@MesoHal-40cm") %>%
  left_join(., detra_asv_all, by = "ASV") %>%
  left_join(., detra_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(detra_mp_corrs)[1:2], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(detra_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
detra.hm.melted <- detra_mp_results %>%
  dplyr::select(taxon, names(detra_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
detra.hm <- ggplot(data = detra.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(detra.hm.melted$taxon), labels = unique(detra.hm.melted$taxon),
                   limits = rev(levels(detra.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
detra.l <- get_legend(detra.hm)
detra.hm.clean <- detra.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
detra.bp.y <- ggplot(data = detra_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(detra_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/DEtra_Multipatt.pdf", width = 8, height = 5)
plot_grid(detra.hm.clean, detra.bp.y, detra.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, OTUs are too rare, not even 1%. So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
detra_tax <- detra$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
detra_genus <- summarize_taxonomy(detra, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
detra_mp <- multipatt(t(detra_genus), 
                      detra$map_loaded$Treatment, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 2)
detra_mp_results <- detra_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.TFM1@MesoHal`, `s.TFM1@MesoHal-40cm`)),
         q.value = qvalue(detra_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(detra_mp$sign)) %>%
  filter(p.value < 0.01) %>%
  filter(stat >= 0.5) %>%
  filter(num_sites <= 2) %>%
  left_join(., detra_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(detra_mp_results)) {
  if (detra_mp_results$num_sites[i] == 1 & detra_mp_results$`s.TFM1@MesoHal`[i] == 1) {
    detra_mp_results$Group[i] <- "TFM1@MesoHal"
  }
}
for (i in 1:nrow(detra_mp_results)) {
  if (detra_mp_results$num_sites[i] == 1 & detra_mp_results$`s.TFM1@MesoHal-40cm`[i] == 1) {
    detra_mp_results$Group[i] <- "TFM1@MesoHal-40cm"
  }
}
table(detra_mp_results$Group)
detra_genus_all <- data.frame("RelAbundance" = round(rowMeans(detra_genus)/min(colSums(detra$data_loaded)) * 100, digits = 4)) %>%
  mutate(Genus = rownames(.))
detra_mp_corrs <- as.data.frame(detra_mp$str) %>%
  dplyr::select(`TFM1@MesoHal`, `TFM1@MesoHal-40cm`) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% detra_mp_results$Genus) %>%
  set_names(c("TFM1@MesoHal", "TFM1@MesoHal-40cm", "Genus"))
# Add corrs and taxonomy
detra_mp_results <- detra_mp_results %>%
  filter(Group == "TFM1@MesoHal" | Group == "TFM1@MesoHal-40cm") %>%
  left_join(., detra_genus_all, by = "Genus") %>%
  left_join(., detra_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(detra_mp_corrs)[1:2], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(detra_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
detra.hm.melted <- detra_mp_results %>%
  dplyr::select(taxon, names(detra_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
detra.hm <- ggplot(data = detra.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(detra.hm.melted$taxon), labels = unique(detra.hm.melted$taxon),
                   limits = rev(levels(detra.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
detra.l <- get_legend(detra.hm)
detra.hm.clean <- detra.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
detra.bp.y <- ggplot(data = detra_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(detra_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/DEtra_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(detra.hm.clean, detra.bp.y, detra.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### ...................................... ####
#### Comparison Overview ####
input_filt_rare <- readRDS("input_filt_rare_comb.rds")
input_filt_rare$map_loaded$Estuary <- factor(input_filt_rare$map_loaded$Estuary,
                                             levels = c("Waccamaw", "Alligator",
                                                        "Delaware", "SF"))
input_filt_rare_abund <- filter_taxa_from_input(input_filt_rare,
                                                filter_thresh = 0.05) # 94802 taxa removed



#### Alpha Diversity ####
#### _OTU rich ####
leveneTest(input_filt_rare$map_loaded$rich ~ input_filt_rare$map_loaded$Estuary)
# Variance not homogeneous (p < 0.05)
m <- aov(input_filt_rare$map_loaded$rich ~ input_filt_rare$map_loaded$Estuary)
shapiro.test(m$residuals)
# Residuals not normally distributed (p < 0.05)
summary(m)

pdf("Figs/CombRichness.pdf", width = 8, height = 4)
ggplot(input_filt_rare$map_loaded, aes(reorder(Info, rich, mean), rich, colour = Estuary)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.75, width = 0.2) +
  labs(x = "Site", y = "Number of OTUs", colour = "Estuary") +
  scale_colour_viridis_d() +
  facet_grid(~ Estuary, scales = "free_x", space = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 6.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.25, "cm"),
        panel.spacing.x = unit(0.05, "cm"))
dev.off()



#### _Shannon ####
leveneTest(input_filt_rare$map_loaded$shannon ~ input_filt_rare$map_loaded$Estuary)
# Variance homogeneous (p > 0.05)
m1 <- aov(input_filt_rare$map_loaded$shannon ~ input_filt_rare$map_loaded$Estuary)
shapiro.test(m1$residuals)
# Residuals not normally distributed (p < 0.05)
summary(m1)

pdf("Figs/CombShannon.pdf", width = 8, height = 4)
ggplot(input_filt_rare$map_loaded, aes(reorder(Info, shannon, mean), shannon, colour = Estuary)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.75, width = 0.2) +
  labs(x = "Site", y = "Shannon diversity", colour = "Estuary") +
  scale_colour_viridis_d() +
  facet_grid(~ Estuary, scales = "free_x", space = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 6.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.6, "cm"),
        panel.spacing.x = unit(0.05, "cm"))
dev.off()



#### Beta Diversity ####
bc <- calc_dm(input_filt_rare$data_loaded)
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_filt_rare$map_loaded$Axis01 <- scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- scores(pcoa)[,2]

# Explore the broad estuary comparison
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Estuary"), find_hull)
pdf("Figs/CombPCoA.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02, colour = Estuary)) +
  geom_polygon(data = micro.hulls, aes(colour = Estuary, fill = Estuary),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Estuary") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
dev.off()

# Interactive, with sample or site info
g1 <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02, shape = Estuary, colour = Info)) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Estuary",
       colour = "Info") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
g1
ggplotly(g1) # Export, Save as Web Page, CombPCoA.html

# Try with filtering out rare taxa. Use the 0.05% cutoff of Wyatt
# Overall pretty similar
bc_abund <- calc_dm(input_filt_rare_abund$data_loaded)
pcoa_abund <- cmdscale(bc_abund, k = nrow(input_filt_rare_abund$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_abund)/sum(eigenvals(pcoa_abund)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_abund)/sum(eigenvals(pcoa_abund)))[2]*100, digits = 1)
input_filt_rare_abund$map_loaded$Axis01 <- scores(pcoa_abund)[,1]
input_filt_rare_abund$map_loaded$Axis02 <- scores(pcoa_abund)[,2]
micro.hulls <- ddply(input_filt_rare_abund$map_loaded, c("Estuary"), find_hull)
pdf("Figs/CombPCoA_0.05.pdf", width = 7, height = 5)
ggplot(input_filt_rare_abund$map_loaded, aes(Axis01, Axis02, colour = Estuary)) +
  geom_polygon(data = micro.hulls, aes(colour = Estuary, fill = Estuary),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Estuary") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
dev.off()



# Stats
set.seed(1150)
adonis2(bc ~ input_filt_rare$map_loaded$Estuary)
set.seed(1150)
adonis2(bc_abund ~ input_filt_rare_abund$map_loaded$Estuary)

anova(betadisper(bc, input_filt_rare$map_loaded$Estuary)) # Dispersion not homogeneous
anova(betadisper(bc_abund, input_filt_rare_abund$map_loaded$Estuary)) # Dispersion not homogeneous


#### Taxa ####
#### _Indicators ####
sim <- simper(t(input_filt_rare_abund$data_loaded), 
              input_filt_rare_abund$map_loaded$Estuary)
s <- summary(sim)
s1 <- as.data.frame(head(s$Delaware_Waccamaw, n = 10)) %>%
  mutate(Comparison = "Delaware_Waccamaw",
         ASV = rownames(.))
s2 <- as.data.frame(head(s$Delaware_Alligator, n = 10)) %>%
  mutate(Comparison = "Delaware_Alligator",
         ASV = rownames(.))
s3 <- as.data.frame(head(s$Delaware_SF, n = 10)) %>%
  mutate(Comparison = "Delaware_SF",
         ASV = rownames(.))
s4 <- as.data.frame(head(s$Waccamaw_Alligator, n = 10)) %>%
  mutate(Comparison = "Waccamaw_Alligator",
         ASV = rownames(.))
s5 <- as.data.frame(head(s$Waccamaw_SF, n = 10)) %>%
  mutate(Comparison = "Waccamaw_SF",
         ASV = rownames(.))
s6 <- as.data.frame(head(s$Alligator_SF, n = 10)) %>%
  mutate(Comparison = "Alligator_SF",
         ASV = rownames(.))
simper_results <- rbind(s1, s2, s3, s4, s5, s6) %>%
  left_join(., input_filt_rare_abund$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
saveRDS(simper_results, "simper_results.RDS")
simper_results <- readRDS("simper_results.RDS")
length(unique(simper_results$ASV)) # 27 unique ASVs
write_xlsx(simper_results, 
           "simper_results_Comb.xlsx",
           format_headers = F)

# MULTIPATT (list ASVs associated with each group)
# Don't plot here because later do at Field Control/Lab/Field Exp. level
set.seed(1202)
mp <- multipatt(t(input_filt_rare_abund$data_loaded), 
                input_filt_rare_abund$map_loaded$Estuary, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp)
View(mp$sign)
multipatt_results <- mp$sign %>%
  filter(p.value == 0.001) %>%
  mutate(num_sites = rowSums(cbind(s.Waccamaw, s.Alligator, s.Delaware, s.SF))) %>%
  filter(num_sites == 1) %>%
  mutate(Group = "NA")
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$s.Waccamaw[i] == 1) {
    multipatt_results$Group[i] <- "Waccamaw"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$s.Alligator[i] == 1) {
    multipatt_results$Group[i] <- "Alligator"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$s.Delaware[i] == 1) {
    multipatt_results$Group[i] <- "Delaware"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$s.SF[i] == 1) {
    multipatt_results$Group[i] <- "SF"
  }
}
table(multipatt_results$Group)



#### _Domain ####
tax_sum_domain <- summarize_taxonomy(input_filt_rare, level = 1, report_higher_tax = F)
plot_ts_heatmap(tax_sum_domain, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_domain,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Domain") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(tax_sum_domain, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')



#### _Phylum ####
tax_sum_phyla <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = F)
plot_ts_heatmap(tax_sum_phyla, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_phyla,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 10,
                       data_only = TRUE)
pdf("Figs/CombPhyla.pdf", width = 7, height = 5)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(tax_sum_phyla, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Look at archaeal phyla
tax_sum_phyla_ar <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = T)
tax_sum_phyla_ar <- tax_sum_phyla_ar[grep("Archaea", rownames(tax_sum_phyla_ar)),]
bars_ar <- plot_taxa_bars(tax_sum_phyla_ar,
                          input_filt_rare$map_loaded,
                          "Estuary",
                          num_taxa = 13,
                          data_only = TRUE)
nb.cols <- nrow(tax_sum_phyla_ar)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_ar, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

#### _Class ####
# Check sulfur reducers
tax_sum_phyla_su <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = T)
tax_sum_phyla_su <- tax_sum_phyla_su[grep("Desulfo", rownames(tax_sum_phyla_su)),]
bars_su <- plot_taxa_bars(tax_sum_phyla_su,
                          input_filt_rare$map_loaded,
                          "Estuary",
                          num_taxa = 13,
                          data_only = TRUE) %>%
  mutate(taxon = substring(taxon, 11))
nb.cols <- nrow(tax_sum_phyla_su)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_su, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))



#### _Family ####
tax_sum_families <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = FALSE)
plot_ts_heatmap(tax_sum_families, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars2 <- plot_taxa_bars(tax_sum_families,
                        input_filt_rare$map_loaded,
                        "Estuary",
                        num_taxa = 10,
                        data_only = TRUE)
ggplot(bars2, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Genus") +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
taxa_summary_by_sample_type(tax_sum_families, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Look at methanogens
tax_sum_families_meth <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = F)
tax_sum_families_meth <- tax_sum_families_meth[grep("Methano", rownames(tax_sum_families_meth)),]
bars_meth <- plot_taxa_bars(tax_sum_families_meth,
                            input_filt_rare$map_loaded,
                            "Estuary",
                            num_taxa = 14,
                            data_only = TRUE)
nb.cols <- nrow(tax_sum_families_meth)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_meth, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))



#### _Genus ####
tax_sum_genera <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = TRUE)
plot_ts_heatmap(tax_sum_genera, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_genera,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Genus") +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
taxa_summary_by_sample_type(tax_sum_genera, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Guilds ####
tax_sum_guilds <- summarize_taxonomy(input_filt_rare, level = 9, report_higher_tax = F)
plot_ts_heatmap(tax_sum_guilds, 
                input_filt_rare$map_loaded, 
                0, 
                'Estuary',
                rev_taxa = T,
                remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_guilds,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 20,
                       data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
pdf("Figs/CombGuilds.pdf", width = 7, height = 5)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(tax_sum_guilds, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')



#### _Venn ####
phy <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = F)
cla <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = F)
ord <- summarize_taxonomy(input_filt_rare, level = 4, report_higher_tax = F)
fam <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = F)
gen <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = F)

input_phylum <- input_filt_rare
input_phylum$data_loaded <- phy
input_class <- input_filt_rare
input_class$data_loaded <- cla
input_order <- input_filt_rare
input_order$data_loaded <- ord
input_family <- input_filt_rare
input_family$data_loaded <- fam
input_genus <- input_filt_rare
input_genus$data_loaded <- gen

pdf("Figs/CombVennOTU.pdf", width = 7, height = 5)
plot_venn_diagram(input_filt_rare,
                  "Estuary",
                  0.00000000000000001)
dev.off()

pdf("Figs/CombVennAll.pdf", width = 15, height = 6.5)
plot_grid(plot_venn_diagram(input_phylum, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_class, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_order, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_family, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_genus, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_filt_rare, "Estuary", 0.00000000000000001),
          labels = c("(a) Phylum", "(b) Class", "(c) Order", 
                     "(d) Family", "(e) Genus", "(f) OTU"))
dev.off()



#### ...................................... ####
#### Comparison Field Control ####
# Only look at field soils (no manipulations or incubations)
# Classify as fresh/oligo/meso
# All sites
# Archaea
# Bacteria

#### __Setup ####
input_filt_rare <- readRDS("input_filt_rare_comb.rds")
input_filt_rare$map_loaded <- input_filt_rare$map_loaded %>%
  mutate(Estuary = factor(Estuary,
                          levels = c("Waccamaw", "Alligator", "Delaware", "SF")),
         sampleID = rownames(.),
         Field = "NA",
         Salt = "NA",
         Depth = gsub("0-5", 0.05, Depth)) %>%
  mutate(Depth = gsub("5-15", 0.15, Depth)) %>%
  mutate(Depth = as.numeric(Depth)) %>%
  mutate(Depth = ifelse(Depth <= 0.05, "< 5 cm", "5 - 15 cm"))

# Add column for if unmanipulated field sample
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "SF") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Delaware" & 
      input_filt_rare$map_loaded$Site[i] == "Soil") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Alligator" & 
      input_filt_rare$map_loaded$Site[i] == "Soil") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Waccamaw" & 
      input_filt_rare$map_loaded$Detail[i] == "Control") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
field <- filter_data(input_filt_rare,
                     filter_cat = "Field",
                     keep_vals = "Field")

# Add column for Salt
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Site[i] == "Sandmound" | field$map_loaded$Site[i] == "West Pond") {
    field$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Site[i] == "Mayberry" | field$map_loaded$Site[i] == "Browns") {
    field$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Site[i] == "Joice" | field$map_loaded$Site[i] == "Rush Ranch") {
    field$map_loaded$Salt[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Site[i] == "Goodyear" | field$map_loaded$Site[i] == "White Slough" |
      field$map_loaded$Site[i] == "Tolay" | field$map_loaded$Site[i] == "China Camp" |
      field$map_loaded$Site[i] == "Muzzi") {
    field$map_loaded$Salt[i] <- "Polyhaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Estuary[i] == "Waccamaw") {
    field$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Estuary[i] == "Alligator") {
    field$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Info[i] == "OligoHal_source") {
    field$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Info[i] == "MesoHal_source") {
    field$map_loaded$Salt[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Info[i] == "TFM1_source" | field$map_loaded$Info[i] == "TFM2_source") {
    field$map_loaded$Salt[i] <- "Freshwater"
  }
}
field$map_loaded <- field$map_loaded %>%
  unite("EstSalt", c(Estuary, Salt), sep = "_", remove = F) %>%
  mutate(EstSalt = factor(EstSalt,
                          levels = c("SF_Freshwater", "Alligator_Freshwater",
                                     "Delaware_Freshwater", "Waccamaw_Freshwater",
                                     "SF_Oligohaline", "Delaware_Oligohaline",
                                     "SF_Mesohaline", "Delaware_Mesohaline",
                                     "SF_Polyhaline")),
         Salt = factor(Salt,
                       levels = c("Freshwater", "Oligohaline", 
                                  "Mesohaline", "Polyhaline")))

#### __Alpha ####
leveneTest(field$map_loaded$rich ~ field$map_loaded$Salt) # Almost homogeneous
m <- aov(rich ~ Estuary + Salt + Depth, data = field$map_loaded)
Anova(m, type = "III") # Estuary, Salt sig. Depth marginal.
m <- aov(rich ~ Salt, data = field$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(field$map_loaded$rich)+(max(field$map_loaded$rich)-min(field$map_loaded$rich))/20)
leveneTest(field$map_loaded$shannon ~ field$map_loaded$Salt) # Homogeneous
m1 <- aov(shannon ~ Estuary + Salt + Depth, data = field$map_loaded)
Anova(m1, type = "III") # All sig
m1 <- aov(shannon ~ Salt, data = field$map_loaded)
shapiro.test(m1$residuals) # Not normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(field$map_loaded$shannon)+(max(field$map_loaded$shannon)-min(field$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- field$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/Comb_Control_Alpha.pdf", width = 8, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth, colour = Estuary)) +
  geom_text(data = label_df, aes(Salt, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs") +
  scale_colour_viridis_d() +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
field_bc <- calc_dm(field$data_loaded)
set.seed(1150)
adonis2(field_bc ~ field$map_loaded$Estuary+field$map_loaded$Salt+field$map_loaded$Depth) # All
anova(betadisper(field_bc, field$map_loaded$Estuary)) # Dispersion not homogeneous
anova(betadisper(field_bc, field$map_loaded$Salt)) # Dispersion not homogeneous
anova(betadisper(field_bc, field$map_loaded$Depth)) # Dispersion homogeneous
field_pcoa <- cmdscale(field_bc, k = nrow(field$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(field_pcoa)/sum(eigenvals(field_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(field_pcoa)/sum(eigenvals(field_pcoa)))[2]*100, digits = 1)
field$map_loaded$Axis01 <- scores(field_pcoa)[,1]
field$map_loaded$Axis02 <- scores(field_pcoa)[,2]
micro.hulls <- ddply(field$map_loaded, c("Salt"), find_hull)
pdf("Figs/Comb_Control_PCoA.pdf", width = 7, height = 5)
ggplot(field$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Estuary)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salinity") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d(guide = guide_legend(override.aes = list(shape = 15))) +
  scale_shape_manual(values = c(16, 17, 18, 3)) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

#### __BC Among Within ####
# Plot BC dissimilarity across site and salt combos
# Is Fresh vs. Poly at same site more different than Fresh vs Fresh at different sites?
# Compare SF_Fresh to the other three fresh and to SF_Oligo, SF_Meso, SF_Poly

# Convert from dist object to matrix object
bac_bray_mat <- as.matrix(field_bc)
# Remove duplicates and diagonal by setting upper triangle and diagonal to NA
bac_bray_mat[upper.tri(bac_bray_mat, diag = TRUE)] <- NA
# Convert from matrix object to dataframe
bac_bray_df <- as.data.frame(bac_bray_mat)
# Make a sample ID column from the row names.The other sample in each pairwise comparison will be the variable column after melting
bac_bray_df$sampleID <- rownames(bac_bray_df)
# Melt with reshape2 package
bac_bray_df_long <- melt(bac_bray_df, id.vars = "sampleID")
# Get rid of the NA's (duplicates and diagonal)
bac_bray_df_long <- na.omit(bac_bray_df_long)
# Note, the length of this dataframe should now equal (n*(n-1))/2
nrow(bac_bray_df_long) == (nrow(field$map_loaded)*(nrow(field$map_loaded)-1))/2 # Good!
# Make sampleID a factor
bac_bray_df_long$sampleID <- as.factor(bac_bray_df_long$sampleID)
# Now add EstSalt, matching to col1 and col2
# This will give the estuary and salinity for each of the two samples in each pairwise comparison
EstSalt <- dplyr::select(field$map_loaded, sampleID, EstSalt)
bac_bray_df_long <- inner_join(bac_bray_df_long, EstSalt, 
                               by = c("sampleID" = "sampleID"))
bac_bray_df_long <- inner_join(bac_bray_df_long, EstSalt, 
                               by = c("variable" = "sampleID"))
# Make new column for comparison and filter to comparisons of interest
bac_bray_df_long$comparison <- "NA"
for (i in 1:nrow(bac_bray_df_long)) {
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Delaware_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Alligator_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Waccamaw_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "Delaware_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Alligator_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "Delaware_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Waccamaw_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "Alligator_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Waccamaw_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "Delaware_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Waccamaw_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "SF_Oligohaline") {
    bac_bray_df_long$comparison[i] <- "SF, Freshwater/Oligohaline"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "SF_Mesohaline") {
    bac_bray_df_long$comparison[i] <- "SF, Freshwater/Mesohaline"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "SF_Polyhaline") {
    bac_bray_df_long$comparison[i] <- "SF, Freshwater/Polyhaline"
  }
}
# Make the new column a factor
bac_bray_df_long$comparison <- as.factor(bac_bray_df_long$comparison)
# Filter unwanted comparisons
bac_bray_df_long <- subset(bac_bray_df_long, comparison != "NA")
# Check the sample sizes
table(bac_bray_df_long$comparison)
# ANOVA
m <- aov(value ~ comparison, data = bac_bray_df_long)
summary(m)
t <- emmeans(object = m, specs = "comparison") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "value",
         y = max(bac_bray_df_long$value)+(max(bac_bray_df_long$value)-min(bac_bray_df_long$value))/20)
bc_plot <- ggplot(data = bac_bray_df_long, aes(comparison, value, colour = comparison)) +
  geom_jitter(size = 0.75, alpha = 0.2, width = 0.3) +
  geom_boxplot(outlier.shape = NA, colour = "black", fill = NA) +
  geom_text(data = t, aes(comparison, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site Comparison",
       y = "Bray-Curtis Dissimilarity") +
  scale_color_brewer(palette = "Set2") +
  ylim(0.55, 1.05) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
bc_plot

# Interactive
ggplotly(bc_plot)

# With density plot
yplot <- ggdensity(bac_bray_df_long, "value", fill = "comparison", 
                   palette = "Set2", size = 0.25) +
  rotate() + 
  clean_theme() + 
  rremove("legend") + 
  xlim(0.55, 1.05)
yplot
plot_final <- insert_yaxis_grob(bc_plot, yplot, position = "right")
pdf("Figs/Comb_Control_BC.pdf", width = 6, height = 5)
ggdraw(plot_final)
dev.off()

#### __Taxa ####
field_phyla <- summarize_taxonomy(field, level = 2, report_higher_tax = F)
plot_ts_heatmap(field_phyla, field$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsP <- plot_taxa_bars(field_phyla, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/Comb_Control_Phyla.pdf", width = 7, height = 5)
ggplot(field_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_phyla, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_class <- summarize_taxonomy(field, level = 3, report_higher_tax = F)
plot_ts_heatmap(field_class, field$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsC <- plot_taxa_bars(field_class, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/Comb_Control_Class.pdf", width = 7, height = 5)
ggplot(field_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_class, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_order <- summarize_taxonomy(field, level = 4, report_higher_tax = F)
plot_ts_heatmap(field_order, field$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsO <- plot_taxa_bars(field_order, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/Comb_Control_Order.pdf", width = 7, height = 5)
ggplot(field_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_order, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_family <- summarize_taxonomy(field, level = 5, report_higher_tax = F)
plot_ts_heatmap(field_family, field$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsF <- plot_taxa_bars(field_family, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/Comb_Control_Family.pdf", width = 7, height = 5)
ggplot(field_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_family, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_genus <- summarize_taxonomy(field, level = 6, report_higher_tax = F)
plot_ts_heatmap(field_genus, field$map_loaded, 0.005, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsG <- plot_taxa_bars(field_genus, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("Figs/Comb_Control_Genus.pdf", width = 7, height = 5)
ggplot(field_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_genus, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_guilds <- summarize_taxonomy(field, level = 9, report_higher_tax = F)
plot_ts_heatmap(field_guilds, field$map_loaded, 0, 'EstSalt', rev_taxa = T, remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsGu <- plot_taxa_bars(field_guilds, field$map_loaded, "EstSalt",
                               num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
pdf("Figs/Comb_Control_Guilds.pdf", width = 7, height = 6)
ggplot(field_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_guilds, field$map_loaded, 'EstSalt', 0.01, 'KW')

#### __Methano ####
# Summarize by family, extract methanogens
tax_sum_family_wTax <- summarize_taxonomy(field, level = 5, report_higher_tax = T, 
                                          relative = FALSE)
methano_wTax <- tax_sum_family_wTax[grep("Methano", rownames(tax_sum_family_wTax)),]
# 20, including 5 NA

# Without higher tax, remove NAs, relative abundance
tax_sum_family <- summarize_taxonomy(field, level = 5, report_higher_tax = FALSE, 
                                     relative = TRUE)
methano <- tax_sum_family[grep("Methano", rownames(tax_sum_family)),]
field_barsMethano <- plot_taxa_bars(methano, field$map_loaded, "EstSalt", 
                              num_taxa = nrow(methano), data_only = T)
nb.cols <- nrow(methano)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
pdf("Figs/Comb_Control_Methano.pdf", width = 7, height = 5)
ggplot(field_barsMethano, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Family") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

#### __Desulfo ####
# Summarize by family, extract sulfate reducers
tax_sum_family_wTax <- summarize_taxonomy(field, level = 5, report_higher_tax = T, 
                                          relative = FALSE)
desulfo_wTax <- tax_sum_family_wTax[grep("Desulfo", rownames(tax_sum_family_wTax)),]
# 49, including some NA. Note, some Desulfobacterota and some Firmicutes
# Too many, try order
tax_sum_order_wTax <- summarize_taxonomy(field, level = 4, report_higher_tax = T, 
                                          relative = TRUE)
desulfo_wTax <- tax_sum_order_wTax[grep("Desulfo", rownames(tax_sum_order_wTax)),]
# 21, including some NA
# Trim off bacteria
rownames(desulfo_wTax) <- substring(rownames(desulfo_wTax), 11)
nb.cols <- nrow(desulfo_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

field_barsDesulfo <- plot_taxa_bars(desulfo_wTax, field$map_loaded, "EstSalt", 
                                    num_taxa = nrow(desulfo_wTax), data_only = T)
pdf("Figs/Comb_Control_Desulfo.pdf", width = 7, height = 7)
ggplot(field_barsDesulfo, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = mycolors) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

#### __Simper ####
field_sim <- simper(t(field$data_loaded), 
                    field$map_loaded$Salt)
field_s <- summary(field_sim)
head(field_s$Freshwater_Oligohaline)
head(field_s$Freshwater_Mesohaline)
head(field_s$Freshwater_Polyhaline)
head(field_s$Oligohaline_Mesohaline)
head(field_s$Oligohaline_Polyhaline)
head(field_s$Mesohaline_Polyhaline)
field_df1 <- head(field_s$Freshwater_Oligohaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Oligohaline")
field_df2 <- head(field_s$Freshwater_Mesohaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline")
field_df3 <- head(field_s$Freshwater_Polyhaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Polyhaline")
field_simper_results <- rbind(field_df1, field_df2, field_df3) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline") %>%
  left_join(., field$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Increase", "Decrease")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanSalt" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanSalt, MeanControl, CumulativeContribution)
write_xlsx(field_simper_results, 
           "simper_results_comb_control.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
field_mp <- multipatt(t(field$data_loaded), 
                      field$map_loaded$Salt, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 1)
field_mp_results <- field_mp$sign %>%
  mutate(q.value = qvalue(field_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(field_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.65) %>%
  left_join(., field$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Freshwater[i] == 1) {
    field_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Oligohaline[i] == 1) {
    field_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Mesohaline[i] == 1) {
    field_mp_results$Group[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Polyhaline[i] == 1) {
    field_mp_results$Group[i] <- "Polyhaline"
  }
}
table(field_mp_results$Group)
field_asv <- summarize_taxonomy(field, level = 8, report_higher_tax = F)
field_asv_all <- data.frame("RelAbundance" = round(rowMeans(field_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
field_mp_corrs <- as.data.frame(field_mp$str) %>%
  dplyr::select(1:length(levels(field$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% field_mp_results$ASV) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesolhaline", "Polyhaline", "ASV"))
# Add corrs and taxonomy
field_mp_results <- field_mp_results %>%
  left_join(., field_asv_all, by = "ASV") %>%
  left_join(., field_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(field_mp_corrs)[1:4], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(field_mp_corrs)[1:4], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
field.hm.melted <- field_mp_results %>%
  dplyr::select(taxon, names(field_mp_corrs)[1:4]) %>%
  melt(., id.vars = c("taxon"))
field.hm <- ggplot(data = field.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(field.hm.melted$taxon), labels = unique(field.hm.melted$taxon),
                   limits = rev(levels(field.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
field.l <- get_legend(field.hm)
field.hm.clean <- field.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
field.bp.y <- ggplot(data = field_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(field_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/Comb_Control_Multipatt.pdf", width = 8, height = 6)
plot_grid(field.hm.clean, field.bp.y, field.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, many OTUs are too rare, not even 1%. So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
field_tax <- field$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
field_genus <- summarize_taxonomy(field, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
field_mp <- multipatt(t(field_genus), 
                      field$map_loaded$Salt, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 1)
field_mp_results <- field_mp$sign %>%
  mutate(q.value = qvalue(field_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(field_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.57) %>%
  left_join(., field_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Freshwater[i] == 1) {
    field_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Oligohaline[i] == 1) {
    field_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Mesohaline[i] == 1) {
    field_mp_results$Group[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Polyhaline[i] == 1) {
    field_mp_results$Group[i] <- "Polyhaline"
  }
}
table(field_mp_results$Group)
field_genus_all <- data.frame("RelAbundance" = round(rowMeans(field_genus)/min(colSums(field$data_loaded)) * 100, digits = 4)) %>%
  mutate(Genus = rownames(.))
field_mp_corrs <- as.data.frame(field_mp$str) %>%
  dplyr::select(1:length(levels(field$map_loaded$Salt))) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% field_mp_results$Genus) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesohaline", "Polyhaline", "Genus"))
# Add corrs and taxonomy
field_mp_results <- field_mp_results %>%
  left_join(., field_genus_all, by = "Genus") %>%
  left_join(., field_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(field_mp_corrs)[1:4], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(field_mp_corrs)[1:4], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
field.hm.melted <- field_mp_results %>%
  dplyr::select(taxon, names(field_mp_corrs)[1:4]) %>%
  melt(., id.vars = c("taxon"))
field.hm <- ggplot(data = field.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(field.hm.melted$taxon), labels = unique(field.hm.melted$taxon),
                   limits = rev(levels(field.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
field.l <- get_legend(field.hm)
field.hm.clean <- field.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
field.bp.y <- ggplot(data = field_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(field_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/Comb_Control_Multipatt_genus.pdf", width = 8, height = 6)
plot_grid(field.hm.clean, field.bp.y, field.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### Comparison Lab Inc ####
# Only look at incubations
# Only look at control and + ASW
# Only look at final time point
# Delaware and North Carolina
# Pull out methanogens and sulfate reducers
#### __Setup ####
input_filt_rare <- readRDS("input_filt_rare_comb.rds")
input_filt_rare$map_loaded <- input_filt_rare$map_loaded %>%
  mutate(Estuary = factor(Estuary,
                          levels = c("Waccamaw", "Alligator", "Delaware", "SF")),
         sampleID = rownames(.),
         Field = "NA",
         Salt = "NA",
         Depth = gsub("0-5", 0.05, Depth)) %>%
  mutate(Depth = gsub("5-15", 0.15, Depth)) %>%
  mutate(Depth = as.numeric(Depth)) %>%
  mutate(Depth = ifelse(Depth <= 0.05, "< 5 cm", "5 - 15 cm"))

# Add column for if unmanipulated field sample
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "SF") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Delaware" & 
      input_filt_rare$map_loaded$Site[i] == "Soil") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Alligator" & 
      input_filt_rare$map_loaded$Site[i] == "Soil") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Waccamaw" & 
      input_filt_rare$map_loaded$Detail[i] == "Control") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
lab <- filter_data(input_filt_rare,
                   filter_cat = "Field",
                   keep_vals = "NA") # 141 samples
table(lab$map_loaded$Estuary)

# Now need to get only final time point and only incubations
# Note that there are 2 incubation experiments
# Both used 5 ppt ASW and went 12 weeks
lab <- filter_data(lab,
                   filter_cat = "Site",
                   keep_vals = c("Soil incubation", "Soil Incubation")) # 67 samples
table(lab$map_loaded$Estuary)

# Note all "Alligator" are 3 months, so get "12 wk" for Delaware
table(lab$map_loaded$Detail)

lab <- filter_data(lab,
                   filter_cat = "Detail",
                   keep_vals = c("5 ppt ASW 12wk",
                                 "Freshwater 12wk",
                                 "5ppt ASW",
                                 "DI_ctrl")) # 26 samples
table(lab$map_loaded$Estuary)

# Make Treatment Column
lab$map_loaded <- lab$map_loaded %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Salt = recode_factor(Detail,
                              "Freshwater 12wk" = "Control",
                              "DI_ctrl" = "Control",
                              "5 ppt ASW 12wk" = "+ASW",
                              "5ppt ASW" = "+ASW"))

# Need combined estuary/salt factor
lab$map_loaded$EstSalt <- paste(lab$map_loaded$Estuary,
                                lab$map_loaded$Salt,
                                sep = "_")

#### __Alpha ####
leveneTest(lab$map_loaded$rich ~ lab$map_loaded$Salt) # Homogeneous
m <- aov(rich ~ Estuary + Salt + Depth, data = lab$map_loaded)
Anova(m, type = "III") # Estuary sig, salt and depth not sig.
shapiro.test(m$residuals) # Almost normal

leveneTest(lab$map_loaded$shannon ~ lab$map_loaded$Salt) # Homogenoeus
m1 <- aov(shannon ~ Estuary + Salt + Depth, data = lab$map_loaded)
Anova(m1, type = "III") # Estuary sig, salt and depth not sig.
shapiro.test(m1$residuals) # Not normal
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- lab$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/Comb_Lab_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, aes(shape = Depth, colour = Estuary)) +
  labs(x = "Site", y = "Number of OTUs") +
  scale_colour_viridis_d() +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()



#### __Beta ####
lab_bc <- calc_dm(lab$data_loaded)
set.seed(1150)
adonis2(lab_bc ~ lab$map_loaded$Estuary+lab$map_loaded$Salt+lab$map_loaded$Depth) # All
adonis2(lab_bc ~ lab$map_loaded$Depth+lab$map_loaded$Salt+lab$map_loaded$Estuary)
anova(betadisper(lab_bc, lab$map_loaded$Estuary)) # Dispersion homogeneous
anova(betadisper(lab_bc, lab$map_loaded$Salt)) # Dispersion homogeneous
anova(betadisper(lab_bc, lab$map_loaded$Depth)) # Dispersion homogeneous
lab_pcoa <- cmdscale(lab_bc, k = nrow(lab$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(lab_pcoa)/sum(eigenvals(lab_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(lab_pcoa)/sum(eigenvals(lab_pcoa)))[2]*100, digits = 1)
lab$map_loaded$Axis01 <- scores(lab_pcoa)[,1]
lab$map_loaded$Axis02 <- scores(lab_pcoa)[,2]
micro.hulls <- ddply(lab$map_loaded, c("EstSalt"), find_hull)
pdf("Figs/Comb_Lab_PCoA.pdf", width = 7, height = 5)
ggplot(lab$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = EstSalt, fill = EstSalt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = EstSalt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salt") +
  scale_fill_manual(values = c("red", "blue", "red", "blue")) +
  scale_colour_manual(values = c("red", "blue", "red", "blue")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   shape = 15)),
         fill = "none") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



#### __Taxa ####
lab_phyla <- summarize_taxonomy(lab, level = 2, report_higher_tax = F)
plot_ts_heatmap(lab_phyla, lab$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsP <- plot_taxa_bars(lab_phyla, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Lab_Phyla.pdf", width = 7, height = 5)
ggplot(lab_barsP, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(lab_phyla, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_class <- summarize_taxonomy(lab, level = 3, report_higher_tax = F)
plot_ts_heatmap(lab_class, lab$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsC <- plot_taxa_bars(lab_class, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Lab_Class.pdf", width = 7, height = 5)
ggplot(lab_barsC, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(lab_class, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_order <- summarize_taxonomy(lab, level = 4, report_higher_tax = F)
plot_ts_heatmap(lab_order, lab$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsO <- plot_taxa_bars(lab_order, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Lab_Order.pdf", width = 7, height = 5)
ggplot(lab_barsO, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(lab_order, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_family <- summarize_taxonomy(lab, level = 5, report_higher_tax = F)
plot_ts_heatmap(lab_family, lab$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsF <- plot_taxa_bars(lab_family, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Lab_Family.pdf", width = 7, height = 5)
ggplot(lab_barsF, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(lab_family, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_genus <- summarize_taxonomy(lab, level = 6, report_higher_tax = F)
plot_ts_heatmap(lab_genus, lab$map_loaded, 0.005, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsG <- plot_taxa_bars(lab_genus, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Lab_Genus.pdf", width = 7, height = 5)
ggplot(lab_barsG, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(lab_genus, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_guilds <- summarize_taxonomy(lab, level = 9, report_higher_tax = F)
plot_ts_heatmap(lab_guilds, lab$map_loaded, 0, 'EstSalt', rev_taxa = T, remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsGu <- plot_taxa_bars(lab_guilds, lab$map_loaded, "EstSalt",
                               num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Lab_Guilds.pdf", width = 7, height = 5)
ggplot(lab_barsGu, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(lab_guilds, lab$map_loaded, 'EstSalt', 0.01, 'KW')

#### __Methano ####
# Summarize by family, extract methanogens
tax_sum_family_wTax <- summarize_taxonomy(lab, level = 5, report_higher_tax = T, 
                                          relative = FALSE)
methano_wTax <- tax_sum_family_wTax[grep("Methano", rownames(tax_sum_family_wTax)),]
# 15, including 4 NA

# Without higher tax, remove NAs, relative abundance
tax_sum_family <- summarize_taxonomy(lab, level = 5, report_higher_tax = FALSE, 
                                     relative = TRUE)
methano <- tax_sum_family[grep("Methano", rownames(tax_sum_family)),]
# 10
lab_barsMethano <- plot_taxa_bars(methano, lab$map_loaded, "EstSalt", 
                                    num_taxa = nrow(methano), data_only = T) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Lab_Methano.pdf", width = 7, height = 5)
ggplot(lab_barsMethano, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Family") +
  scale_fill_brewer(palette = "Paired") +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

#### __Desulfo ####
# Summarize by family, extract sulfate reducers
tax_sum_family_wTax <- summarize_taxonomy(lab, level = 5, report_higher_tax = T, 
                                          relative = FALSE)
desulfo_wTax <- tax_sum_family_wTax[grep("Desulfo", rownames(tax_sum_family_wTax)),]
# 34, including some NA. Note, some Desulfobacterota and some Firmicutes
# Too many, try order
tax_sum_order_wTax <- summarize_taxonomy(lab, level = 4, report_higher_tax = T, 
                                         relative = TRUE)
desulfo_wTax <- tax_sum_order_wTax[grep("Desulfo", rownames(tax_sum_order_wTax)),]
# 17, including 2 NA
# Trim off bacteria
rownames(desulfo_wTax) <- substring(rownames(desulfo_wTax), 11)
nb.cols <- nrow(desulfo_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

lab_barsDesulfo <- plot_taxa_bars(desulfo_wTax, lab$map_loaded, "EstSalt", 
                                    num_taxa = nrow(desulfo_wTax), data_only = T) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Lab_Desulfo.pdf", width = 7, height = 5)
ggplot(lab_barsDesulfo, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

#### __Simper ####
lab_sim <- simper(t(lab$data_loaded), 
                    lab$map_loaded$Salt)
lab_s <- summary(lab_sim)
head(lab_s$`Control_+ASW`)
lab_df1 <- head(lab_s$`Control_+ASW`, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Control/+ASW") %>%
  left_join(., lab$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Increase", "Decrease")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanSalt" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanSalt, MeanControl, CumulativeContribution)
write_xlsx(lab_df1, 
           "simper_results_comb_lab.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
lab_mp <- multipatt(t(lab$data_loaded), 
                      lab$map_loaded$Salt, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 1)
lab_mp_results <- lab_mp$sign %>%
  mutate(q.value = qvalue(lab_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(lab_mp$sign)) %>%
  filter(p.value < 0.01) %>%
  filter(stat >= 0.5) %>%
  left_join(., lab$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(lab_mp_results)) {
  if (lab_mp_results$s.Control[i] == 1) {
    lab_mp_results$Group[i] <- "Control"
  }
}
for (i in 1:nrow(lab_mp_results)) {
  if (lab_mp_results$`s.+ASW`[i] == 1) {
    lab_mp_results$Group[i] <- "+ASW"
  }
}
table(lab_mp_results$Group)
lab_asv <- summarize_taxonomy(lab, level = 8, report_higher_tax = F)
lab_asv_all <- data.frame("RelAbundance" = round(rowMeans(lab_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
lab_mp_corrs <- as.data.frame(lab_mp$str) %>%
  dplyr::select(1:length(levels(lab$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% lab_mp_results$ASV) %>%
  set_names(c("Control", "+ASW", "ASV"))
# Add corrs and taxonomy
lab_mp_results <- lab_mp_results %>%
  left_join(., lab_asv_all, by = "ASV") %>%
  left_join(., lab_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(lab_mp_corrs)[1:2], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(lab_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
lab.hm.melted <- lab_mp_results %>%
  dplyr::select(taxon, names(lab_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
lab.hm <- ggplot(data = lab.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(lab.hm.melted$taxon), labels = unique(lab.hm.melted$taxon),
                   limits = rev(levels(lab.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
lab.l <- get_legend(lab.hm)
lab.hm.clean <- lab.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)), size = 5),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
lab.bp.y <- ggplot(data = lab_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(lab_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/Comb_Lab_Multipatt.pdf", width = 8, height = 10)
plot_grid(lab.hm.clean, lab.bp.y, lab.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### Comparison Field Exp ####
# Only look at field manipulations
# South Carolina and Delaware
# Delaware just look at transplant to oligo which is similar to seawater addition
# Just look at surface in Delaware, similar to SC depth
# "Soil Field plots", and "Soil mesocosm"
#### __Setup ####
input_filt_rare <- readRDS("input_filt_rare_comb.rds")
input_filt_rare$map_loaded <- input_filt_rare$map_loaded %>%
  mutate(Estuary = factor(Estuary,
                          levels = c("Waccamaw", "Alligator", "Delaware", "SF")),
         sampleID = rownames(.),
         Field = "NA",
         Salt = "NA",
         Depth = gsub("0-5", 0.05, Depth)) %>%
  mutate(Depth = gsub("5-15", 0.15, Depth)) %>%
  mutate(Depth = as.numeric(Depth)) %>%
  mutate(Depth = ifelse(Depth <= 0.05, "< 5 cm", "5 - 15 cm"))

exp <- filter_data(input_filt_rare,
                   filter_cat = "Site",
                   keep_vals = c("Soil Field plots", "Soil mesocosm")) # 84 samples
table(exp$map_loaded$Estuary)

# Get surface transplants to TFM and OligoHal from Delaware
# Get freshwater amended and saltwater amended from SC
exp <- filter_data(exp,
                   filter_cat = "Detail",
                   keep_vals = c("TFM1@TFM2",
                                 "TFM1@OligoHal",
                                 "Freshwater amended",
                                 "Saltwater amended")) # 28 samples
table(exp$map_loaded$Estuary) # 20 SC, 8 Del

# Make Treatment Column
exp$map_loaded <- exp$map_loaded %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Salt = recode_factor(Detail,
                              "TFM1@TFM2" = "Control",
                              "Freshwater amended" = "Control",
                              "TFM1@OligoHal" = "+ASW",
                              "Saltwater amended" = "+ASW"))

# Need combined estuary/salt factor
exp$map_loaded$EstSalt <- paste(exp$map_loaded$Estuary,
                                exp$map_loaded$Salt,
                                sep = "_")

#### __Alpha ####
leveneTest(exp$map_loaded$rich ~ exp$map_loaded$Salt) # Homogeneous
m <- aov(rich ~ Estuary + Salt + Depth, data = exp$map_loaded)
Anova(m, type = "III") # Salt sig, estuary and depth not sig
shapiro.test(m$residuals) # Normal

leveneTest(exp$map_loaded$shannon ~ exp$map_loaded$Salt) # Homogenoeus
m1 <- aov(shannon ~ Estuary + Salt + Depth, data = exp$map_loaded)
Anova(m1, type = "III") # Salt sig, estuary and depth not sig
shapiro.test(m1$residuals) # Not normal
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- exp$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/Comb_Exp_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, aes(shape = Depth, colour = Estuary)) +
  labs(x = "Site", y = "Number of OTUs") +
  scale_colour_viridis_d() +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
exp_bc <- calc_dm(exp$data_loaded)
set.seed(1150)
adonis2(exp_bc ~ Estuary + Salt + Depth, data = exp$map_loaded) # All
adonis2(exp_bc ~ Depth + Salt + Estuary, data = exp$map_loaded) # All
anova(betadisper(exp_bc, exp$map_loaded$Estuary)) # Dispersion homogeneous
anova(betadisper(exp_bc, exp$map_loaded$Salt)) # Dispersion homogeneous
anova(betadisper(exp_bc, exp$map_loaded$Depth)) # Dispersion homogeneous
exp_pcoa <- cmdscale(exp_bc, k = nrow(exp$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(exp_pcoa)/sum(eigenvals(exp_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(exp_pcoa)/sum(eigenvals(exp_pcoa)))[2]*100, digits = 1)
exp$map_loaded$Axis01 <- scores(exp_pcoa)[,1]
exp$map_loaded$Axis02 <- scores(exp_pcoa)[,2]
micro.hulls <- ddply(exp$map_loaded, c("EstSalt"), find_hull)
pdf("Figs/Comb_Exp_PCoA.pdf", width = 7, height = 5)
ggplot(exp$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = EstSalt, fill = EstSalt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = EstSalt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salt") +
  scale_fill_manual(values = c("red", "blue", "red", "blue")) +
  scale_colour_manual(values = c("red", "blue", "red", "blue")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   shape = 15)),
         fill = "none") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

#### __Taxa ####
exp_phyla <- summarize_taxonomy(exp, level = 2, report_higher_tax = F)
plot_ts_heatmap(exp_phyla, exp$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsP <- plot_taxa_bars(exp_phyla, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Exp_Phyla.pdf", width = 7, height = 5)
ggplot(exp_barsP, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(exp_phyla, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_class <- summarize_taxonomy(exp, level = 3, report_higher_tax = F)
plot_ts_heatmap(exp_class, exp$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsC <- plot_taxa_bars(exp_class, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Exp_Class.pdf", width = 7, height = 5)
ggplot(exp_barsC, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(exp_class, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_order <- summarize_taxonomy(exp, level = 4, report_higher_tax = F)
plot_ts_heatmap(exp_order, exp$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsO <- plot_taxa_bars(exp_order, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Exp_Order.pdf", width = 7, height = 5)
ggplot(exp_barsO, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(exp_order, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_family <- summarize_taxonomy(exp, level = 5, report_higher_tax = F)
plot_ts_heatmap(exp_family, exp$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsF <- plot_taxa_bars(exp_family, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Exp_Family.pdf", width = 7, height = 5)
ggplot(exp_barsF, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(exp_family, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_genus <- summarize_taxonomy(exp, level = 6, report_higher_tax = F)
plot_ts_heatmap(exp_genus, exp$map_loaded, 0.005, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsG <- plot_taxa_bars(exp_genus, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Exp_Genus.pdf", width = 7, height = 5)
ggplot(exp_barsG, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(exp_genus, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_guilds <- summarize_taxonomy(exp, level = 9, report_higher_tax = F)
plot_ts_heatmap(exp_guilds, exp$map_loaded, 0, 'EstSalt', rev_taxa = T, remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsGu <- plot_taxa_bars(exp_guilds, exp$map_loaded, "EstSalt",
                             num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Exp_Guilds.pdf", width = 7, height = 5)
ggplot(exp_barsGu, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative Abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(exp_guilds, exp$map_loaded, 'EstSalt', 0.01, 'KW')

#### __Methano ####
# Summarize by family, extract methanogens
tax_sum_family_wTax <- summarize_taxonomy(exp, level = 5, report_higher_tax = T, 
                                          relative = FALSE)
methano_wTax <- tax_sum_family_wTax[grep("Methano", rownames(tax_sum_family_wTax)),]
# 16, including 4 NA

# Without higher tax, remove NAs, relative abundance
tax_sum_family <- summarize_taxonomy(exp, level = 5, report_higher_tax = FALSE, 
                                     relative = TRUE)
methano <- tax_sum_family[grep("Methano", rownames(tax_sum_family)),]
# 11
exp_barsMethano <- plot_taxa_bars(methano, exp$map_loaded, "EstSalt", 
                                  num_taxa = nrow(methano), data_only = T) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Exp_Methano.pdf", width = 7, height = 5)
ggplot(exp_barsMethano, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Family") +
  scale_fill_brewer(palette = "Paired") +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

#### __Desulfo ####
# Summarize by family, extract sulfate reducers
tax_sum_family_wTax <- summarize_taxonomy(exp, level = 5, report_higher_tax = T, 
                                          relative = FALSE)
desulfo_wTax <- tax_sum_family_wTax[grep("Desulfo", rownames(tax_sum_family_wTax)),]
# 41, including some NA. Note, some Desulfobacterota and some Firmicutes
# Too many, try order
tax_sum_order_wTax <- summarize_taxonomy(exp, level = 4, report_higher_tax = T, 
                                         relative = TRUE)
desulfo_wTax <- tax_sum_order_wTax[grep("Desulfo", rownames(tax_sum_order_wTax)),]
# 20, including 4 NA
# Trim off bacteria
rownames(desulfo_wTax) <- substring(rownames(desulfo_wTax), 11)
nb.cols <- nrow(desulfo_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

exp_barsDesulfo <- plot_taxa_bars(desulfo_wTax, exp$map_loaded, "EstSalt", 
                                  num_taxa = nrow(desulfo_wTax), data_only = T) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("Figs/Comb_Exp_Desulfo.pdf", width = 7, height = 5)
ggplot(exp_barsDesulfo, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative Abundance", fill = "Order") +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

#### __Simper ####
exp_sim <- simper(t(exp$data_loaded), 
                  exp$map_loaded$Salt)
exp_s <- summary(exp_sim)
head(exp_s$`Control_+ASW`)
exp_df1 <- head(exp_s$`Control_+ASW`, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Control/+ASW") %>%
  left_join(., exp$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Increase", "Decrease")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanSalt" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanSalt, MeanControl, CumulativeContribution)
write_xlsx(exp_df1, 
           "simper_results_comb_exp.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
exp_mp <- multipatt(t(exp$data_loaded), 
                    exp$map_loaded$Salt, 
                    func = "r.g", 
                    control = how(nperm=999),
                    max.order = 1)
exp_mp_results <- exp_mp$sign %>%
  mutate(q.value = qvalue(exp_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(exp_mp$sign)) %>%
  filter(p.value < 0.01) %>%
  filter(stat >= 0.55) %>%
  left_join(., exp$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(exp_mp_results)) {
  if (exp_mp_results$s.Control[i] == 1) {
    exp_mp_results$Group[i] <- "Control"
  }
}
for (i in 1:nrow(exp_mp_results)) {
  if (exp_mp_results$`s.+ASW`[i] == 1) {
    exp_mp_results$Group[i] <- "+ASW"
  }
}
table(exp_mp_results$Group)
exp_asv <- summarize_taxonomy(exp, level = 8, report_higher_tax = F)
exp_asv_all <- data.frame("RelAbundance" = round(rowMeans(exp_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
exp_mp_corrs <- as.data.frame(exp_mp$str) %>%
  dplyr::select(1:length(levels(exp$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% exp_mp_results$ASV) %>%
  set_names(c("Control", "+ASW", "ASV"))
# Add corrs and taxonomy
exp_mp_results <- exp_mp_results %>%
  left_join(., exp_asv_all, by = "ASV") %>%
  left_join(., exp_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(exp_mp_corrs)[1:2], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(exp_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
exp.hm.melted <- exp_mp_results %>%
  dplyr::select(taxon, names(exp_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
exp.hm <- ggplot(data = exp.hm.melted, 
                 aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(exp.hm.melted$taxon), labels = unique(exp.hm.melted$taxon),
                   limits = rev(levels(exp.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
exp.l <- get_legend(exp.hm)
exp.hm.clean <- exp.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)), size = 5),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
exp.bp.y <- ggplot(data = exp_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(exp_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("Figs/Comb_Exp_Multipatt.pdf", width = 8, height = 10)
plot_grid(exp.hm.clean, exp.bp.y, exp.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

#### End Script ####