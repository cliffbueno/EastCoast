# NTI

# But, what we really should do is NTI for each experiment.
saveRDS(sf, "sf.rds")
saveRDS(sc, "sc.rds")
saveRDS(detra, "detra.rds")
saveRDS(deinc, "deinc.rds")
saveRDS(nc, "nc.rds")

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
saveRDS(sf_NTI, "sf_NTI.rds")

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
saveRDS(sc_NTI, "sc_NTI.rds")

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
saveRDS(detra_NTI, "detra_NTI.rds")

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
saveRDS(deinc_NTI, "deinc_NTI.rds")

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
saveRDS(nc_NTI, "nc_NTI.rds")