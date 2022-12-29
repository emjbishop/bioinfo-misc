
# Look at HTO counts and number of doublet/singlet/negative cells post demux

##############################
## Load libraries and input ##
##############################

library(tidyverse)
library(Seurat)
library(hdf5r)
library(stringr)
library(knitr)
library(Rgb)
library(reshape2)

my_h5 <- "" 
ref_gtf <- ""

#############
## Load h5 ##
#############
h5 <- Read10X_h5(my_h5)

# Raw counts for each ADT and HTO
rowSums(h5$`Antibody Capture`) %>%
  as.data.frame() %>%
  kable()

# How many cells have expression above 0
rowSums(h5$`Antibody Capture`> 0) %>%
  as.data.frame() %>%
  tail(8) %>%
  kable()

# How many cells have expression above 100
rowSums(h5$`Antibody Capture`> 100) %>%
  as.data.frame() %>%
  tail(8) %>%
  kable()

# Take a look at three cells' HTO counts
as.data.frame(h5$`Antibody Capture`)[21:28, 1:3] 

##################
## Rename genes ##
##################
#  - Use ensembl gene name (not id) when available
#  - Append "MT-" to front of mitochondrial gene names

gtf_raw <- read.gtf(ref_gtf, attr = "split")

# Get gene names
ensembl <- gtf_raw %>%
  select(gene_id, gene_name, gene_biotype) %>%
  rename(id = gene_id,
         biotype = gene_biotype) %>%
  distinct() %>%
  # If gene_name is NA use id
  mutate(gene_name = coalesce(gene_name, id)) %>%
  # Append "MT-"
  mutate(name = case_when(
    biotype == "Mt_tRNA" | biotype == "Mt_rRNA" ~ paste0("MT-", gene_name),
    TRUE                                        ~ gene_name
  )) %>%
  melt(id.var = c('biotype', 'name'),
       variable.name = 'identifier') %>%
  select(value, biotype, name) %>%
  distinct()

# Remove unneeded object
rm(gtf_raw)
gc()

# Update names in Seurat object (takes a minute)
for (i in 1:nrow(h5$`Gene Expression`)) {
  current <-  rownames(h5$`Gene Expression`)[i]
  updated <- ensembl[ensembl$value == current, ]$name
  if (length(updated) > 0) {
    rownames(h5$`Gene Expression`)[i] <- updated
  } else {
    rownames(h5$`Gene Expression`)[i] <- current
  }
}

# Confirm mitochondrial genes present
# TODO: just have this output T/F
rownames(h5$`Gene Expression`)[grepl("^MT-", rownames(h5$`Gene Expression`))] %>% head()

#################################
## Separate HTO vs rest of ADT ##
#################################

ADT <- !grepl("HTO", h5$`Antibody Capture`@Dimnames[[1]])
ADT_counts <- h5$`Antibody Capture`[ADT,]

HTO <- grepl("HTO", h5$`Antibody Capture`@Dimnames[[1]])
HTO_counts <- h5$`Antibody Capture`[HTO,]

##########################
## Set up Seurat object ##
##########################

dat <- CreateSeuratObject(counts = h5$`Gene Expression`) #, project = "HAARVIVAC_10X_S1", min.cells = 3)
dat[["ADT"]] <-  CreateAssayObject(counts = ADT_counts)
dat[["HTO"]] <-  CreateAssayObject(counts = HTO_counts)

# Number of cells
length(dat$orig.ident)

########################
## Demultiplex by HTO ##
########################

# Normalize
dat_norm <- NormalizeData(dat, assay = "HTO", normalization.method = "CLR")

# HTODemux
dat_demux <- HTODemux(dat_norm, assay = "HTO", positive.quantile = 0.99)

# MULTISeqDemux
dat_multi_demux <- MULTIseqDemux(dat_norm, assay = "HTO", autoThresh = TRUE)
# Assign global classification
i1 <- grepl("HTO", dat_multi_demux@meta.data$MULTI_ID)
i2 <- grepl("Doublet", dat_multi_demux@meta.data$MULTI_ID)
i3 <- grepl("Negative", dat_multi_demux@meta.data$MULTI_ID)
dat_multi_demux@meta.data$MULTI_ID_global <- NULL
dat_multi_demux@meta.data$MULTI_ID_global[i1] <- "Singlet"
dat_multi_demux@meta.data$MULTI_ID_global[i2] <- "Doublet"
dat_multi_demux@meta.data$MULTI_ID_global[i3] <- "Negative"
dat_multi_demux@meta.data$MULTI_ID_global <- factor(dat_multi_demux@meta.data$MULTI_ID_global, 
                                                    levels = c("Doublet", "Negative", "Singlet"))

# Compare the two:
table(dat_demux$HTO_classification.global)  # HTODemux
table(dat_multi_demux@meta.data$MULTI_ID_global)  # MULTISeqDemux

###############
## Summarize ##
###############

# Compare demux methods
table(dat_demux$HTO_classification.global)  # HTODemux
table(dat_multi_demux@meta.data$MULTI_ID_global)  # MULTISeqDemux

# How many cells have expression above 0
rowSums(h5$`Antibody Capture`> 0) %>%
  as.data.frame() %>%
  tail(8) %>%
  kable()

# How many cells have expression above 100
rowSums(h5$`Antibody Capture`> 100) %>%
  as.data.frame() %>%
  tail(8) %>%
  kable()

# Take a look at three cells' HTO counts
as.data.frame(h5$`Antibody Capture`)[21:28, 1:3] 

