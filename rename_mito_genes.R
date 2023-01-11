
# Append "MT-" to front of mitochondrial gene names and prefer using gene names 
# over ensemble IDs in our GEX data

##############################
## Load libraries and input ##
##############################

library(Seurat)
library(Rgb)
library(dplyr)
library(reshape2)

gex <- Read10X_h5("~/filtered_feature_bc_matrix.h5")

# Output from `cellranger mkgtf`. This gtf was used to make our reference.
gtf_raw <- read.gtf("~/Macaca_mulatta.Mmul_10.107.filtered_with_mt.gtf", 
                    attr = "split")

###########################
## Extract info from GTF ##
###########################

# Make a data frame of gene names and ensembl IDs
ensembl <- gtf_raw %>%
  select(gene_id, gene_name, gene_biotype) %>%
  rename(id = gene_id,
         biotype = gene_biotype) %>%
  distinct() %>%
  # If gene_name is NA use ensemble id instead
  mutate(gene_name = coalesce(gene_name, id)) %>%
  # Append "MT-" to front of mitochondrial gene names
  mutate(name = case_when(
    biotype == "Mt_tRNA" | biotype == "Mt_rRNA" ~ paste0("MT-", gene_name),
    TRUE ~ gene_name
  )) %>%
  melt(id.var = c('biotype', 'name'),
       variable.name = 'identifier') %>%
  select(value, biotype, name) %>%
  distinct()

# Remove large, unneeded object
rm(gtf_raw)
gc()

# Confirm mitochondrial genes are "labelled" correctly in our data frame
ensembl %>% filter(grepl("^MT-", name)) %>% head(3)

########################################
## Update gene names in Seurat object ##
########################################

# Takes a few minutes
for (i in 1:nrow(gex)) {
  current <-  rownames(gex)[i]
  updated <- ensembl[ensembl$value == current, ]$name
  if (length(updated) > 0) {
    rownames(gex)[i] <- updated
  } else {
    rownames(gex)[i] <- current
  }
}

# Confirm there are mitochondrial genes in our GEX data
rownames(gex)[grepl("^MT-", rownames(gex))] %>% head()

