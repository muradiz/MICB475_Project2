#required files
#anemia_metadata.txt
#taxonomy.tsv
#table_250.tsv
#tree.nwk

#loading packages 
library(tidyverse)
library(phyloseq)
library(dplyr)
library(ape)
library(vegan)
library(indicspecies)

anemia_metadata <- read_delim("anemia_metadata.txt", delim="\t")
tax <- read_delim("taxonomy.tsv", delim="\t")

#only keep healthy infant samples
filtered_meta_healthy <- filter(anemia_metadata, anemia == "normal", parasites == "N")
# filter out samples with no CRP value
filtered_meta_healthy_crp <- filtered_meta_healthy[!is.na(as.numeric(filtered_meta_healthy$crp)),]
#make metadata file for only 6 month infant samples
meta_6m <-  filter(filtered_meta_healthy_crp, age_months == 6)
#make metadata file for only 12 month infant samples
meta_12m <- filter(filtered_meta_healthy_crp, age_months == 12)

#determine median value for crp at 6m
median(meta_6m$crp)
#determine median value for crp at 12m
median(meta_12m$crp)

#create column based on median crp value - 6m
meta_6m$crp_median <- ifelse(meta_6m$crp >= 0.54,"Above", "Below")
#create column based on median crp value - 12m
meta_12m$crp_median <- ifelse(meta_12m$crp >= 0.95,"Above", "Below")

# Select columns
meta_6m_filt <- select(meta_6m, "#SampleID", "host_subject_id", "sex", "agp", "agp_status", "crp", "crp_status", "infection_status", "crp_median")
meta_12m_filt <- select(meta_12m, "#SampleID", "host_subject_id", "sex", "agp", "agp_status", "crp", "crp_status", "infection_status", "crp_median")

# Save as RData files
save(meta_6m_filt, file = "sorted_metadata_6m.RData")
save(meta_12m_filt, file = "sorted_metadata_12m.RData")

# Save as txt files - NOT WORKING
write.table(meta_12m_filt, file = "sorted_metadata_12m.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(meta_6m_filt, file = "sorted_metadata_6m.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Group by and summarize
# group_by(meta_6m_filt, crp_median)

#load in metadata and taxonomy tsvs
otu <- read_delim(file="table_250.tsv", delim = "\t", skip=1)
tax <- read_delim(file = "taxonomy.tsv", delim="\t")
meta_12m_crp <- read_delim("sorted_metadata_12m.txt", delim="\t")
meta_6m_crp <- read_delim("sorted_metadata_6m.txt", delim="\t")
phylotree <- read.tree(file = "tree.nwk")

#rename first column of metadata
names(meta_12m_crp)[names(meta_12m_crp) == '#SampleID'] <- 'sampleid'
names(meta_6m_crp)[names(meta_6m_crp) == '#SampleID'] <- 'sampleid'

#filtering otu and tax files
otu_filt_12m <- otu %>% select("#OTU ID", one_of(meta_12m_crp$sampleid))
otu_filt_6m <- otu %>% select("#OTU ID", one_of(meta_6m_crp$sampleid))

tax_sep <- tax %>%
  separate(col=Taxon, sep=";", into = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species"))

#### Format OTU table ####
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)


#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df_12m_crp <- as.data.frame(meta_12m_crp[,-1])
samp_df_6m_crp <- as.data.frame(meta_6m_crp[,-1])

# Make sampleids the rownames
rownames(samp_df_12m_crp)<- meta_12m_crp$sampleid
rownames(samp_df_6m_crp) <- meta_6m_crp$sampleid

# Make phyloseq sample data with sample_data() function
SAMP_12m_crp <- sample_data(samp_df_12m_crp)
SAMP_6m_crp <- sample_data(samp_df_6m_crp)
class(SAMP_12m_crp)
class(SAMP_6m_crp)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence) %>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$"Feature ID"
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

class(phylotree)

#### Create phyloseq object ####
# Merge all into a phyloseq object
infant_12m_crp <- phyloseq(OTU, SAMP_12m_crp, TAX, phylotree)
infant_6m_crp <- phyloseq(OTU, SAMP_6m_crp, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
#otu_table(infant_12m)
#sample_data(infant_12m)
#tax_table(infant_12m)
#phy_tree(infant_12m)

### Indicator Species Analysis

#Setting random seed as 1 
set.seed(1)

crp_genus_12m <- tax_glom(infant_12m_crp, "Genus", NArm= FALSE)
crp_genus_RA_12m <- transform_sample_counts(crp_genus_12m, fun=function(x) x/sum(x))

#ISA CRP 12 months
isa_crp_12m <- multipatt(t(otu_table(crp_genus_RA_12m)), cluster = sample_data(crp_genus_RA_12m)$crp_median)
summary(isa_crp_12m)
taxtable_12m_crp <- tax_table(infant_12m_crp) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_crp_12m$sign %>% rownames_to_column(var="ASV") %>%
  left_join(taxtable_12m_crp) %>% filter(p.value<0.05) %>% View()

#ISA CRP 6 months
crp_genus_6m <- tax_glom(infant_6m_crp, "Genus", NArm= FALSE)
crp_genus_RA_6m <- transform_sample_counts(crp_genus_6m, fun=function(x) x/sum(x))

isa_crp_6m <- multipatt(t(otu_table(crp_genus_RA_6m)), cluster = sample_data(crp_genus_RA_6m)$crp_median)
summary(isa_crp_6m)
taxtable_6m_crp <- tax_table(infant_6m_crp) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_crp_6m$sign %>% rownames_to_column(var="ASV") %>%
  left_join(taxtable_6m_crp) %>% filter(p.value<0.05) %>% View()

