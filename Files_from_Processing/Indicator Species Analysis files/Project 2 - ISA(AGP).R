#required files
#anemia_metadata.txt
#taxonomy.tsv
#table_250.tsv
#tree.nwk

#load packages
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(dplyr)

anemia_metadata <- read_delim("anemia_metadata.txt", delim="\t")
tax <- read_delim("taxonomy.tsv", delim="\t")

#only keep healthy infant samples
filtered_meta_healthy <- filter(anemia_metadata, anemia == "normal", parasites == "N")
# filter out samples with no agp value
filtered_meta_healthy_agp <- filtered_meta_healthy[!is.na(as.numeric(filtered_meta_healthy$agp)),]

#make metadata file for only 6 month infant samples
meta_6m_agp <-  filter(filtered_meta_healthy_agp, age_months == 6)
#make metadata file for only 12 month infant samples
meta_12m_agp <- filter(filtered_meta_healthy_agp, age_months == 12)

##create column based on clinical AGP value - 6m
meta_6m_agp$agp_clin <- ifelse(meta_6m_agp$agp >= 1,"Above", "Below")
#create column based on clinical AGP value - 12m
meta_12m_agp$agp_clin <- ifelse(meta_12m_agp$agp >= 1,"Above", "Below")

# Select columns
meta_6m_filt_agp <- select(meta_6m_agp, "#SampleID", "host_subject_id", "sex", "agp", "agp_status", "agp_clin", "infection_status")
meta_12m_filt_agp <- select(meta_12m_agp, "#SampleID", "host_subject_id", "sex", "agp", "agp_status", "agp_clin", "infection_status")

#save files
save(meta_6m_filt_agp, file = "agp_sorted_metadata_6m.RData")
save(meta_12m_filt_agp, file = "agp_sorted_metadata_12m.RData")

# Save as txt files - NOT WORKING
write.table(meta_12m_filt_agp, file = "agp_sorted_metadata_12m.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(meta_6m_filt_agp, file = "agp_sorted_metadata_6m.txt", sep="\t", quote=FALSE, row.names=FALSE)


#load in metadata and taxonomy tsvs
otu <- read_delim(file="table_250.tsv", delim = "\t", skip=1)
tax <- read_delim(file = "taxonomy.tsv", delim="\t")
meta_12m <- read_delim("agp_sorted_metadata_12m.txt", delim="\t")
meta_6m <- read_delim("agp_sorted_metadata_6m.txt", delim="\t")
phylotree <- read.tree(file = "tree.nwk")

#rename first column of metadata for 12 months and 6 months
names(meta_12m)[names(meta_12m) == '#SampleID'] <- 'sampleid'
names(meta_6m)[names(meta_6m) == '#SampleID'] <- 'sampleid'

#filtering tax files
tax_sep <- tax %>%
  separate(col=Taxon, sep=";", into = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species"))

#filtering OTU for 12 months and 6 months 
otu_filt_12m <- otu %>% select("#OTU ID", one_of(meta_12m$sampleid))
otu_filt_6m <- otu %>% select("#OTU ID", one_of(meta_6m$sampleid))


#### Format OTU table ####

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame for 12 months and 6 months
samp_df_12m <- as.data.frame(meta_12m[,-1])
samp_df_6m <- as.data.frame(meta_6m[,-1])

# Make sampleids the rownames for 12 months and 6 months
rownames(samp_df_12m)<- meta_12m$sampleid
rownames(samp_df_6m) <- meta_6m$sampleid

# Make phyloseq sample data with sample_data() function
SAMP_12m <- sample_data(samp_df_12m)
SAMP_6m <- sample_data(samp_df_6m)
class(SAMP_12m)
class(SAMP_6m)

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
agp_infant_12m <- phyloseq(OTU, SAMP_12m, TAX, phylotree)
agp_infant_6m <- phyloseq(OTU, SAMP_6m, TAX, phylotree)

otu_table(agp_infant_12m)
sample_data(agp_infant_12m)
tax_table(agp_infant_12m)
phy_tree(agp_infant_12m)

otu_table(agp_infant_6m)
sample_data(agp_infant_6m)
tax_table(agp_infant_6m)
phy_tree(agp_infant_6m)

library(indicspecies)

#Setting random seed as 1 
set.seed(1)

agp_genus_12m <- tax_glom(agp_infant_12m, "Genus", NArm= FALSE)
agp_genus_RA_12m <- transform_sample_counts(agp_genus_12m, fun=function(x) x/sum(x))

#ISA AGP 12 months
isa_agp_12m <- multipatt(t(otu_table(agp_genus_RA_12m)), cluster = sample_data(agp_genus_RA_12m)$agp_clin)
summary(isa_agp_12m)
taxtable_12m <- tax_table(agp_infant_12m) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_agp_12m$sign %>% rownames_to_column(var="ASV") %>%
 left_join(taxtable_12m) %>% filter(p.value<0.05) %>% View()

#ISA AGP 6 months
agp_genus_6m <- tax_glom(agp_infant_6m, "Genus", NArm= FALSE)
agp_genus_RA_6m <- transform_sample_counts(agp_genus_6m, fun=function(x) x/sum(x))

isa_agp_6m <- multipatt(t(otu_table(agp_genus_RA_6m)), cluster = sample_data(agp_genus_RA_6m)$agp_clin)
summary(isa_agp_6m)
taxtable_6m <- tax_table(agp_infant_6m) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_agp_6m$sign %>% rownames_to_column(var="ASV") %>%
  left_join(taxtable_6m) %>% filter(p.value<0.05) %>% View()

