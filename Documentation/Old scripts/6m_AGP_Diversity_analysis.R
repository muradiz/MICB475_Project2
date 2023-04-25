#load packages
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggforce)

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
meta <- read_delim("agp_sorted_metadata_6m.txt", delim="\t")
phylotree <- read.tree(file = "tree.nwk")

#rename first column of metadata
names(meta)[names(meta) == '#SampleID'] <- 'sampleid'

#filtering otu and tax files
otu_filt <- otu %>% select("#OTU ID", one_of(meta$sampleid))

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
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$sampleid
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)



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
agp_infant_6m <- phyloseq(OTU, SAMP, TAX, phylotree)

otu_table(agp_infant_6m)
sample_data(agp_infant_6m)
tax_table(agp_infant_6m)
phy_tree(agp_infant_6m)


#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(agp_infant_6m)
# colnames(sample_data(infant_6m))
get_variable(agp_infant_6m, c("agp", "agp_status")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# What is the sum of reads in each sample?
sample_sums(agp_infant_6m)
# Save the sample names of the 3 samples with the most reads
getsamps <- names(sort(sample_sums(agp_infant_6m), decreasing = TRUE)[1:3])
# filter to see taxa abundances for each sample
get_taxa(agp_infant_6m, getsamps) 


## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(agp_infant_6m)
# How many taxa do we have?
ntaxa(agp_infant_6m)
# What is the total read count for each taxa?
taxa_sums(agp_infant_6m)
# Let's find the top 3 most abundant taxa
gettaxa <- names(sort(taxa_sums(agp_infant_6m), decreasing = TRUE)[1:3] )
get_sample(agp_infant_6m, gettaxa)

######### ANALYZE ##########
# Remove non-bacterial sequences, if any
agp_infant_6m_filt <- subset_taxa(agp_infant_6m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
agp_infant_6m_filt_nolow <- filter_taxa(agp_infant_6m_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
agp_infant_6m_filt_nolow_samps <- prune_samples(sample_sums(agp_infant_6m_filt_nolow)>100, agp_infant_6m_filt_nolow)
# Remove samples where agp is na
agp_infant_6m_final <- subset_samples(agp_infant_6m_filt_nolow_samps, !is.na(agp) )


# Rarefy samples #########
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
rarecurve(t(as.data.frame(otu_table(agp_infant_6m_final))), cex=0.1)
agp_infant_6m_rare <- rarefy_even_depth(agp_infant_6m_final, rngseed = 1, sample.size = 15000)

agp_infant_6m_rare


##### Saving #####
save(agp_infant_6m_final, file="agp_infant_6m_final.RData")
save(agp_infant_6m_rare, file="agp_infant_6m_rare.RData")

#Alpha Diversity Plots
plot_richness(agp_infant_6m_rare)
# select certain alpha diversity metrics

plot_richness(agp_infant_6m_rare, measures = c("Shannon","Chao1"))
# Add ggplot layers if desired to adjust visuals
plot_richness(agp_infant_6m_rare, x = "agp_clin", measures = c("Shannon","Chao1")) +
  xlab("agp_clin") +
  geom_boxplot()

##Beta Diversity Plots
bc_dm <- distance(agp_infant_6m_rare, method="bray")
pcoa_bc <- ordinate(agp_infant_6m_rare, method="PCoA", distance=bc_dm)

plot_ordination(agp_infant_6m_rare, pcoa_bc, color = "agp", shape="agp_clin") +
  scale_color_gradient(low="darkgreen", high="lightblue") +
  labs(pch="AGP Clincal Levels", col = "AGP (g/L)")


#create a taxa summaries plot
agp_infant_6m_RA <- transform_sample_counts(agp_infant_6m_rare, function(x) x/sum(x))
plot_bar(agp_infant_6m_RA, fill="Phylum") 
agp_infant_6m_phylum <- tax_glom(agp_infant_6m_RA, taxrank = "Phylum", NArm=FALSE)
plot_bar(agp_infant_6m_phylum, fill="Phylum") +
  facet_wrap(.~agp_clin, scales = "free_x")
