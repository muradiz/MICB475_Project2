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
group_by(meta_6m_filt, crp_median)

#load in metadata and taxonomy tsvs
otu <- read_delim(file="table_250.tsv", delim = "\t", skip=1)
tax <- read_delim(file = "taxonomy.tsv", delim="\t")
meta <- read_delim("sorted_metadata_6m.txt", delim="\t")
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
infant_6m <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(infant_6m)
sample_data(infant_6m)
tax_table(infant_6m)
phy_tree(infant_6m)


#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(infant_6m)
# colnames(sample_data(infant_12m))
get_variable(infant_6m, c("crp", "crp_status")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# What is the sum of reads in each sample?
sample_sums(infant_6m)
# Save the sample names of the 3 samples with the most reads
getsamps <- names(sort(sample_sums(infant_6m), decreasing = TRUE)[1:3])
# filter to see taxa abundances for each sample
get_taxa(infant_6m, getsamps) 


## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(infant_6m)
# How many taxa do we have?
ntaxa(infant_6m)
# What is the total read count for each taxa?
taxa_sums(infant_6m)
# Let's find the top 3 most abundant taxa
gettaxa <- names(sort(taxa_sums(infant_6m), decreasing = TRUE)[1:3] )
get_sample(infant_6m, gettaxa)


######### ANALYZE ##########
# Remove non-bacterial sequences, if any
infant_6m_filt <- subset_taxa(infant_6m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
infant_6m_filt_nolow <- filter_taxa(infant_6m_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
infant_6m_filt_nolow_samps <- prune_samples(sample_sums(infant_6m_filt_nolow)>100, infant_6m_filt_nolow)
# Remove samples where ph is na
infant_6m_final <- subset_samples(infant_6m_filt_nolow_samps, !is.na(crp) )


# Rarefy samples #########
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
rarecurve(t(as.data.frame(otu_table(infant_6m_final))), cex=0.1)
infant_6m_rare <- rarefy_even_depth(infant_6m_final, rngseed = 1, sample.size = 15000)

infant_6m_rare


##### Saving #####
save(infant_6m_final, file="infant_6m_final.RData")
save(infant_6m_rare, file="infant_6m_rare.RData")


#Alpha Diversity Plots
plot_richness(infant_6m_rare)
# select certain alpha diversity metrics

plot_richness(infant_6m_rare, measures = c("Shannon","Chao1"))
# Add ggplot layers if desired to adjust visuals
plot_richness(infant_6m_rare, x = "crp_median", measures = c("Shannon","Chao1")) +
  xlab("crp median") +
  geom_boxplot()

##Beta Diversity Plots
bc_dm <- distance(infant_6m_rare, method="bray")
pcoa_bc <- ordinate(infant_6m_rare, method="PCoA", distance=bc_dm)

plot_ordination(infant_6m_rare, pcoa_bc, color = "crp", shape="crp_median") +
  scale_color_gradient(low="darkgreen", high="lightblue") +
  labs(pch="CRP status", col = "CRP (mg/L)")

#create a taxa summaries plot
infant_6m_RA <- transform_sample_counts(infant_6m_rare, function(x) x/sum(x))
plot_bar(infant_6m_RA, fill="Phylum") 
infant_6m_phylum <- tax_glom(infant_6m_RA, taxrank = "Phylum", NArm=FALSE)
plot_bar(infant_6m_phylum, fill="Phylum") +
  facet_wrap(.~crp_median, scales = "free_x")
