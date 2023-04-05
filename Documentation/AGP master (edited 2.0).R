#load packages
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggforce)
library(DESeq2)
library(indicspecies)

anemia_metadata <- read_delim("Files_from_Processing/anemia_metadata.txt", delim="\t")
tax <- read_delim("Files_from_Processing/Exported_files/taxonomy.tsv", delim="\t")

#only keep healthy infant samples
filtered_meta_healthy <- filter(anemia_metadata, anemia == "normal", parasites == "N")
# filter out samples with no agp value
filtered_meta_healthy_agp <- filtered_meta_healthy[!is.na(as.numeric(filtered_meta_healthy$agp)),]

#make metadata file for only 6 month infant samples
meta_6m_agp <-  filter(filtered_meta_healthy_agp, age_months == 6)
#make metadata file for only month infant samples
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
otu <- read_delim(file="Files_from_Processing/Exported_files/table_250.tsv", delim = "\t", skip=1)
tax <- read_delim(file = "Files_from_Processing/Exported_files/taxonomy.tsv", delim="\t")
meta_6m <- read_delim("agp_sorted_metadata_6m.txt", delim="\t")
meta_12m <- read_delim("agp_sorted_metadata_12m.txt", delim="\t")
phylotree <- read.tree(file = "Files_from_Processing/Exported_files/tree.nwk")

#rename first column of metadata
names(meta_6m)[names(meta_6m) == '#SampleID'] <- 'sampleid'
names(meta_12m)[names(meta_12m) == '#SampleID'] <- 'sampleid'

#filtering otu and tax files
otu_filt_6m <- otu %>% select("#OTU ID", one_of(meta_6m$sampleid))
otu_filt_12m <- otu %>% select("#OTU ID", one_of(meta_12m$sampleid))

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
samp_df_6m <- as.data.frame(meta_6m[,-1])
samp_df_12m <- as.data.frame(meta_12m[,-1])

# Make sampleids the rownames
rownames(samp_df_6m)<- meta_6m$sampleid
rownames(samp_df_12m) <- meta_12m$sampleid

# Make phyloseq sample data with sample_data() function
SAMP_6m <- sample_data(samp_df_6m)
SAMP_12m <- sample_data(samp_df_12m)
class(SAMP_6m)
class(SAMP_12m)

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
agp_infant_6m <- phyloseq(OTU, SAMP_6m, TAX, phylotree)
agp_infant_12m <- phyloseq(OTU, SAMP_12m, TAX, phylotree)

otu_table(agp_infant_6m)
sample_data(agp_infant_6m)
tax_table(agp_infant_6m)
phy_tree(agp_infant_6m)

otu_table(agp_infant_12m)
sample_data(agp_infant_12m)
tax_table(agp_infant_12m)
phy_tree(agp_infant_12m)

#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(agp_infant_6m)
sample_variables(agp_infant_12m)
# colnames(sample_data(infant_6m))
get_variable(agp_infant_6m, c("agp", "agp_status")) # equivalent to "select" in tidyverse
get_variable(agp_infant_12m, c("agp", "agp_status"))
## Let's say we want to filter OTU table by sample. 
# What is the sum of reads in each sample?
sample_sums(agp_infant_6m)
sample_sums(agp_infant_12m)
# Save the sample names of the 3 samples with the most reads
getsamps_6m <- names(sort(sample_sums(agp_infant_6m), decreasing = TRUE)[1:3])
getsamps_12m <- names(sort(sample_sums(agp_infant_12m), decreasing = TRUE)[1:3])

# filter to see taxa abundances for each sample
get_taxa(agp_infant_6m, getsamps_6m) 
get_taxa(agp_infant_12m, getsamps_12m)

## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(agp_infant_6m)
taxa_names(agp_infant_12m)

# How many taxa do we have?
ntaxa(agp_infant_6m)
ntaxa(agp_infant_12m)

# What is the total read count for each taxa?
taxa_sums(agp_infant_6m)
taxa_sums(agp_infant_12m)

# Let's find the top 3 most abundant taxa
gettaxa_6m <- names(sort(taxa_sums(agp_infant_6m), decreasing = TRUE)[1:3] )
get_sample(agp_infant_6m, gettaxa_6m)

gettaxa_12m <- names(sort(taxa_sums(agp_infant_6m), decreasing = TRUE)[1:3] )
get_sample(agp_infant_12m, gettaxa_12m)
######### ANALYZE ##########
# Remove non-bacterial sequences, if any
agp_infant_6m_filt <- subset_taxa(agp_infant_6m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
agp_infant_12m_filt <- subset_taxa(agp_infant_12m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Remove ASVs that have less than 5 counts total
agp_infant_6m_filt_nolow <- filter_taxa(agp_infant_6m_filt, function(x) sum(x)>5, prune = TRUE)
agp_infant_12m_filt_nolow <- filter_taxa(agp_infant_12m_filt, function(x) sum(x)>5, prune = TRUE)

# Genus Level
agp_infant_6m_genus <- tax_glom(agp_infant_6m_filt_nolow, "Genus")
agp_infant_12m_genus <- tax_glom(agp_infant_12m_filt_nolow, "Genus")

# Remove samples with less than 100 reads
agp_infant_6m_filt_nolow_samps <- prune_samples(sample_sums(agp_infant_6m_genus)>100, agp_infant_6m_filt_nolow)
agp_infant_12m_filt_nolow_samps <- prune_samples(sample_sums(agp_infant_12m_genus)>100, agp_infant_12m_filt_nolow)

# Remove samples where agp is na
agp_infant_6m_final <- subset_samples(agp_infant_6m_filt_nolow_samps, !is.na(agp) )
agp_infant_12m_final <- subset_samples(agp_infant_12m_filt_nolow_samps, !is.na(agp) )

# Rarefy samples #########
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
#6 months
rarecurve(t(as.data.frame(otu_table(agp_infant_6m_final))), cex=0.1)
agp_infant_6m_rare <- rarefy_even_depth(agp_infant_6m_final, rngseed = 1, sample.size = 15000)

agp_infant_6m_rare

#12 months
rarecurve(t(as.data.frame(otu_table(agp_infant_12m_final))), cex=0.1)
agp_infant_12m_rare <- rarefy_even_depth(agp_infant_12m_final, rngseed = 1, sample.size = 15000)

agp_infant_12m_rare

##### Saving #####
#6 months
save(agp_infant_6m_final, file="agp_infant_6m_final.RData")
save(agp_infant_6m_rare, file="agp_infant_6m_rare.RData")

#12 months
save(agp_infant_12m_final, file="agp_infant_12m_final.RData")
save(agp_infant_12m_rare, file="agp_infant_12m_rare.RData")

########################### DIVERSITY ANALYSIS ################################
#6 months
#Alpha Diversity Plots
plot_richness(agp_infant_6m_rare)
# select certain alpha diversity metrics

plot_richness(agp_infant_6m_rare, measures = c("Shannon","Chao1"))
# Add ggplot layers if desired to adjust visuals
plot_richness(agp_infant_6m_rare, x = "agp_clin", measures = c("Shannon","Chao1")) +
  xlab("agp_clin") +
  geom_boxplot()

##Beta Diversity Plots
bc_dm_6m <- phyloseq::distance(agp_infant_6m_rare, method="bray")
pcoa_bc_6m <- ordinate(agp_infant_6m_rare, method="PCoA", distance=bc_dm_6m)

plot_ordination(agp_infant_6m_rare, pcoa_bc_6m, color = "agp", shape="agp_clin") +
  scale_color_gradient(low="darkgreen", high="lightblue") +
  labs(pch="AGP Clincal Levels", col = "AGP (g/L)")

#create a taxa summaries plot
agp_infant_6m_RA <- transform_sample_counts(agp_infant_6m_rare, function(x) x/sum(x))
plot_bar(agp_infant_6m_RA, fill="Phylum") 
agp_infant_6m_phylum <- tax_glom(agp_infant_6m_RA, taxrank = "Phylum", NArm=FALSE)
plot_bar(agp_infant_6m_phylum, fill="Phylum") +
  facet_wrap(.~agp_clin, scales = "free_x")

#12 months
#Alpha Diversity Plots
plot_richness(agp_infant_12m_rare)
# select certain alpha diversity metrics

plot_richness(agp_infant_12m_rare, measures = c("Shannon","Chao1"))
# Add ggplot layers if desired to adjust visuals
plot_richness(agp_infant_12m_rare, x = "agp_clin", measures = c("Shannon","Chao1")) +
  xlab("agp_clin") +
  geom_boxplot()

##Beta Diversity Plots
bc_dm_12m <- phyloseq::distance(agp_infant_12m_rare, method="bray")
pcoa_bc_12m <- ordinate(agp_infant_12m_rare, method="PCoA", distance=bc_dm_12m)

plot_ordination(agp_infant_12m_rare, pcoa_bc_12m, color = "agp", shape="agp_clin") +
  scale_color_gradient(low="darkgreen", high="lightblue") +
  labs(pch="AGP Clincal Levels", col = "AGP (g/L)")

#create a taxa summaries plot
agp_infant_12m_RA <- transform_sample_counts(agp_infant_12m_rare, function(x) x/sum(x))
plot_bar(agp_infant_12m_RA, fill="Phylum") 
agp_infant_12m_phylum <- tax_glom(agp_infant_12m_RA, taxrank = "Phylum", NArm=FALSE)
plot_bar(agp_infant_12m_phylum, fill="Phylum") +
  facet_wrap(.~agp_clin, scales = "free_x")

################### STATISTICAL ANALYSIS #########################
#6 months
#load files
load("agp_infant_6m_rare.RData")

#t-test analysis
######## Comparison of two means with t-test (Parametric) ##########
# Let's do very simple plot with t-test
plot_richness(agp_infant_6m_rare, x = "agp_clin", measures="Shannon")
# Need to extract information
alphadiv_6m <- estimate_richness(agp_infant_6m_rare)
samp_dat_6m <- sample_data(agp_infant_6m_rare)
samp_dat_wdiv_6m <- data.frame(samp_dat_6m, alphadiv_6m)
# These are equivalent:
# t.test()
t.test(samp_dat_wdiv_6m$Shannon ~ samp_dat_wdiv_6m$agp_clin)
t.test(Shannon ~ agp_clin, data=samp_dat_wdiv_6m)

# Note: you can set variances to be equal for a "classic" t-test
t.test(Shannon ~ agp_clin, data=samp_dat_wdiv_6m, var.equal=TRUE)


#### Microbial count data is generally NON-NORMAL ####
# In fact, it is even more complex because microbial data is usually in RELATIVE ABUNDANCE
allCounts_6m <- as.vector(otu_table(agp_infant_6m_rare))
allCounts_6m <- allCounts_6m[allCounts_6m>0]
hist(allCounts_6m)
hist(log(allCounts_6m))

#remove 0 values
samp_dat_wdiv_6m %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=agp_clin, y=Shannon))+
  geom_point() 

#check distribution:
ggplot(samp_dat_wdiv_6m) +
  geom_histogram(aes(x=Shannon), bins=25)

#logging everything to deal with non-parametric distribution
# (1) Transform your data (usually with a log function)
ggplot(samp_dat_wdiv_6m) +
  geom_histogram(aes(x=log(Shannon)), bins=25)
t.test(log(Shannon) ~ agp_clin, data=samp_dat_wdiv_6m)
# Let's see what transformed data looks like:
samp_dat_wdiv_6m %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=agp_clin, y=log(Shannon)))+
  geom_boxplot() +
  geom_jitter()

#Wilcoxon Rank Sum Test
wilcox.test(Shannon ~ agp_clin, data=samp_dat_wdiv_6m)
wilcox.test(log(Shannon) ~ agp_clin, data=samp_dat_wdiv_6m)

#PERMANOVA Analysis
##load data
samp_dat_wdiv_6m <- data.frame(sample_data(agp_infant_6m_rare), estimate_richness(agp_infant_6m_rare))

### PERMANOVA (Permutational ANOVA) ####
# non-parametric version of ANOVA
# Takes a distance matrix, which can be calculated with any kind of metric you want
# e.g. Bray, Jaccard, Unifrac
# Need the package, "vegan"
# Use phyloseq to calculate weighted Unifrac distance matrix
?UniFrac
dm_unifrac_6m <- UniFrac(agp_infant_6m_rare, weighted=TRUE)
?adonis2
adonis2(dm_unifrac_6m ~ agp_clin, data=samp_dat_wdiv_6m)

# Also use other metrics: for example, the vegan package includes bray and jaccard
dm_bray_6m <- vegdist(t(otu_table(agp_infant_6m_rare)), method="bray")
adonis2(dm_bray_6m ~ agp_clin, data=samp_dat_wdiv_6m)

dm_jaccard_6m <- vegdist(t(otu_table(agp_infant_6m_rare)), method="jaccard")
adonis2(dm_jaccard_6m ~ agp_clin, data=samp_dat_wdiv_6m)

#12 months
#load files
load("agp_infant_12m_rare.RData")

#t-test analysis
######## Comparison of two means with t-test (Parametric) ##########
# Let's do very simple plot with t-test
plot_richness(agp_infant_12m_rare, x = "agp_clin", measures="Shannon")
# Need to extract information
alphadiv_12m <- estimate_richness(agp_infant_12m_rare)
samp_dat_12m <- sample_data(agp_infant_12m_rare)
samp_dat_wdiv_12m <- data.frame(samp_dat_12m, alphadiv_12m)
# These are equivalent:
# t.test()
t.test(samp_dat_wdiv_12m$Shannon ~ samp_dat_wdiv_12m$agp_clin)
t.test(Shannon ~ agp_clin, data=samp_dat_wdiv_12m)

# Note: you can set variances to be equal for a "classic" t-test
t.test(Shannon ~ agp_clin, data=samp_dat_wdiv_12m, var.equal=TRUE)

#### Microbial count data is generally NON-NORMAL ####
# In fact, it is even more complex because microbial data is usually in RELATIVE ABUNDANCE
allCounts_12m <- as.vector(otu_table(agp_infant_12m_rare))
allCounts_12m <- allCounts_12m[allCounts_12m>0]
hist(allCounts_12m)
hist(log(allCounts_12m))

#remove 0 values
samp_dat_wdiv_12m %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=agp_clin, y=Shannon))+
  geom_point() 

#check distribution:
ggplot(samp_dat_wdiv_12m) +
  geom_histogram(aes(x=Shannon), bins=25)

#logging everything to deal with non-parametric distribution
# (1) Transform your data (usually with a log function)
ggplot(samp_dat_wdiv_12m) +
  geom_histogram(aes(x=log(Shannon)), bins=25)
t.test(log(Shannon) ~ agp_clin, data=samp_dat_wdiv_12m)
# Let's see what transformed data looks like:
samp_dat_wdiv_12m %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=agp_clin, y=log(Shannon)))+
  geom_boxplot() +
  geom_jitter()

#Wilcoxon Rank Sum Test
wilcox.test(Shannon ~ agp_clin, data=samp_dat_wdiv_12m)
wilcox.test(log(Shannon) ~ agp_clin, data=samp_dat_wdiv_12m)

#PERMANOVA Analysis
##load data
samp_dat_wdiv_12m <- data.frame(sample_data(agp_infant_12m_rare), estimate_richness(agp_infant_12m_rare))

### PERMANOVA (Permutational ANOVA) ####
# non-parametric version of ANOVA
# Takes a distance matrix, which can be calculated with any kind of metric you want
# e.g. Bray, Jaccard, Unifrac
# Need the package, "vegan"
# Use phyloseq to calculate weighted Unifrac distance matrix
?UniFrac
dm_unifrac_12m <- UniFrac(agp_infant_12m_rare, weighted=TRUE)
?adonis2
adonis2(dm_unifrac_12m ~ agp_clin, data=samp_dat_wdiv_12m)

# Also use other metrics: for example, the vegan package includes bray and jaccard
dm_bray_12m <- vegdist(t(otu_table(agp_infant_12m_rare)), method="bray")
adonis2(dm_bray_12m ~ agp_clin, data=samp_dat_wdiv_12m)

dm_jaccard_12m <- vegdist(t(otu_table(agp_infant_12m_rare)), method="jaccard")
adonis2(dm_jaccard_12m ~ agp_clin, data=samp_dat_wdiv_12m)

##################### DESEQ ANALYSIS #############################
#6 months
#setting random seed 
set.seed(1)

#Using filtered data so taxglom() at GENUS level
agp_infant_6m_genus

#DESEq (No need to rarefaction [use infant_6m_final])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
infant_6m_plus1 <- transform_sample_counts(agp_infant_6m_genus, function(x) x+1)
infant_6m_deseq <- phyloseq_to_deseq2(infant_6m_plus1, ~agp_clin)
DESEQ_infant_6m <- DESeq(infant_6m_deseq)

#Viewing Results
res_6m <- results(DESEQ_infant_6m, tidy=TRUE)
View(res_6m)

#Filtering out the ASV results with N.A. values for p adjusted values
filtered_res_6m <- drop_na(res_6m, padj)
view(filtered_res_6m)

#### Visualizing the results of DESEq ####
#Chose p value <0.05, and significant fold change >2

##Volcano Plot##
vol_plot_6m <- ggplot(filtered_res_6m, aes(x=log2FoldChange, y=-log(padj))) +
  geom_point()
vol_plot_sig_6m <- filtered_res_6m |>
  mutate(significant = padj<0.05 & abs(log2FoldChange) >2) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), col= significant)) +
  geom_point()

vol_plot_6m            
vol_plot_sig_6m

#Significant ASVs table
sigASVs_6m <- filtered_res_6m |>
  filter(padj <0.05 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

save(sigASVs_6m, file= "Files_from_Processing/DESeq files/DESEq_sigASVs_6m")

view(sigASVs_6m)

#BarPlot
#removed any Genus that has NA

sigASVs_vec_6m <- sigASVs_6m |>
  pull(ASV)

infant_6m_DESeq <- prune_taxa(sigASVs_vec_6m, agp_infant_6m_genus)

sigASVs_6m <- tax_table(infant_6m_DESeq) |>
  as.data.frame() |>
  rownames_to_column(var="ASV") |>
  right_join(sigASVs_6m) |>
  arrange(log2FoldChange) |>
  mutate(Genus = make.unique(Genus)) |>
  mutate(Genus = factor(Genus, levels= unique(Genus)))|>
  drop_na(Genus) |>
  filter(Genus != "NA.2" & Genus != "NA.1")

view(sigASVs_6m)

bar_plot_6m <- ggplot(sigASVs_6m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat ="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
bar_plot_6m

#12 months
#setting random seed 
set.seed(1)

#Using filtered data so taxglom() at GENUS level
agp_infant_12m_genus

#DESEq (No need to rarefaction [use infant_6m_final])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
infant_12m_plus1 <- transform_sample_counts(agp_infant_12m_genus, function(x) x+1)
infant_12m_deseq <- phyloseq_to_deseq2(infant_12m_plus1, ~agp_clin)
DESEQ_infant_12m <- DESeq(infant_12m_deseq)

#Viewing Results
res_12m <- results(DESEQ_infant_12m, tidy=TRUE)
View(res_12m)

#Filtering out the ASV results with N.A. values for p adjusted values
filtered_res_12m <- drop_na(res_12m, padj)
view(filtered_res_12m)

#### Visualizing the results of DESEq ####
#Chose p value <0.05, and significant fold change >2

##Volcano Plot##
vol_plot_12m <- ggplot(filtered_res_12m, aes(x=log2FoldChange, y=-log(padj))) +
  geom_point()
vol_plot_sig_12m <- filtered_res_12m |>
  mutate(significant = padj<0.05 & abs(log2FoldChange) >2) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), col= significant)) +
  geom_point()

vol_plot_12m            
vol_plot_sig_12m

#Significant ASVs table
sigASVs_12m <- filtered_res_12m |>
  filter(padj <0.05 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

save(sigASVs_12m, file= "Files_from_Processing/DESeq files/DESEq_sigASVs_12m")

view(sigASVs_12m)

#BarPlot
#removed any Genus that has NA

sigASVs_vec_12m <- sigASVs_12m |>
  pull(ASV)

infant_12m_DESeq <- prune_taxa(sigASVs_vec_12m, agp_infant_12m_genus)

sigASVs_12m <- tax_table(infant_12m_DESeq) |>
  as.data.frame() |>
  rownames_to_column(var="ASV") |>
  right_join(sigASVs_12m) |>
  arrange(log2FoldChange) |>
  mutate(Genus = make.unique(Genus)) |>
  mutate(Genus = factor(Genus, levels= unique(Genus)))|>
  drop_na(Genus) |>
  filter(Genus != "NA.2" & Genus != "NA.1")

view(sigASVs_12m)

bar_plot_12m <- ggplot(sigASVs_12m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat ="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
bar_plot_12m

#Merged BarPlot
sigASVs_12m_mod <- sigASVs_12m |> add_column(age = "12")
sigASVs_6m_mod <- sigASVs_6m |> add_column(age = "6")
sigASVs_merged <- rbind(sigASVs_12m_mod, sigASVs_6m_mod)
view(sigASVs_merged)
bar_plot_merged <- ggplot(sigASVs_merged) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill=age), stat ="identity", position = "dodge") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

bar_plot_merged
################## INDICATOR SPECIES ANALYSIS ####################
#ISA 6 months
agp_genus_6m <- tax_glom(agp_infant_6m, "Genus", NArm= FALSE)
agp_genus_RA_6m <- transform_sample_counts(agp_genus_6m, fun=function(x) x/sum(x))

isa_agp_6m <- multipatt(t(otu_table(agp_genus_RA_6m)), cluster = sample_data(agp_genus_RA_6m)$agp_clin)
summary(isa_agp_6m)
taxtable_6m_agp <- tax_table(agp_infant_6m) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_agp_6m$sign %>% rownames_to_column(var="ASV") %>%
  left_join(taxtable_6m_agp) %>% filter(p.value<0.05) %>% View()

#ISA 12 months
agp_genus_12m <- tax_glom(agp_infant_12m, "Genus", NArm= FALSE)
agp_genus_RA_12m <- transform_sample_counts(agp_genus_12m, fun=function(x) x/sum(x))

isa_agp_12m <- multipatt(t(otu_table(agp_genus_RA_12m)), cluster = sample_data(agp_genus_RA_12m)$agp_clin)
summary(isa_agp_12m)
taxtable_12m_agp <- tax_table(agp_infant_12m) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_agp_12m$sign %>% rownames_to_column(var="ASV") %>%
  left_join(taxtable_12m_agp) %>% filter(p.value<0.05) %>% View()

