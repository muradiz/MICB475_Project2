### REQUIRED FILES:
# anemia-metadata.txt 
# tree.nwk
# taxonomy.tsv
# table_250.tsv.txt

library(tidyverse)
library(phyloseq)
library(dplyr)
library(ape)
library(vegan)
library(ggplot2)
library(ggforce)
library(DESeq2)
library(indicspecies)

anemia_metadata <- read_delim("Files_from_Processing/anemia_metadata.txt", delim="\t")
tax <- read_delim("Files_from_Processing/Exported_files/taxonomy.tsv", delim="\t")

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
save(meta_6m_filt, file = "sorted_metadata_crp_6m.RData")
save(meta_12m_filt, file = "sorted_metadata_crp_12m.RData")

# Save as txt files - NOT WORKING
write.table(meta_12m_filt, file = "sorted_metadata_crp_12m.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(meta_6m_filt, file = "sorted_metadata_crp_6m.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Group by and summarize
group_by(meta_6m_filt, crp_median)
group_by(meta_12m_filt, crp_median)

#load in metadata and taxonomy tsvs
otu <- read_delim(file="Files_from_Processing/Exported_files/table_250.tsv", delim = "\t", skip=1)
tax <- read_delim(file = "Files_from_Processing/Exported_files/taxonomy.tsv", delim="\t")
meta_6m <- read_delim("sorted_metadata_crp_6m.txt", delim="\t")
meta_12m <- read_delim("sorted_metadata_crp_12m.txt", delim="\t")
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

# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`

# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df_6m<- as.data.frame(meta_6m[,-1])
samp_df_12m <- as.data.frame(meta_12m[,-1])
# Make sampleids the rownames
rownames(samp_df_6m)<- meta_6m$sampleid
rownames(samp_df_12m)<- meta_12m$sampleid

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
crp_infant_6m <- phyloseq(OTU, SAMP_6m, TAX, phylotree)
crp_infant_12m <- phyloseq(OTU, SAMP_12m, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(crp_infant_6m)
sample_data(crp_infant_6m)
tax_table(crp_infant_6m)
phy_tree(crp_infant_6m)

otu_table(crp_infant_12m)
sample_data(crp_infant_12m)
tax_table(crp_infant_12m)
phy_tree(crp_infant_12m)

#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(crp_infant_6m)
sample_variables(crp_infant_12m)
# colnames(sample_data(crp_infant_12m))
get_variable(crp_infant_6m, c("crp", "crp_status")) # equivalent to "select" in tidyverse
get_variable(crp_infant_12m, c("crp", "crp_status"))

## Let's say we want to filter OTU table by sample. 
# What is the sum of reads in each sample?
sample_sums(crp_infant_6m)
sample_sums(crp_infant_12m)

# Save the sample names of the 3 samples with the most reads
getsamps_6m <- names(sort(sample_sums(crp_infant_6m), decreasing = TRUE)[1:3])
getsamps_12m <- names(sort(sample_sums(crp_infant_12m), decreasing = TRUE)[1:3])

# filter to see taxa abundances for each sample
get_taxa(crp_infant_6m, getsamps_6m) 
get_taxa(crp_infant_12m, getsamps_12m)

## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(crp_infant_6m)
taxa_names(crp_infant_12m)

# How many taxa do we have?
ntaxa(crp_infant_6m)
ntaxa(crp_infant_12m)

# What is the total read count for each taxa?
taxa_sums(crp_infant_6m)
taxa_sums(crp_infant_12m)

# Let's find the top 3 most abundant taxa
gettaxa_6m <- names(sort(taxa_sums(crp_infant_6m), decreasing = TRUE)[1:3] )
get_sample(crp_infant_6m, gettaxa_6m)

gettaxa_12m <- names(sort(taxa_sums(crp_infant_12m), decreasing = TRUE)[1:3] )
get_sample(crp_infant_12m, gettaxa_12m)


######### ANALYZE ##########
# Remove non-bacterial sequences, if any
crp_infant_6m_filt <- subset_taxa(crp_infant_6m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
crp_infant_12m_filt <- subset_taxa(crp_infant_12m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Remove ASVs that have less than 5 counts total
crp_infant_6m_filt_nolow <- filter_taxa(crp_infant_6m_filt, function(x) sum(x)>5, prune = TRUE)
crp_infant_12m_filt_nolow <- filter_taxa(crp_infant_12m_filt, function(x) sum(x)>5, prune = TRUE)


# Remove samples with less than 100 reads
crp_infant_6m_filt_nolow_samps <- prune_samples(sample_sums(crp_infant_6m_filt_nolow)>100, crp_infant_6m_filt_nolow)
crp_infant_12m_filt_nolow_samps <- prune_samples(sample_sums(crp_infant_12m_filt_nolow)>100, crp_infant_12m_filt_nolow)

# Remove samples where ph is na
crp_infant_6m_final <- subset_samples(crp_infant_6m_filt_nolow_samps, !is.na(crp))
crp_infant_12m_final <- subset_samples(crp_infant_12m_filt_nolow_samps, !is.na(crp))

# Rarefy samples #########
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
rarecurve(t(as.data.frame(otu_table(crp_infant_6m_final))), cex=0.1)
crp_infant_6m_rare <- rarefy_even_depth(crp_infant_6m_final, rngseed = 1, sample.size = 15000)

rarecurve(t(as.data.frame(otu_table(crp_infant_12m_final))), cex=0.1)
crp_infant_12m_rare <- rarefy_even_depth(crp_infant_12m_final, rngseed = 1, sample.size = 15000)

#view
crp_infant_6m_rare
crp_infant_12m_rare

##### Saving #####
save(crp_infant_6m_final, file="crp_infant_6m_final.RData")
save(crp_infant_6m_rare, file="crp_infant_6m_rare.RData")

save(crp_infant_12m_final, file="crp_infant_12m_final.RData")
save(crp_infant_12m_rare, file="crp_infant_12m_rare.RData")


###### DIVERSITY ###########
#6 months
##Alpha Diversity Plots
plot_richness(crp_infant_6m_rare)
# select certain alpha diversity metrics

plot_richness(crp_infant_6m_rare, measures = c("Shannon","Chao1"))
# Add ggplot layers if desired to adjust visuals
"6M_CRP_alpha" <- plot_richness(crp_infant_6m_rare, x = "crp_median", measures = c("Shannon")) +
  xlab("CRP Levels") +
  geom_boxplot()+geom_signif(comparisons = list(c("High","Low")),annotation = c("p=0.6915")) + ylim(0,5) + theme_bw()

ggsave("6M_CRP_alpha.png", height = 4, width = 3)

##Beta Diversity Plots
bc_dm_6m <- phyloseq::distance(crp_infant_6m_rare, method="bray")
pcoa_bc_6m <- ordinate(crp_infant_6m_rare, method="PCoA", distance=bc_dm_6m)

sixMCRPbeta <- plot_ordination(crp_infant_6m_rare, pcoa_bc_6m, color = "crp_median") +
  labs( col = "CRP Levels") + stat_ellipse()
ggsave("sixMCRPbeta.png", width = 5, height  = 3.5)


#create a taxa summaries plot
crp_infant_6m_RA <- transform_sample_counts(crp_infant_6m_rare, function(x) x/sum(x))
plot_bar(crp_infant_6m_RA, fill="Phylum") 
crp_infant_6m_phylum <- tax_glom(crp_infant_6m_RA, taxrank = "Phylum", NArm=FALSE)
plot_bar(crp_infant_6m_phylum, fill="Phylum") +
  facet_wrap(.~crp_median, scales = "free_x")

#12 months
##Alpha Diversity Plots
plot_richness(crp_infant_12m_rare)
# select certain alpha diversity metrics

plot_richness(crp_infant_12m_rare, measures = c("Shannon","Chao1"))
# Add ggplot layers if desired to adjust visuals
"12M_CRP_alpha" <- plot_richness(crp_infant_12m_rare, x = "crp_median", measures = c("Shannon")) +
  xlab("CRP Levels") +
  geom_boxplot() + geom_signif(comparisons = list(c("High","Low")),annotation = c("p=0.1781")) + ylim(0,5) + theme_bw()

ggsave("12M_CRP_alpha.png", height = 4, width = 3)

##Beta Diversity Plots
bc_dm_12m <- phyloseq::distance(crp_infant_12m_rare, method="bray")
pcoa_bc_12m <- ordinate(crp_infant_12m_rare, method="PCoA", distance=bc_dm_12m)

twelveMbeta <- plot_ordination(crp_infant_12m_rare, pcoa_bc_12m, color = "crp_median") +
  labs( col = "CRP Levels") + stat_ellipse()

ggsave("twelveMCRPbeta.png", width = 5, height  = 3.5)

#create a taxa summaries plot
crp_infant_12m_RA <- transform_sample_counts(crp_infant_12m_rare, function(x) x/sum(x))
plot_bar(crp_infant_12m_RA, fill="Phylum") 
crp_infant_12m_phylum <- tax_glom(crp_infant_12m_RA, taxrank = "Phylum", NArm=FALSE)
plot_bar(crp_infant_12m_phylum, fill="Phylum") +
  facet_wrap(.~crp_median, scales = "free_x")


########## STATISTICAL ANALYSIS ##############
#6 months
#load files
load("crp_infant_6m_rare.RData")

#t-test analysis
######## Comparison of two means with t-test (Parametric) ##########
# Let's do very simple plot with t-test
plot_richness(crp_infant_6m_rare, x = "crp_median", measures="Shannon")
# Need to extract information
alphadiv_6m <- estimate_richness(crp_infant_6m_rare)
samp_dat_6m <- sample_data(crp_infant_6m_rare)
samp_dat_wdiv_6m <- data.frame(samp_dat_6m, alphadiv_6m)
# These are equivalent:
# t.test()
t.test(samp_dat_wdiv_6m$Shannon ~ samp_dat_wdiv_6m$crp_median)
t.test(Shannon ~ crp_median, data=samp_dat_wdiv_6m)

# Note: you can set variances to be equal for a "classic" t-test
t.test(Shannon ~ crp_median, data=samp_dat_wdiv_6m, var.equal=TRUE)


#### Microbial count data is generally NON-NORMAL ####
# In fact, it is even more complex because microbial data is usually in RELATIVE ABUNDANCE
allCounts_6m <- as.vector(otu_table(crp_infant_6m_rare))
allCounts_6m <- allCounts_6m[allCounts_6m>0]
hist(allCounts_6m)
hist(log(allCounts_6m))


#remove 0 values
samp_dat_wdiv_6m %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=crp_median, y=Shannon))+
  geom_point() 

#check distribution:
ggplot(samp_dat_wdiv_6m) +
  geom_histogram(aes(x=Shannon), bins=25)

#logging everything to deal with non-parametric distribution
# (1) Transform your data (usually with a log function)
ggplot(samp_dat_wdiv_6m) +
  geom_histogram(aes(x=log(Shannon)), bins=25)
t.test(log(Shannon) ~ crp_median, data=samp_dat_wdiv_6m)
# Let's see what transformed data looks like:
samp_dat_wdiv_6m %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=crp_median, y=log(Shannon)))+
  geom_boxplot() +
  geom_jitter()


#Wilcoxon Rank Sum Test
wilcox.test(Shannon ~ crp_median, data=samp_dat_wdiv_6m)
wilcox.test(log(Shannon) ~ crp_median, data=samp_dat_wdiv_6m)

#PERMANOVA Analysis
##load data
samp_dat_wdiv_6m <- data.frame(sample_data(crp_infant_6m_rare), estimate_richness(crp_infant_6m_rare))

### PERMANOVA (Permutational ANOVA) ####
# non-parametric version of ANOVA
# Takes a distance matrix, which can be calculated with any kind of metric you want
# e.g. Bray, Jaccard, Unifrac
# Need the package, "vegan"
# Use phyloseq to calculate weighted Unifrac distance matrix
?UniFrac
dm_unifrac_6m <- UniFrac(crp_infant_6m_rare, weighted=TRUE)
?adonis2
adonis2(dm_unifrac_6m ~ crp_median, data=samp_dat_wdiv_6m)

# Also use other metrics: for example, the vegan package includes bray and jaccard
dm_bray_6m <- vegdist(t(otu_table(crp_infant_6m_rare)), method="bray")
adonis2(dm_bray_6m ~ crp_median, data=samp_dat_wdiv_6m)

dm_jaccard_6m <- vegdist(t(otu_table(crp_infant_6m_rare)), method="jaccard")
adonis2(dm_jaccard_6m ~ crp_median, data=samp_dat_wdiv_6m)


#12 months
#load files
load("crp_infant_12m_rare.RData")

#t-test analysis
######## Comparison of two means with t-test (Parametric) ##########
# Let's do very simple plot with t-test
plot_richness(crp_infant_12m_rare, x = "crp_median", measures="Shannon")
# Need to extract information
alphadiv_12m <- estimate_richness(crp_infant_12m_rare)
samp_dat_12m <- sample_data(crp_infant_12m_rare)
samp_dat_wdiv_12m <- data.frame(samp_dat_12m, alphadiv_12m)
# These are equivalent:
# t.test()
t.test(samp_dat_wdiv_12m$Shannon ~ samp_dat_wdiv_12m$crp_median)
t.test(Shannon ~ crp_median, data=samp_dat_wdiv_12m)

# Note: you can set variances to be equal for a "classic" t-test
t.test(Shannon ~ crp_median, data=samp_dat_wdiv_12m, var.equal=TRUE)


#### Microbial count data is generally NON-NORMAL ####
# In fact, it is even more complex because microbial data is usually in RELATIVE ABUNDANCE
allCounts_12m <- as.vector(otu_table(crp_infant_12m_rare))
allCounts_12m <- allCounts_12m[allCounts_12m>0]
hist(allCounts_12m)
hist(log(allCounts_12m))


#remove 0 values
samp_dat_wdiv_12m %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=crp_median, y=Shannon))+
  geom_point() 

#check distribution:
ggplot(samp_dat_wdiv_12m) +
  geom_histogram(aes(x=Shannon), bins=25)

#logging everything to deal with non-parametric distribution
# (1) Transform your data (usually with a log function)
ggplot(samp_dat_wdiv_12m) +
  geom_histogram(aes(x=log(Shannon)), bins=25)
t.test(log(Shannon) ~ crp_median, data=samp_dat_wdiv_12m)
# Let's see what transformed data looks like:
samp_dat_wdiv_12m %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=crp_median, y=log(Shannon)))+
  geom_boxplot() +
  geom_jitter()


#Wilcoxon Rank Sum Test
wilcox.test(Shannon ~ crp_median, data=samp_dat_wdiv_12m)
wilcox.test(log(Shannon) ~ crp_median, data=samp_dat_wdiv_12m)

#PERMANOVA Analysis
##load data
samp_dat_wdiv_12m <- data.frame(sample_data(crp_infant_12m_rare), estimate_richness(crp_infant_12m_rare))

### PERMANOVA (Permutational ANOVA) ####
# non-parametric version of ANOVA
# Takes a distance matrix, which can be calculated with any kind of metric you want
# e.g. Bray, Jaccard, Unifrac
# Need the package, "vegan"
# Use phyloseq to calculate weighted Unifrac distance matrix
?UniFrac
dm_unifrac_12m <- UniFrac(crp_infant_12m_rare, weighted=TRUE)
?adonis2
adonis2(dm_unifrac_12m ~ crp_median, data=samp_dat_wdiv_12m)

# Also use other metrics: for example, the vegan package includes bray and jaccard
dm_bray_12m <- vegdist(t(otu_table(crp_infant_12m_rare)), method="bray")
adonis2(dm_bray_12m ~ crp_median, data=samp_dat_wdiv_12m)

dm_jaccard_12m <- vegdist(t(otu_table(crp_infant_12m_rare)), method="jaccard")
adonis2(dm_jaccard_12m ~ crp_median, data=samp_dat_wdiv_12m)


############ DESEQ ANALYSIS ##################
#6 months
#setting random seed 
set.seed(1)

#Filtering Genus-level, no need to rarefy

# Genus Level
#crp_infant_6m_genus <- tax_glom(crp_infant_6m_filt_nolow, "Genus")
#crp_infant_12m_genus <- tax_glom(crp_infant_12m_filt_nolow, "Genus")

# Remove samples with less than 100 reads
crp_infant_6m_filt_nolow_samps_DESeq <- prune_samples(sample_sums(crp_infant_6m_final)>100, crp_infant_6m_final)
crp_infant_12m_filt_nolow_samps_DESeq <- prune_samples(sample_sums(crp_infant_12m_final)>100, crp_infant_12m_final)

# Remove samples where crp is na
crp_infant_6m_final_DESeq <- subset_samples(crp_infant_6m_filt_nolow_samps_DESeq, !is.na(crp))
crp_infant_12m_final_DESeq <- subset_samples(crp_infant_12m_filt_nolow_samps_DESeq, !is.na(crp))


#DESEq (No need to rarefaction [use crp_infant_6m_final_DESeq])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
crp_infant_6m_plus1 <- transform_sample_counts(crp_infant_6m_final_DESeq, function(x) x+1)
crp_infant_6m_deseq <- phyloseq_to_deseq2(crp_infant_6m_plus1, ~crp_median)
crp_DESEQ_infant_6m <- DESeq(crp_infant_6m_deseq)

#Viewing Results
res_6m <- results(crp_DESEQ_infant_6m, tidy=TRUE, contrast= c("crp_median", "Above", "Below"))
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
crp_sigASVs_6m <- filtered_res_6m |>
  filter(padj <0.05 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

save(crp_sigASVs_6m, file= "Files_from_Processing/DESeq files/DESEq_sigASVs_crp_6m_final")

view(crp_sigASVs_6m)

#BarPlot
#removed any Genus that has NA

sigASVs_vec_6m <- crp_sigASVs_6m |>
  pull(ASV)

crp_infant_6m_DESeq <- prune_taxa(sigASVs_vec_6m, crp_infant_6m_final_DESeq)

crp_sigASVs_6m <- tax_table(crp_infant_6m_DESeq) |>
  as.data.frame() |>
  rownames_to_column(var="ASV") |>
  right_join(crp_sigASVs_6m) |>
  arrange(log2FoldChange) |>
  mutate(Genus = make.unique(Genus)) |>
  mutate(Genus = factor(Genus, levels= unique(Genus)))|>
  drop_na(Genus) |>
  filter(Genus != "NA.2" & Genus != "NA.1")

view(crp_sigASVs_6m)

bar_plot_6m <- ggplot(crp_sigASVs_6m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat ="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
bar_plot_6m

#12 months
#setting random seed 
set.seed(1)

#DESEq (No need to rarefaction [use crp_infant_12m_final_DESeq])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
crp_infant_12m_plus1 <- transform_sample_counts(crp_infant_12m_final_DESeq, function(x) x+1)
crp_infant_12m_deseq <- phyloseq_to_deseq2(crp_infant_12m_plus1, ~crp_median)
crp_DESEQ_infant_12m <- DESeq(crp_infant_12m_deseq)

#Viewing Results
res_12m <- results(crp_DESEQ_infant_12m, tidy=TRUE, contrast= c("crp_median", "Above", "Below"))
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
crp_sigASVs_12m <- filtered_res_12m |>
  filter(padj <0.05 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

save(crp_sigASVs_12m, file= "Files_from_Processing/DESeq files/DESEq_sigASVs_crp_12m_final")

view(crp_sigASVs_12m)

#BarPlot
#removed any Genus that has NA

sigASVs_vec_12m <- crp_sigASVs_12m |>
  pull(ASV)

crp_infant_12m_DESeq <- prune_taxa(sigASVs_vec_12m, crp_infant_12m_final_DESeq)

crp_sigASVs_12m <- tax_table(crp_infant_12m_DESeq) |>
  as.data.frame() |>
  rownames_to_column(var="ASV") |>
  right_join(crp_sigASVs_12m) |>
  arrange(log2FoldChange) |>
  mutate(Genus = make.unique(Genus)) |>
  mutate(Genus = factor(Genus, levels= unique(Genus)))|>
  drop_na(Genus) |>
  filter(Genus != "NA.2" & Genus != "NA.1")

view(crp_sigASVs_12m)

bar_plot_12m <- ggplot(crp_sigASVs_12m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat ="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
bar_plot_12m

#Merged BarPlot
#crp_sigASVs_12m_mod <- crp_sigASVs_12m |> add_column(age = "12")
#crp_sigASVs_6m_mod <- crp_sigASVs_6m |> add_column(age = "6")
#crp_sigASVs_merged <- rbind(crp_sigASVs_12m_mod, crp_sigASVs_6m_mod)
#view(crp_sigASVs_merged)
#bar_plot_merged <- ggplot(crp_sigASVs_merged) +
#  geom_bar(aes(x=Genus, y=log2FoldChange, fill=age), stat ="identity", position = "dodge") +
#  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
#  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#bar_plot_merged

#Merged BarPlot (New Code)
sigASVs_merged <- crp_sigASVs_12m |>
  rbind(crp_sigASVs_6m) |>
  pull(ASV) |>
  unique()

all_lfc <- res_12m |>
  mutate(age ='12') |>
  rbind(res_6m |> mutate(age='6')) |>
  filter(row %in% sigASVs_merged)

names(all_lfc)[1] = 'ASV'
#view(all_lfc)              

tax_tbl_crp_12m <- crp_infant_12m_final@tax_table |>
  as.matrix() |>
  as.data.frame() 

tax_tbl_crp_6m <- crp_infant_6m_final@tax_table |>
  as.matrix() |>
  as.data.frame()

tax_tbl_merged <- rbind(tax_tbl_crp_6m, tax_tbl_crp_12m) |>
  unique()


all_lfc = all_lfc %>% left_join(tax_table(crp_infant_12m_final) %>% as.data.frame() %>% rownames_to_column('ASV') %>% select(ASV,Genus))

all_lfc = all_lfc %>% mutate(Significant = ifelse(padj<0.05,T,F)) |>
  mutate(fold_change = ifelse(log2FoldChange<0, "Increased in Low CRP", "Increased in High CRP")) |>
  filter(Genus != "NA") |>
  filter(Significant == TRUE) |>
  unique()

#order = all_lfc %>% filter(age=="12") %>% arrange(log2FoldChange) %>% mutate(Genus = factor(Genus, levels = unique(.$Genus)))
#all_lfc = all_lfc %>% mutate(Genus = factor(Genus, levels = unique(order$Genus))) |> filter(Genus != "NA") |>

crp_otu_key <- all_lfc |> select(ASV,Genus) %>% unique() %>% group_by(Genus) %>%
  mutate(Genus2 = paste(Genus,row_number(),sep='.')) %>% ungroup()

all_lfc <- all_lfc |> left_join(crp_otu_key)
view(all_lfc)

view(all_lfc)

bar_plot_mergedv2 <- all_lfc %>% 
  ggplot(aes(reorder(Genus2, -log2FoldChange),log2FoldChange, fill=fold_change)) +
  scale_alpha_manual(values = c(0.3,2)) +
  geom_col() +
  geom_hline(yintercept=0, linewidth = 0.75)+
  geom_errorbar(aes(x=Genus2, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        strip.text.y.right = element_text(angle=0), text = element_text(size = 16)) +
  labs(fill = "Fold Change", x = "\n ASVs mapped to Genus", y ="Log2FoldChange (High/Low CRP)") +
  facet_grid(rows = vars(age))

bar_plot_mergedv2
################ INDICATOR SPECIES ANALYSIS ###############
#6 months
#Setting random seed as 1 
set.seed(1)

crp_genus_6m <- tax_glom(crp_infant_6m, "Genus", NArm= FALSE)
crp_genus_RA_6m <- transform_sample_counts(crp_genus_6m, fun=function(x) x/sum(x))

isa_crp_6m <- multipatt(t(otu_table(crp_genus_RA_6m)), cluster = sample_data(crp_genus_RA_6m)$crp_median)
summary(isa_crp_6m)
taxtable_6m_crp <- tax_table(crp_infant_6m) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_crp_6m$sign %>% rownames_to_column(var="ASV") %>%
  left_join(taxtable_6m_crp) %>% filter(p.value<0.05) %>% View()

#12 months
#Setting random seed as 1 
set.seed(1)
crp_genus_12m <- tax_glom(crp_infant_12m, "Genus", NArm= FALSE)
crp_genus_RA_12m <- transform_sample_counts(crp_genus_12m, fun=function(x) x/sum(x))

#ISA CRP 12 months
isa_crp_12m <- multipatt(t(otu_table(crp_genus_RA_12m)), cluster = sample_data(crp_genus_RA_12m)$crp_median)
summary(isa_crp_12m)
taxtable_12m_crp <- tax_table(crp_infant_12m) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_crp_12m$sign %>% rownames_to_column(var="ASV") %>%
  left_join(taxtable_12m_crp) %>% filter(p.value<0.05) %>% View()


