
### REQUIRED FILES to generate phyloseq object
# anemia-metadata.txt 
# tree.nwk
# taxonomy.tsv
# table_250.tsv.txt

#Code from line 8-166 to make phyloseq object 
library(tidyverse)
library(phyloseq)
library(dplyr)
library(ape)
library(vegan)

anemia_metadata <- read_delim("Files_from_Processing/DESeq files/anemia_metadata.txt", delim="\t")
tax <- read_delim("Files_from_Processing/DESeq files/taxonomy.tsv", delim="\t")

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
otu <- read_delim(file="Files_from_Processing/DESeq files/table_250.tsv", delim = "\t", skip=1)
tax <- read_delim(file = "Files_from_Processing/DESeq files/taxonomy.tsv", delim="\t")
meta <- read_delim("sorted_metadata_12m.txt", delim="\t")
phylotree <- read.tree(file = "Files_from_Processing/DESeq files/tree.nwk")

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
infant_12m <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(infant_12m)
sample_data(infant_12m)
tax_table(infant_12m)
phy_tree(infant_12m)

library(vegan)
library(phyloseq)
library(tidyverse)

#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(infant_12m)
# colnames(sample_data(infant_12m))
get_variable(infant_12m, c("crp", "crp_status")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# What is the sum of reads in each sample?
sample_sums(infant_12m)
# Save the sample names of the 3 samples with the most reads
getsamps <- names(sort(sample_sums(infant_12m), decreasing = TRUE)[1:3])
# filter to see taxa abundances for each sample
get_taxa(infant_12m, getsamps) 


## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(infant_12m)
# How many taxa do we have?
ntaxa(infant_12m)
# What is the total read count for each taxa?
taxa_sums(infant_12m)
# Let's find the top 3 most abundant taxa
gettaxa <- names(sort(taxa_sums(infant_12m), decreasing = TRUE)[1:3] )
get_sample(infant_12m, gettaxa)


######### ANALYZE ##########
# Remove non-bacterial sequences, if any
infant_12m_filt <- subset_taxa(infant_12m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
infant_12m_filt_nolow <- filter_taxa(infant_12m_filt, function(x) sum(x)>5, prune = TRUE)
# Genus Level
infant_12m_genus <- tax_glom (infant_12m_filt_nolow, "Genus")
# Remove samples with less than 100 reads
infant_12m_filt_nolow_samps <- prune_samples(sample_sums(infant_12m_genus)>100, infant_12m_genus)
# Remove samples where ph is na
infant_12m_final <- subset_samples(infant_12m_filt_nolow_samps, !is.na(crp) )

####DESEq ANALYSIS#### 

#!/usr/bin/env Rscript

#Loading the libraries
library(tidyverse)
library(phyloseq)
library(DESeq2)

#setting random seed 
set.seed(1)

#Need to filter anything????
infant_12m_final

#DESEq (No need to rarefaction [use infant_12m_final])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
infant_12m_plus1 <- transform_sample_counts(infant_12m_final, function(x) x+1)
infant_12m_deseq <- phyloseq_to_deseq2(infant_12m_plus1, ~crp_median)
DESEQ_infant_12m <- DESeq(infant_12m_deseq)

#Viewing Results
res <- results(DESEQ_infant_12m, tidy=TRUE)
View(res)

#Filtering out the ASV results with N.A. values for p adjusted values
filtered_res <- drop_na(res, padj)
view(filtered_res)

#### Visualizing the results of DESEq ####
#Chose p value <0.05, and significant fold change >2


##Volcano Plot##
vol_plot <- ggplot(filtered_res, aes(x=log2FoldChange, y=-log(padj))) +
            geom_point()
vol_plot_sig <- filtered_res |>
            mutate(significant = padj<0.05 & abs(log2FoldChange) >2) |>
            ggplot(aes(x=log2FoldChange, y=-log10(padj), col= significant)) +
            geom_point()

vol_plot            
vol_plot_sig

#Significant ASVs table
sigASVs <- filtered_res |>
  filter(padj <0.05 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

save(sigASVs, file= "Files_from_Processing/DESeq files/DESEq_sigASVs")

view(sigASVs)

#BarPlot
#removed any Genus that has NA

sigASVs_vec <- sigASVs |>
                pull(ASV)

infant_12m_DESeq <- prune_taxa(sigASVs_vec, infant_12m_final)

sigASVs <- tax_table(infant_12m_DESeq) |>
            as.data.frame() |>
            rownames_to_column(var="ASV") |>
            right_join(sigASVs) |>
            arrange(log2FoldChange) |>
            mutate(Genus = make.unique(Genus)) |>
            mutate(Genus = factor(Genus, levels= unique(Genus)))|>
            drop_na(Genus) |>
            filter(Genus != "NA.2" & Genus != "NA.1") 

view(sigASVs)

bar_plot <- ggplot(sigASVs) +
            geom_bar(aes(x=Genus, y=log2FoldChange), stat ="identity") +
            geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
            theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
bar_plot
