#load packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(dplyr)


#load files
load("agp_infant_12m_rare.RData")

#t-test analysis
######## Comparison of two means with t-test (Parametric) ##########
# Let's do very simple plot with t-test
plot_richness(agp_infant_12m_rare, x = "agp_clin", measures="Shannon")
# Need to extract information
alphadiv <- estimate_richness(agp_infant_12m_rare)
samp_dat <- sample_data(agp_infant_12m_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)
# These are equivalent:
# t.test()
t.test(samp_dat_wdiv$Shannon ~ samp_dat_wdiv$agp_clin)
t.test(Shannon ~ agp_clin, data=samp_dat_wdiv)

# Note: you can set variances to be equal for a "classic" t-test
t.test(Shannon ~ agp_clin, data=samp_dat_wdiv, var.equal=TRUE)


#### Microbial count data is generally NON-NORMAL ####
# In fact, it is even more complex because microbial data is usually in RELATIVE ABUNDANCE
allCounts <- as.vector(otu_table(agp_infant_12m_rare))
allCounts <- allCounts[allCounts>0]
hist(allCounts)
hist(log(allCounts))


#remove 0 values
samp_dat_wdiv %>%
  filter(!is.na(agp)) %>%
  ggplot(aes(x=agp_clin, y=agp))+
  geom_point() 

#check distribution:
ggplot(samp_dat_wdiv) +
  geom_histogram(aes(x=agp), bins=25)

#logging everything to deal with non-parametric distribution
# (1) Transform your data (usually with a log function)
ggplot(samp_dat_wdiv) +
  geom_histogram(aes(x=log(agp)), bins=25)
t.test(log(agp) ~ agp_clin, data=samp_dat_wdiv)
# Let's see what transformed data looks like:
samp_dat_wdiv %>%
  filter(!is.na(agp)) %>%
  ggplot(aes(x=agp_clin, y=log(agp)))+
  geom_boxplot() +
  geom_jitter()


#Wilcoxon Rank Sum Test
wilcox.test(agp ~ agp_clin, data=samp_dat_wdiv)
wilcox.test(log(agp) ~ agp_clin, data=samp_dat_wdiv)



#PERMANOVA Analysis
##load data
samp_dat_wdiv <- data.frame(sample_data(agp_infant_12m_rare), estimate_richness(agp_infant_12m_rare))

### PERMANOVA (Permutational ANOVA) ####
# non-parametric version of ANOVA
# Takes a distance matrix, which can be calculated with any kind of metric you want
# e.g. Bray, Jaccard, Unifrac
# Need the package, "vegan"
# Use phyloseq to calculate weighted Unifrac distance matrix
?UniFrac
dm_unifrac <- UniFrac(agp_infant_12m_rare, weighted=TRUE)
?adonis2
adonis2(dm_unifrac ~ agp_clin*agp, data=samp_dat_wdiv)

# Also use other metrics: for example, the vegan package includes bray and jaccard
dm_bray <- vegdist(t(otu_table(agp_infant_12m_rare)), method="bray")
adonis2(dm_bray ~ agp_clin*agp, data=samp_dat_wdiv)

dm_jaccard <- vegdist(t(otu_table(agp_infant_12m_rare)), method="jaccard")
adonis2(dm_jaccard ~ agp_clin*agp, data=samp_dat_wdiv)
