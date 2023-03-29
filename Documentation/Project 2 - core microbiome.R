################## CORE MICROBIOME ANALYSIS ####################
#load libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#load files
agp_infant_6m_final <- load("agp_infant_6m_final.RData")
agp_infant_12m_final <- load("agp_infant_12m_final.RData")
crp_infant_6m_final <- load("crp_infant_6m_final.RData")
crp_infant_12m_final <- load("crp_infant_12m_final.RData")

#Converting to relative abundance
#6 months 
agp_abun_6m <- transform_sample_counts(agp_infant_6m_final, fun=function(x) x/sum(x))
crp_abun_6m <- transform_sample_counts(crp_infant_6m_final, fun=function(x) x/sum(x))

#12 months
agp_abun_12m <- transform_sample_counts(agp_infant_12m_final, fun=function(x) x/sum(x))
crp_abun_12m <- transform_sample_counts(crp_infant_12m_final, fun=function(x) x/sum(x))

#Subsetting 
#6 months
agp_6m_above <- subset_samples(agp_infant_6m_final, agp_clin=="Above")
agp_6m_below <- subset_samples(agp_infant_6m_final, agp_clin=="Below")
crp_6m_above <- subset_samples(crp_infant_6m_final, crp_median=="Above")
crp_6m_below <- subset_samples(crp_infant_6m_final, crp_median=="Below")

#12 months
agp_12m_above <- subset_samples(agp_infant_12m_final, agp_clin=="Above")
agp_12m_below <- subset_samples(agp_infant_12m_final, agp_clin=="Below")
crp_12m_above <- subset_samples(crp_infant_12m_final, crp_median=="Above")
crp_12m_below <- subset_samples(crp_infant_12m_final, crp_median=="Below")

#Core microbiome
#6 months
core_agp_6m_above <- core_members(agp_6m_above, detection=0, prevalence=0.8)
core_agp_6m_below <- core_members(agp_6m_below, detection=0, prevalence=0.8)
core_crp_6m_above <- core_members(crp_6m_above, detection=0, prevalence=0.8)
core_crp_6m_below <- core_members(crp_6m_below, detection=0, prevalence=0.8)
  
#12 months
core_agp_12m_above <- core_members(agp_12m_above, detection=0, prevalence=0.8)
core_agp_12m_below <- core_members(agp_12m_below, detection=0, prevalence=0.8)
core_crp_12m_above <- core_members(crp_12m_above, detection=0, prevalence=0.8)
core_crp_12m_below <- core_members(crp_12m_below, detection=0, prevalence=0.8)

#Visualization
#4-way 6 months
four_way_6m <- ggVennDiagram(x=list(core_agp_6m_above, core_agp_6m_below, core_crp_6m_above, core_crp_6m_below), category.names = c("AGP, Above","AGP, Below","CRP, Above", "CRP, Below"), label = "count")
four_way_6m + scale_x_continuous(expand = expansion(mult = .2)) + labs(title = "Four-way Venn Diagram for 6-month-old infants")

four_way_12m <- ggVennDiagram(x=list(core_agp_12m_above, core_agp_12m_below, core_crp_12m_above, core_crp_12m_below), category.names = c("AGP, Above","AGP, Below","CRP, Above", "CRP, Below"), label = "count")
four_way_12m + scale_x_continuous(expand = expansion(mult = .2)) + labs(title = "Four-way Venn Diagram for 12-month-old infants")
