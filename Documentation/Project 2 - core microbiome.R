################## CORE MICROBIOME ANALYSIS ####################
#load libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggeasy)

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
agp_6m_above <- subset_samples(agp_abun_6m, agp_clin=="Above")
agp_6m_below <- subset_samples(agp_abun_6m, agp_clin=="Below")
crp_6m_above <- subset_samples(crp_abun_6m, crp_median=="Above")
crp_6m_below <- subset_samples(crp_abun_6m, crp_median=="Below")

#12 months
agp_12m_above <- subset_samples(agp_abun_12m, agp_clin=="Above")
agp_12m_below <- subset_samples(agp_abun_12m, agp_clin=="Below")
crp_12m_above <- subset_samples(crp_abun_12m, crp_median=="Above")
crp_12m_below <- subset_samples(crp_abun_12m, crp_median=="Below")

#Core microbiome for Venn Diagram
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
four_way_6m <- ggVennDiagram(x=list(core_agp_6m_above, core_agp_6m_below, core_crp_6m_above, core_crp_6m_below), category.names = c("High AGP","Low AGP","High CRP", "Low CRP"), label = "count")
four_way_6m + scale_x_continuous(expand = expansion(mult = .2)) + labs(title = "Core Taxa in 6-month-old Infants", fill = "# Genera") + theme(plot.title = element_text(size = 17, vjust = 8), plot.margin = margin(r = 20)) + easy_center_title()

four_way_12m <- ggVennDiagram(x=list(core_agp_12m_above, core_agp_12m_below, core_crp_12m_above, core_crp_12m_below), category.names = c("High AGP","Low AGP","High CRP", "Low CRP"), label = "count")
four_way_12m + scale_x_continuous(expand = expansion(mult = .2)) + labs(title = "Core Taxa in 12-month-old Infants", fill = "# Genera") + theme(plot.title = element_text(size = 17, vjust = 8), plot.margin = margin(r = 20)) + easy_center_title()

#Core phyloseq
#6 months
core_highagp_6m <- core(agp_6m_above, detection=0, prevalence=0.8)
core_lowagp_6m <- core(agp_6m_below, detection=0, prevalence=0.8)
core_highcrp_6m <- core(crp_6m_above, detection=0, prevalence=0.8)
core_lowcrp_6m <- core(crp_6m_below, detection=0, prevalence=0.8)

#12 months
core_highagp_12m <- core(agp_12m_above, detection=0, prevalence=0.8)
core_lowagp_12m <- core(agp_12m_below, detection=0, prevalence=0.8)
core_highcrp_12m <- core(crp_12m_above, detection=0, prevalence=0.8)
core_lowcrp_12m <- core(crp_12m_below, detection=0, prevalence=0.8)

#Identifying Core Taxa
#6 months
taxa_highagp_6m <- as.data.frame(tax_table(core_highagp_6m))
taxa_lowagp_6m <- as.data.frame(tax_table(core_lowagp_6m))
taxa_highcrp_6m <- as.data.frame(tax_table(core_highcrp_6m))
taxa_lowcrp_6m <- as.data.frame(tax_table(core_lowcrp_6m))

#12 months
taxa_highagp_12m <- as.data.frame(tax_table(core_highagp_12m))
taxa_lowagp_12m <- as.data.frame(tax_table(core_lowagp_12m))
taxa_highcrp_12m <- as.data.frame(tax_table(core_highcrp_12m))
taxa_lowcrp_12m <- as.data.frame(tax_table(core_lowcrp_12m))
