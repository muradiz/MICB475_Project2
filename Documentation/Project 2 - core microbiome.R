################## CORE MICROBIOME ANALYSIS ####################
#load libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggeasy)
library(RColorBrewer)
library(ggpubr)

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
core_agp_6m_above <- core_members(agp_6m_above, detection=0, prevalence=0.7)
core_agp_6m_below <- core_members(agp_6m_below, detection=0, prevalence=0.7)
core_crp_6m_above <- core_members(crp_6m_above, detection=0, prevalence=0.7)
core_crp_6m_below <- core_members(crp_6m_below, detection=0, prevalence=0.7)
  
#12 months
core_agp_12m_above <- core_members(agp_12m_above, detection=0, prevalence=0.7)
core_agp_12m_below <- core_members(agp_12m_below, detection=0, prevalence=0.7)
core_crp_12m_above <- core_members(crp_12m_above, detection=0, prevalence=0.7)
core_crp_12m_below <- core_members(crp_12m_below, detection=0, prevalence=0.7)

#Visualization

#4-way 6 months 
four_way_6m <- ggVennDiagram(x=list(core_agp_6m_above, core_agp_6m_below, core_crp_6m_above, core_crp_6m_below), category.names = c("High AGP","Low AGP","High CRP", "Low CRP"), legend.position = "right", label = "count", label_size = 4)
core_6m <- four_way_6m_count + scale_x_continuous(expand = expansion(mult = .2)) + labs(caption = "6-month-old infants", fill = "# Genera") + theme(plot.caption = element_text(size = 12, hjust=0.5), plot.margin = margin(r = 20)) + scale_fill_gradient(low = "#EBE8FC", high = "#756BB1") + border() 
core_6m
ggsave("four_way_6m.png", plot = last_plot(), height = 3.5, width = 5, dpi = 300)

#4-way 12 month (count)
four_way_12m_count <- ggVennDiagram(x=list(core_agp_12m_above, core_agp_12m_below, core_crp_12m_above, core_crp_12m_below), category.names = c("High AGP","Low AGP","High CRP", "Low CRP"), legend.position = "right", label = "count", label_size = 4)
core_12m <- four_way_12m_count + scale_x_continuous(expand = expansion(mult = .2)) + labs(caption = "12-month-old infants", fill = "# Genera") + theme(plot.caption = element_text(size = 12, hjust=0.5), plot.margin = margin(r = 20)) + scale_fill_gradient(low = "#EBE8FC", high = "#756BB1") + border()
core_12m
ggsave("four_way_12m.png", plot = last_plot(), height = 3.5, width = 5, dpi = 300)

ggarrange(core_6m, core_12m, ncol = 2, nrow = 1) + bgcolor("white")

ggsave("combinedcore_plot.png", plot = last_plot(), height = 7, width = 15, dpi = 300)

#Core phyloseq
#6 months
core_highagp_6m <- core(agp_6m_above, detection=0, prevalence=0.7)
core_lowagp_6m <- core(agp_6m_below, detection=0, prevalence=0.7)
core_highcrp_6m <- core(crp_6m_above, detection=0, prevalence=0.7)
core_lowcrp_6m <- core(crp_6m_below, detection=0, prevalence=0.7)

#12 months
core_highagp_12m <- core(agp_12m_above, detection=0, prevalence=0.7)
core_lowagp_12m <- core(agp_12m_below, detection=0, prevalence=0.7)
core_highcrp_12m <- core(crp_12m_above, detection=0, prevalence=0.7)
core_lowcrp_12m <- core(crp_12m_below, detection=0, prevalence=0.7)

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
