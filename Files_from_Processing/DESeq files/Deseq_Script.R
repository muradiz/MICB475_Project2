#!/usr/bin/env Rscript

library(tidyverse)
library(phyloseq)
library(DESeq2)

#Loading the data (FIX THIS CODE)
load("infant_12m_final")

#Need to filter anything????

#DESEq
infant_12m_deseq <- phyloseq_to_deseq2()