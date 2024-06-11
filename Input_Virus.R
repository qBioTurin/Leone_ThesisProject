library(phyloseq)
library("ape")
library(ggplot2)
library(readxl)
library(dplyr)
library(magrittr)
library("ggtree")
library(decontam)
library(gridExtra)
library(patchwork)
library(microbiome)
library("ggVennDiagram")
library(compositions)
library('hilldiv')
library(tidyverse)
library(metagenomeSeq)
library(functional)
setwd("C:/Users/alele/iCloudDrive/Tirocinio/Virus")
#setwd("/Users/alessandro/Library/Mobile Documents/com~apple~CloudDocs/Tirocinio/Virus/")
biom = import_biom("BOV.biom")
Cohort = read_excel("CohortDescription.xlsx")
sample = read_excel("csf.xlsx", range = "A1:BZ85", col_names = TRUE, .name_repair = "universal", na = "")
biom@tax_table@.Data = substring(biom@tax_table@.Data, 4)
colnames(biom@tax_table@.Data)= c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
col = c("Sex","HIV_CSF","HIV_plasma","CSF_Escape","Barrier","Intrathec_synthesis","Regimen","NNRTI","PI","INSTI","Altered.Tau","Altered.ptau","Altered.BetAm","Altered.S100","Altered_Neopt","ETNIA","FdR","Pregr_Lue","Nadir.below.200","Nadir.below100","CSF_cell","TEST_NPS","Neurocogn.Impairment","Symptomatic.NCI","Depression.any","Depression.moderate.severe","Depression.severe","Detectable.CXCL13","Batch")
sample[,col] = lapply(sample[,col], as.factor) 
sample$Data.prelievo = as.POSIXct(sample$Data.prelievo, format = "%Y-%m-%dT%H:%M")
sample$Age=round(sample$Age)
rm(col)
names(sample) = gsub("\\.$", "%", names(sample)) 
sample_names(biom)[c(83, 54, 55)] = c('Neg1', 'Neg2', 'Neg3')
sample_names(biom) = sample_names(biom) %>% 
  gsub("sample", "", .) %>% 
  sub("_.*", "", .)
Cohort[, c("age", "CSF", "plasma")] = sapply(Cohort[, c("age", "CSF", "plasma")], as.numeric) %>% round()
sample$HIV.Diagnosis_da_mesi = round(sample$HIV.Diagnosis_da_mesi)
Cohort$Sample[c(82, 83, 84)] = c('Neg1', 'Neg2', 'Neg3')
sample$Paziente[c(82, 83, 84)] = c('Neg1', 'Neg2', 'Neg3')
sample= sample[-6,] # Rimuovo il campione 172 
Cohort = Cohort[-24,] # Same di riga 37 
sample = left_join(data.frame(Paziente=sample_names(biom)), sample, by="Paziente")
Cohort = Cohort[match(sample$Paziente, Cohort$Sample),]
sample$Quantification = Cohort$quantification
sample$SoC = Cohort$Sample_or_Control
samplephylo = phyloseq(sample_data(sample))
sample_names(samplephylo) = sample_names(biom)
samplephylo = samplephylo[!(samplephylo$Paziente %in% c(406, 178, 174, 184, 420, 464, 386)),]
sample = sample[!(sample$Paziente %in% c(406, 178, 174, 184, 420, 464, 386)),]
Neg = c("Neg1", "Neg2", "Neg3")
sampleNN = prune_samples(!(sample_names(samplephylo)%in% Neg), samplephylo)
OTU = as.data.frame(biom@otu_table@.Data)
OTU = OTU[ , ! colnames(OTU) %in% c('406', '178', '174', '184', '420', '464', '386')]
OTUnn = OTU[ , ! colnames(OTU) %in% c('406', '178', '174', '184', '420', '464', '386', 'Neg1', 'Neg2', 'Neg3')]
samplephylo$Batch[c(76,52,53)] = c(1,2,3)
s1nn = sampleNN[sampleNN$Batch==1 ]
s2nn = sampleNN[sampleNN$Batch==2 ]
s3nn = sampleNN[sampleNN$Batch==3 ]
s1 = samplephylo[samplephylo$Batch==1 ]
s2 = samplephylo[samplephylo$Batch==2 ]
s3 = samplephylo[samplephylo$Batch==3 ]

tax = as.data.frame(biom@tax_table@.Data) 
OTU = OTU[rownames(OTU) %in% rownames(tax), ]
OTUnn = OTUnn[rownames(OTUnn) %in% rownames(tax), ]
OTU1 = OTU[,s1$Paziente] 
OTU2 = OTU[,s2$Paziente]
OTU3 = OTU[,s3$Paziente]
OTU1nn = OTUnn[,s1nn$Paziente]
OTU2nn = OTUnn[,s2nn$Paziente]
OTU3nn = OTUnn[,s3nn$Paziente]
Phyloseq1 = phyloseq(otu_table(OTU1, taxa_are_rows = TRUE), tax_table(biom), sample_data(s1))
Phyloseq2 = phyloseq(otu_table(OTU2, taxa_are_rows = TRUE), tax_table(biom), sample_data(s2))
Phyloseq3 = phyloseq(otu_table(OTU3, taxa_are_rows = TRUE), tax_table(biom), sample_data(s3))
samplephylo = phyloseq(otu_table(OTU, taxa_are_rows = TRUE), tax_table(biom), sample_data(samplephylo))
samplephylonn = phyloseq(otu_table(OTUnn, taxa_are_rows = TRUE), tax_table(biom), sample_data(samplephylo))
Phyloseq1nn = phyloseq(otu_table(OTU1nn, taxa_are_rows = TRUE), tax_table(biom), sample_data(s1nn))
Phyloseq2nn = phyloseq(otu_table(OTU2nn, taxa_are_rows = TRUE), tax_table(biom), sample_data(s2nn))
Phyloseq3nn = phyloseq(otu_table(OTU3nn, taxa_are_rows = TRUE), tax_table(biom), sample_data(s3nn))
#rm(s1,s2,s3,s1nn,s2nn,s3nn,Neg)


glomphy = tax_glom(samplephylo, taxrank = "Class")
glomOTU = as.data.frame(glomphy@otu_table@.Data)
abOTU = as.data.frame(abundances(glomOTU, transform = 'compositional'))
abOTU = data.matrix(t(abOTU))

plot_bar(Phyloseq1, fill = "Class")
plot_bar(Phyloseq2, fill = "Class")
ggplot(phyloseq(otu_table(glomOTU, taxa_are_rows = TRUE))) 
species_counts <- data.frame(
  Sample = colnames(OTU),
  TotalSpecies = colSums(OTU))
ggplot(data = species_counts, aes(x = Sample, y = TotalSpecies)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18),  
    legend.title = element_text(size = 18), 
    legend.key.size = unit(1.5, "lines"))  

totalread = data.frame(tr =colSums(OTU))
abOTU_df = as.data.frame(abOTU)
abOTU_df$sample_name = rownames(abOTU_df)
abOTU_long = pivot_longer(
  data = abOTU_df,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU_long$Virus_Class = tax[abOTU_long$Virus_Class, 'Class']
ggplot(data = abOTU_long, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample Name", y = "Relative Abundance (%)") + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=18)) +
  theme(
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    legend.text = element_text(size = 30),  
    legend.title = element_text(size = 30), 
    legend.key.size = unit(1.5, "lines"))    


  





# Plot library size 
df = as.data.frame(sample_data(samplephylo))
df$LibrarySize = sample_sums(samplephylo)
df = df[order(df$LibrarySize),]
df$Paziente <- factor(df$Paziente, levels = df$Paziente[order(df$LibrarySize)])
ggplot(data=df, aes(x=Paziente, y=LibrarySize, color=SoC)) + geom_point() + 
  geom_point(size = 6) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 21),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    legend.text = element_text(size = 30),  
    legend.title = element_text(size = 30), 
    legend.key.size = unit(1.5, "lines"))    
                                 
      