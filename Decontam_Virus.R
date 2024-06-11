dectm_output = isContaminant(
  normphylo,
  neg = normphylo@sam_data$SoC=='Control',
  method = "prevalence",
  threshold = 0.3, 
  batch = normphylo@sam_data$Batch ,
  normalize = TRUE,
  detailed = TRUE
)  
ggplot(dectm_output, aes(x=p)) + 
  geom_density(fill="#69b3a2", color="#e9ecef") 

ggplot(dectm_output, aes(x=p,fill=as.factor(prev)) )+ 
  geom_histogram(color="#e9ecef") + 
  geom_density(fill="#69b3a2", color="#e9ecef",alpha=0.4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


dectm_output$prevalence_category <- ifelse(dectm_output$prev == 2, "2",
                                ifelse(dectm_output$prev >= 3 & dectm_output$prev <= 5, "3-5",
                                       ifelse(dectm_output$prev >= 6 & dectm_output$prev <= 10, "6-10",
                                              "11+")))


ggplot(dectm_output, aes(x = p, fill = prevalence_category))  +
  geom_histogram(color = "#e9ecef", alpha = 0.5) +
  geom_density(fill="#69b3a2", color="#e9ecef",alpha=0.4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25),  
    legend.title = element_text(size = 25), 
    legend.key.size = unit(1.5, "lines"))
  
  
  
  
  
  
  
  
  
  
  

iscont = subset(dectm_output, dectm_output$contaminant==FALSE)
iscont = as.data.frame(tax[rownames(tax) %in% rownames(iscont), ])
OTUcontdectm =  normOTU[row.names(iscont),]
phylobrutecontdectm = phyloseq(otu_table(OTUcontdectm, taxa_are_rows = TRUE), tax_table(biom), sample_data(data))
glomphycontdectm = tax_glom(phylobrutecontdectm, taxrank = "Class")
glomOTUcontdectm = as.data.frame(glomphycontdectm@otu_table@.Data)
abOTUcontdectm = as.data.frame(abundances(glomOTUcontdectm, transform = 'compositional'))
abOTUcontdectm = data.matrix(t(abOTUcontdectm))
abOTU_dfcontdectm = as.data.frame(abOTUcontdectm)
abOTU_dfcontdectm$sample_name = rownames(abOTU_dfcontdectm)
abOTU_longcontdectm = pivot_longer(
  data = abOTU_dfcontdectm,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU_longcontdectm$Virus_Class = tax[abOTU_longcontdectm$Virus_Class, 'Class']

ggplot(data = abOTU_longcontdectm, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
  geom_bar(stat = "identity") + 
  labs(x = "Sample Name", y = "Relative Abundance (%)") + 
  scale_y_continuous(labels = scales::percent) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 21),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    legend.text = element_text(size = 30),  
    legend.title = element_text(size = 30), 
    legend.key.size = unit(1.5, "lines"))    




# Create a new column in dectm_output that maps prev values to categories

isnot = isNotContaminant(
  normphylo,
  neg = normphylo@sam_data$SoC == "Control",
  method = "prevalence",
  threshold = 0.5,
  normalize = FALSE,
  detailed = TRUE
)
ggplot(isnot, aes(x=not.contaminant)) + 
  geom_bar(fill="#69b3a2", color="#e9ecef") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
isnotcont = subset(isnot, isnot$not.contaminant==TRUE)
isnotcont = as.data.frame(tax[rownames(tax) %in% rownames(isnotcont), ])
OTUdectm =  normOTU[row.names(isnotcont),]
OTUdectm =  OTUdectm[ , colSums(OTUdectm) > 0]
phylobrutedectm = phyloseq(otu_table(OTUdectm, taxa_are_rows = TRUE), tax_table(biom), sample_data(data))


glomphydectm = tax_glom(phylobrutedectm, taxrank = "Class")
glomOTUdectm = as.data.frame(glomphydectm@otu_table@.Data)
abOTUdectm = as.data.frame(abundances(glomOTUdectm, transform = 'compositional'))
abOTUdectm = data.matrix(t(abOTUdectm))
abOTU_dfdectm = as.data.frame(abOTUdectm)
abOTU_dfdectm$sample_name = rownames(abOTU_dfdectm)
abOTU_longdectm = pivot_longer(
  data = abOTU_dfdectm,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU_longdectm$Virus_Class = tax[abOTU_longdectm$Virus_Class, 'Class']

ggplot(data = abOTU_longdectm, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
  geom_bar(stat = "identity") + 
  labs(x = "Sample Name", y = "Relative Abundance (%)") + 
  scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 21),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    legend.text = element_text(size = 30),  
    legend.title = element_text(size = 30), 
    legend.key.size = unit(1.5, "lines"))  +
scale_fill_manual(values = c("#cd9500", "#00be68", "#00a8ff", rep("gray", length(unique(abOTU_longdectm$Virus_Class)) - 3)))



ggplot(data = FilabOTU_long3, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
  geom_bar(stat = "identity") + 
  labs(x = "Sample Name", y = "Relative Abundance (%)") + 
  scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25),  
    legend.title = element_text(size = 25), 
    legend.key.size = unit(1.5, "lines")
  ) +
  scale_fill_manual(values = c("#cd9500", "#00be68", rep("gray", length(unique(FilabOTU_long3$Virus_Class)) - 2)))




















ggplot(data = abOTU_longdectm, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
  geom_bar(stat = "identity") + 
  labs(x = "Sample Name", y = "Relative Abundance (%)") + 
  scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25),  
    legend.title = element_text(size = 25), 
    legend.key.size = unit(1.5, "lines")) +
  








## parte vecchia del codice con il biom batteri 
genusnc = unique(isnotcont$Genus)
Burkholderia = nrow(subset(tax, tax$Genus == 'Burkholderia'))/nrow(tax)*100
Paraburkholderia = nrow(subset(tax, tax$Genus == 'Paraburkholderia'))/nrow(tax)*100
Acidovorax = nrow(subset(tax, tax$Genus == 'Acidovorax'))/nrow(tax)*100
Mesorhizobium = nrow(subset(tax, tax$Genus == 'Mesorhizobium'))/nrow(tax)*100
Lawsonella = nrow(subset(tax, tax$Genus == 'Lawsonella'))/nrow(tax)*100
genus = as.data.frame(Burkholderia, Paraburkholderia, Acidovorax, Mesorhizobium, Lawsonella)
names(genus) = c("Percentage")
rownames(genus) = c("Burkholderia", "Paraburkholderia", "Acidovorax", "Mesorhizobium", "Lawsonella")
table(dectm_output$contaminant)
table(isnot$not.contaminant)
dectm_species = subset(dectm_output, dectm_output$contaminant==TRUE)
taxa_to_keep = subset(dectm_output, dectm_output$contaminant==FALSE)
sample_filter_dectm = prune_taxa(rownames(taxa_to_keep), samplenn)
dectm_species = as.data.frame(tax[rownames(tax) %in% rownames(dectm_species), ])

 
order1$Var1 = as.factor(order1$Var1 )
order2$Var1 = as.factor(order2$Var1)
order3$Var1 = as.factor(order3$Var1)
dectm_species$Species = as.factor(dectm_species$Species)
ven = list(order1$Var1,order2$Var1,order3$Var1, dectm_species$Species)
names(ven) = c("Order1","Order2","Order3","Dectm")
ggVennDiagram(ven) 

df <- as.data.frame(sample_data(normphylo))
df$LibrarySize <- sample_sums(normphylo)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=SoC)) + geom_point()


contamdf.freq <- isContaminant(normphylo, method="frequency", conc="Quantification")
head(contamdf.freq)
set.seed(100)
plot_frequency(normphylo, taxa_names(normphylo)[sample(which(contamdf.freq$contaminant),3)], conc="Quantification") 





















