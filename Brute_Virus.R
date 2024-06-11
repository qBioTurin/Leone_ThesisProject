tax$rn = rownames(tax)

brut_sum1 = normOTU1[normOTU1$Neg1 != 0, ]
brut_sum1$rn = rownames(brut_sum1)
merged_df1 <- inner_join(tax, brut_sum1, by = "rn")
rownames(merged_df1) = merged_df1$rn
brut_sum1 = subset(brut_sum1, select = -rn) 
merged_df1 = merged_df1[,1:7]
OTUbrute1 = normOTU1[normOTU1$Neg1 == 0, ]
OTUbrute1 =  OTUbrute1[, !names(OTUbrute1) %in% "Neg1"]
phylobrute1 = phyloseq(otu_table(OTUbrute1, taxa_are_rows = TRUE), tax_table(biom), sample_data(s1))


brut_sum2 = normOTU2[normOTU2$Neg2 != 0, ]
brut_sum2$rn = rownames(brut_sum2)
merged_df2 <- inner_join(tax, brut_sum2, by = "rn")
rownames(merged_df2) = merged_df2$rn
merged_df2 = merged_df2[,1:7]
OTUbrute2 = normOTU2[normOTU2$Neg2 == 0, ]
OTUbrute2 =  OTUbrute2[, !names(OTUbrute2) %in% "Neg2"]
phylobrute2 = phyloseq(otu_table(OTUbrute2, taxa_are_rows = TRUE), tax_table(biom), sample_data(s2))

brut_sum3 = normOTU3[normOTU3$Neg3 != 0, ]
brut_sum3$rn = rownames(brut_sum3)
merged_df3 <- inner_join(tax, brut_sum3, by = "rn")
rownames(merged_df3) = merged_df3$rn
merged_df3 = merged_df3[,1:7]
OTUbrute3 = normOTU3[normOTU3$Neg3 == 0, ]
OTUbrute3 =  OTUbrute3[, !names(OTUbrute3) %in% "Neg3"]
phylobrute3 = phyloseq(otu_table(OTUbrute3, taxa_are_rows = TRUE), tax_table(biom), sample_data(s3))

ven = list(rownames(merged_df1), rownames(merged_df2), rownames(merged_df3))
ggVennDiagram(ven) + ggplot2::scale_fill_gradient(low="#99BDD5",high = "#044E7E")
intersect((intersect(row.names(merged_df3), row.names(merged_df2))), row.names(merged_df1))

noncontbrutOTU1 = OTUbrute1[,colSums(OTUbrute1)!=0 ]
noncontbrutOTU2 = OTUbrute2[,colSums(OTUbrute2)!=0 ]
noncontbrutOTU3 = OTUbrute3[,colSums(OTUbrute3)!=0 ]
noncontbrutOTU1phylo = phyloseq(otu_table(noncontbrutOTU1, taxa_are_rows = TRUE), tax_table(biom), sample_data(s1))
noncontbrutOTU2phylo = phyloseq(otu_table(noncontbrutOTU2, taxa_are_rows = TRUE), tax_table(biom), sample_data(s2))
noncontbrutOTU3phylo = phyloseq(otu_table(noncontbrutOTU3, taxa_are_rows = TRUE), tax_table(biom), sample_data(s3))



mergedOTU1 = normOTU1NN[,!colSums(OTUbrute1!=0) ]
mergedOTU2 = normOTU2NN[,!colSums(OTUbrute2!=0) ]
mergedOTU3 = normOTU3NN[,!colSums(OTUbrute3!=0) ]
mergedphylo1 = phyloseq(otu_table(mergedOTU1, taxa_are_rows = TRUE), tax_table(biom), sample_data(s1))
mergedphylo2 = phyloseq(otu_table(mergedOTU2, taxa_are_rows = TRUE), tax_table(biom), sample_data(s2))
mergedphylo3 = phyloseq(otu_table(mergedOTU3, taxa_are_rows = TRUE), tax_table(biom), sample_data(s3))

# Plot batch1
glomphy1 = tax_glom(phylobrute1, taxrank = "Class")
glomOTU1 = as.data.frame(glomphy1@otu_table@.Data)
abOTU1 = as.data.frame(abundances(glomOTU1, transform = 'compositional'))
abOTU1 = data.matrix(t(abOTU1))
abOTU_df1 = as.data.frame(abOTU1)
abOTU_df1$sample_name = rownames(abOTU_df1)
abOTU_long1 = pivot_longer(
  data = abOTU_df1,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU_long1$Virus_Class = tax[abOTU_long1$Virus_Class, 'Class']

ggplot(data = abOTU_long1, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
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
    legend.key.size = unit(1.5, "lines"))

# Plot batch2
glomphy2 = tax_glom(phylobrute2, taxrank = "Class")
glomOTU2 = as.data.frame(glomphy2@otu_table@.Data)
abOTU2 = as.data.frame(abundances(glomOTU2, transform = 'compositional'))
abOTU2 = data.matrix(t(abOTU2))
abOTU_df2 = as.data.frame(abOTU2)
abOTU_df2$sample_name = rownames(abOTU_df2)
abOTU_long2 = pivot_longer(
  data = abOTU_df2,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU_long2$Virus_Class = tax[abOTU_long2$Virus_Class, 'Class']

ggplot(data = abOTU_long2, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
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
  scale_fill_brewer(palette = "Set2")

# Plot batch3

glomphy3 = tax_glom(phylobrute3, taxrank = "Class")
glomOTU3 = as.data.frame(glomphy3@otu_table@.Data)
abOTU3 = as.data.frame(abundances(glomOTU3, transform = 'compositional'))
abOTU3 = data.matrix(t(abOTU3))
abOTU_df3 = as.data.frame(abOTU3)
abOTU_df3$sample_name = rownames(abOTU_df3)
abOTU_long3 = pivot_longer(
  data = abOTU_df3,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU_long3$Virus_Class = tax[abOTU_long3$Virus_Class, 'Class']

ggplot(data = abOTU_long3, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
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
    legend.key.size = unit(1.5, "lines"))

## I sample vuoti sono sia quelli che hanno 0 reads che quelli che sono stati rimossi da bruteforce
## Plot mergedotu 

#1
glomphymerged1 = tax_glom(noncontbrutOTU1phylo, taxrank = "Class")
glomOTU1merged = as.data.frame(glomphymerged1@otu_table@.Data)
abOTU1merged = as.data.frame(abundances(glomOTU1merged, transform = 'compositional'))
abOTU1merged = data.matrix(t(abOTU1merged))
abOTU1merged_df1 = as.data.frame(abOTU1merged)
abOTU1merged_df1$sample_name = rownames(abOTU1merged)
abOTU1merged_long1 = pivot_longer(
  data = abOTU1merged_df1,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU1merged_long1$Virus_Class = tax[abOTU1merged_long1$Virus_Class, 'Class']

ggplot(data = abOTU1merged_long1, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
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
    legend.key.size = unit(1.5, "lines"))

# 2
glomphymerged2 = tax_glom(noncontbrutOTU2phylo, taxrank = "Class")
glomOTU2merged = as.data.frame(glomphymerged2@otu_table@.Data)
abOTU2merged = as.data.frame(abundances(glomOTU2merged, transform = 'compositional'))
abOTU2merged = data.matrix(t(abOTU2merged))
abOTU2merged_df2 = as.data.frame(abOTU2merged)
abOTU2merged_df2$sample_name = rownames(abOTU2merged)
abOTU2merged_long2 = pivot_longer(
  data = abOTU2merged_df2,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU2merged_long2$Virus_Class = tax[abOTU2merged_long2$Virus_Class, 'Class']

ggplot(data = abOTU2merged_long2, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
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
    legend.key.size = unit(1.5, "lines"))
# 3
glomphymerged3 = tax_glom(noncontbrutOTU3phylo, taxrank = "Class")
glomOTU3merged = as.data.frame(glomphymerged3@otu_table@.Data)
abOTU3merged = as.data.frame(abundances(glomOTU3merged, transform = 'compositional'))
abOTU3merged = data.matrix(t(abOTU3merged))
abOTU3merged_df3 = as.data.frame(abOTU3merged)
abOTU3merged_df3$sample_name = rownames(abOTU3merged)
abOTU3merged_long3 = pivot_longer(
  data = abOTU3merged_df3,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
abOTU3merged_long3$Virus_Class = tax[abOTU3merged_long3$Virus_Class, 'Class']

ggplot(data = abOTU3merged_long3, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
  geom_bar(stat = "identity",s ) + 
  labs(x = "Sample Name", y = "Relative Abundance (%)") + 
  scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25),  
    legend.title = element_text(size = 25), 
    legend.key.size = unit(1.5, "lines"))
