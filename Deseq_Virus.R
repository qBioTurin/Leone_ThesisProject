# differential gene expression (deseq2)
library( "DESeq2" )
library(ggplot2)
countData = as.matrix(OTU)
data = sample
data$SoCnum[data$SoC == "Sample"] = "FALSE"
data$SoCnum[data$SoC == "Control"] = "TRUE"
data$SoCnum = as.factor(data$SoCnum)
data = data[!(sample$Paziente %in% c(172,406, 178, 174, 184, 420, 464, 386)),]
data$Batch[c(76,52,53)] = c(1,2,3)
rownames(data) = data$Paziente
dds = DESeqDataSetFromMatrix(countData=OTU, 
                            colData=data, 
                             design=~SoC)

dds = DESeq(dds)
res = results(dds)
summary(res)
vsdata = varianceStabilizingTransformation (dds, blind=FALSE)
plotPCA(vsdata, intgroup="SoC")
normOTU = as.data.frame(counts(dds, normalized=TRUE))
normphylo = phyloseq(otu_table(normOTU, taxa_are_rows = TRUE), sample_data(data) , tax_table(biom@tax_table))
rsnormOTU = data.frame(colSums(normOTU))
rsnormOTU=rownames_to_column(rsnormOTU, var = "RowNames")
rsnormOTU_long <- rsnormOTU %>%
  pivot_longer(cols = -RowNames, names_to = "ColumnNames", values_to = "Values")
ggplot(data = rsnormOTU_long, aes(x = RowNames, y = Values)) +
  geom_point() +
  facet_wrap(~ColumnNames, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
rsnormOTU_sorted <- rsnormOTU_long %>%
  arrange(Values)


library(ggplot2)

ggplot(data = rsnormOTU_long, aes(x = RowNames, y = Values)) +
  geom_point(size = 6, color = "#3793D1") +  # Increase dot size and change color
  scale_y_continuous(trans = 'log10') +
  labs(x = "Sample Names", y = "Log(Normalized Reads)") +  # Set axis labels
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.title.x = element_text(size = 18),  # Increase x-axis title size
    axis.title.y = element_text(size = 18)   # Increase y-axis title size
  )

library(ggplot2)

# Create a new color column based on RowNames
rsnormOTU_long$Color <- ifelse(rsnormOTU_long$RowNames %in% c("Neg1", "Neg3", "Neg2"), "orange", "#3793D1")

# Reorder RowNames based on Values
rsnormOTU_long$RowNames <- factor(rsnormOTU_long$RowNames, levels = rsnormOTU_long$RowNames[order(rsnormOTU_long$Values)])

# Plot using ggplot2 with ordered RowNames and custom colors
ggplot(data = rsnormOTU_long, aes(x = RowNames, y = Values, color = Color)) +
  geom_point(size = 6) +  # Increase dot size
  scale_y_continuous(trans = 'log10') +
  labs(x = "Sample Names", y = "Log(Normalized Reads)") +  # Set axis labels
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18)) +
  
  scale_color_identity()  # Use the color values directly +
theme(
  axis.text.y = element_text(size = 30),
  axis.title = element_text(size = 30),
  legend.text = element_text(size = 30),  
  legend.title = element_text(size = 30), 
  legend.key.size = unit(1.5, "lines"))    



























ggplot(data = rsnormOTU_long_sorted, aes(x = RowNames, y = Values)) +
  geom_point() +
  scale_y_continuous(trans = 'log10') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# separo i tre phyloseq nei batch 

normOTU1 = normOTU[,s1$Paziente] 
normOTU2 = normOTU[,s2$Paziente]
normOTU3 = normOTU[,s3$Paziente]
normphylo1 = phyloseq(otu_table(normOTU1, taxa_are_rows = TRUE), sample_data(s1) , tax_table(biom@tax_table))                      
normphylo2 = phyloseq(otu_table(normOTU2, taxa_are_rows = TRUE), sample_data(s2) , tax_table(biom@tax_table))
normphylo3 = phyloseq(otu_table(normOTU3, taxa_are_rows = TRUE), sample_data(s3) , tax_table(biom@tax_table))
normOTU1NN = subset(normOTU1, select = -23) 
normOTU2NN = subset(normOTU2, select = -27)
normOTU3NN = subset(normOTU3, select = -26)
normphylo1NN = phyloseq(otu_table(normOTU1NN, taxa_are_rows = TRUE), sample_data(s1) , tax_table(biom@tax_table))
normphylo2NN = phyloseq(otu_table(normOTU2NN, taxa_are_rows = TRUE), sample_data(s2) , tax_table(biom@tax_table))
normphylo3NN = phyloseq(otu_table(normOTU3NN, taxa_are_rows = TRUE), sample_data(s3) , tax_table(biom@tax_table))

normphyloab = abundances(normphylo, transform = 'compositional')
