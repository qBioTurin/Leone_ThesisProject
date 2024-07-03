# 1) Calcolare la relative abundance dei virus per la intera otu 
tanOTU =abundances(samplephylonn, transform = 'compositional')
tanOTUbin = tanOTU
tanOTUbin[tanOTUbin > 0.001] <- 1
tanOTUbin[tanOTUbin <= 0.001] <- 0
rs = data.frame(Sum=rowSums(tanOTUbin),rn=rownames(tanOTUbin) )
rowsumtanOTUbin = rs[rs$Sum != 0, ]
tanOTU = OTUnn[rownames(rowsumtanOTUbin),]
tanOTU1 = tanOTU[,s1nn$Paziente]
tanOTU2 = tanOTU[,s2nn$Paziente]
tanOTU3 = tanOTU[,s3nn$Paziente]
tanphylo1 = phyloseq(otu_table(tanOTU1, taxa_are_rows = TRUE), sample_data(s1nn) , tax_table(biom@tax_table))
tanphylo2 = phyloseq(otu_table(tanOTU2, taxa_are_rows = TRUE), sample_data(s2nn) , tax_table(biom@tax_table))
tanphylo3 = phyloseq(otu_table(tanOTU3, taxa_are_rows = TRUE), sample_data(s3nn) , tax_table(biom@tax_table))

# 2) Calcolo la prevalenza e filtro quelli con prevalenza maggiore del 25% Consiglio di filtrare dopo
Pv1 = data.frame(Prev = prevalence(tanphylo1),tax_table(tanphylo1))
Pv1tf = Pv1[Pv1$Prev > 0.25,]
Pv2 = data.frame(Prev = prevalence(tanphylo2),tax_table(tanphylo2))
Pv2tf = Pv2[Pv2$Prev > 0.25,]
Pv3 = data.frame(Prev = prevalence(tanphylo3),tax_table(tanphylo3))
Pv3tf = Pv3[Pv3$Prev > 0.25,]

# 3) Calcolo il rapporto tra prev massima e minima
dfrank = data.frame(row.names = row.names(tanOTU), Pv1 = as.numeric(Pv1$Prev), Pv2 = as.numeric(Pv2$Prev), Pv3 = as.numeric(Pv3$Prev))
dfrank$Max = pmax(dfrank$Pv1, dfrank$Pv2, dfrank$Pv3)
dfrank$Min = pmin(dfrank$Pv1, dfrank$Pv2, dfrank$Pv3)
dfrankfilt = dfrank[dfrank$Min == 0 & dfrank$Max > 0.25, ]
dfrankfilt$foldmaxmin = dfrankfilt$Max / dfrankfilt$Min
dfrank$foldmaxmin = dfrank$Max / dfrank$Min
dfrankfinite = dfrank[is.finite(dfrank$foldmaxmin), ]
dfrankfinite = dfrankfinite[dfrankfinite$foldmaxmin >= 2, ]
dfrankfinitefilt = rbind(dfrankfinite, dfrankfilt)

# se min e med sono 0 e il max è maggiore di zero li tengo 


# Interseco con Pvtf e prendo solo quelli che hanno prevalenza maggiore di 0.25
#contPv1 = dfrankfinitefilt[row.names(dfrankfinitefilt) %in% row.names(Pv1tf), ] # qui ho preso quelli che hanno maggiore di 0.25 e che hanno il doppio della fold (max/min), essendo che tutte le maxx appartengono a pv1
#contPv2 = dfrankfinitefilt[row.names(dfrankfinitefilt) %in% row.names(Pv2tf), ]
#contPv3 = dfrankfinitefilt[row.names(dfrankfinitefilt) %in% row.names(Pv3tf), ]
rntr1 = c("2956710","2956237","2169967")
rntr3 = c("202910")
NcontPv1= Pv1[!(row.names(Pv1) %in% rntr1),]
NcontPv2= Pv2
NcontPv3 = Pv3[!(row.names(Pv3) %in% rntr3),]
contOTU1 = OTU1[rntr1,]
contOTU3 = OTU3[rntr3,]

clr1 = as.data.frame(clr(OTU1nn[row.names(NcontPv1),])) 
clr3 = as.data.frame(clr(OTU3nn[row.names(NcontPv3),]))

clr1cont = as.data.frame(clr(contOTU1))
clr3cont =  bind_rows((clr(contOTU3))) # qui li ho fatti cosi perche non so quale classe sia meglio utilizzare
clr1cont = subset(clr1cont, select = -Neg1)
clr3cont = subset(clr3cont, select = -Neg3)
## ciclo


significant_rows <- data.frame(otu = character(0), estimate = numeric(0), p.value = numeric(0))


for (i in 1:nrow(clr3)) {
  a <- as.numeric(clr3[i,])
  
  
  for (x in 1:nrow(clr3cont)) { 
    b <- as.numeric(clr3cont[x,])
    
    
    if (sd(a) != 0 & sd(b) != 0) {
      
      result <- cor.test(a, b, method = "spearman", exact = FALSE)
      f <-as.numeric(result$estimate)
      g <- result$p.value
      
      if (f > 0.7 & g < 0.05) {
      significant_rows_tmp = data.frame(otu = rownames(clr1)[i], estimate = f, p.value = g)
      significant_rows <- rbind(significant_rows, significant_rows_tmp)
      #}
      
    }
    
  }
  
}

print(significant_rows)
# la condizione sd(a) e sd(b) serve perche se no errore Error in if (f > 0.5 & g < 0.05) { : argument is of length zero In addition: Warning message: In cor(rank(x), rank(y)) : the standard deviation is zero

# When calculating correlations (such as Spearman rank correlation), it’s essential to ensure that both vectors being compared have variability (i.e., they are not constant). If a vector has zero standard deviation (meaning all its values are the same), it cannot be used for correlation calculations because it lacks variability

significant_rows3 <- character(0)
for (k in 1:nrow(clr3)) {
  a3 <- as.numeric(clr3[k,])
  
  for (l in 1:nrow(clr3cont)) {
    b3 <- as.numeric(clr3cont[l,])
    if (sd(a3) != 0 & sd(b3) != 0) {
      result3 <- cor.test(a3, b3, method = "spearman", exact = FALSE)
      f3 <- na.omit(as.numeric(result3$estimate))
      g3 <- na.omit(result3$p.value)
      if (f3 > 0.7 & g3 < 0.05) {
        significant_rows3 <- c(significant_rows3, rownames(clr3)[k])
      }
    }
  }
}

print(significant_rows3)

# la condizione sd(a) e sd(b) serve perche se no errore Error in if (f > 0.5 & g < 0.05) { : argument is of length zero In addition: Warning message: In cor(rank(x), rank(y)) : the standard deviation is zero

# When calculating correlations (such as Spearman rank correlation), it’s essential to ensure that both vectors being compared have variability (i.e., they are not constant). If a vector has zero standard deviation (meaning all its values are the same), it cannot be used for correlation calculations because it lacks variability

significant_rows2 <- character(0)
for (k in 1:nrow(clr2)) {
  a2 <- as.numeric(clr2[k,])
  
  for (l in 1:nrow(clr2cont)) {
    b2 <- as.numeric(clr2cont[l,])
    if (sd(a2) != 0 & sd(b2) != 0) {
      result2 <- cor.test(a2, b2, method = "spearman", exact = FALSE)
      f2 <- na.omit(as.numeric(result2$estimate))
      g2 <- na.omit(result2$p.value)
      if (f2 > 0.7 & g2 < 0.05) {
        significant_rows2 <- c(significant_rows2, rownames(clr2)[k])
      }
    }
  }
}

print(significant_rows2)






# la condizione sd(a) e sd(b) serve perche se no errore Error in if (f > 0.5 & g < 0.05) { : argument is of length zero In addition: Warning message: In cor(rank(x), rank(y)) : the standard deviation is zero

# When calculating correlations (such as Spearman rank correlation), it’s essential to ensure that both vectors being compared have variability (i.e., they are not constant). If a vector has zero standard deviation (meaning all its values are the same), it cannot be used for correlation calculations because it lacks variability

significant_rows1 <- character(0)
for (k in 1:nrow(clr1)) {
  a1 <- as.numeric(clr1[k,])
  
  for (l in 1:nrow(clr1cont)) {
    b1 <- as.numeric(clr1cont[l,])
    if (sd(a1) != 0 & sd(b1) != 0) {
      result1 <- cor.test(a1, b1, method = "spearman", exact = FALSE)
      f1 <- na.omit(as.numeric(result1$estimate))
      g1 <- na.omit(result1$p.value)
      if (f1 > 0.7 & g1 < 0.05) {
        significant_rows1 <- c(significant_rows1, rownames(clr1)[k])
      }
    }
  }
}

print(significant_rows1)
tanfinal1 = normOTU1[row.names(normOTU1) %in% row.names(rowsumtanOTUbin),]
tanfinal1 = tanfinal1[!rownames(tanfinal1) %in% rownames(contOTU1),]
tanfinal1 = tanfinal1[!rownames(tanfinal1) %in% significant_rows1,]




tanfinal2 = normOTU2[row.names(rowsumtanOTUbin),]

tanfinal3 = normOTU3[row.names(OTU3) %in% row.names(rowsumtanOTUbin),]
tanfinal3 = tanfinal3[!rownames(tanfinal3) %in% rownames(contOTU3),]
tanfinal3 = tanfinal3[!rownames(tanfinal3) %in% significant_rows3,]


normOTUNN = round(normOTUNN)
filOTUnn = normOTUNN[!rowSums(normOTU > 100),]
filOTUnn = filOTUnn[rowSums(filOTUnn)!=0,]
contaminantitan2 = normOTU2[row.names(normOTU2) %in% row.names(filOTUnn),]
contaminantitan3 = normOTU3[row.names(normOTU3) %in% row.names(filOTUnn),]
contaminantitan3 = rbind(contaminantitan3, normOTU3['202910',])
contaminantitan3 = rbind(contaminantitan3, normOTU3['2843645',])
contamiantitan1 = normOTU1[row.names(normOTU1) %in% row.names(filOTUnn),]
contamiantitan1 = rbind(contamiantitan1, normOTU1['2956237',])
contamiantitan1 = rbind(contamiantitan1, normOTU1['2169967',])
contamiantitan1 = contamiantitan1[rowSums(contamiantitan1)!=0,]
contaminantitan2 = contaminantitan2[rowSums(contaminantitan2)!=0,]
contaminantitan3 = contaminantitan3[rowSums(contaminantitan3)!=0,]
  

filOTU1nn = subset(filOTU1nn, select = -Neg1)
filOTU2nn = subset(filOTU2nn, select = -Neg2)
filOTU3nn = subset(filOTU3nn, select = -Neg3)
filOTU1nnNorm = normOTU1NN[!rownames(normOTU1) %in% rownames(contamiantitan1),]
filOTU2nnNorm = normOTU2NN[!rownames(normOTU2) %in% rownames(contaminantitan2),]
filOTU3nnNorm = normOTU3NN[!rownames(normOTU3) %in% rownames(contaminantitan3),]

FilPhyloseq1nn = phyloseq(otu_table(filOTU1nn, taxa_are_rows = TRUE), tax_table(biom), sample_data(s1nn))
FilPhyloseq2nn = phyloseq(otu_table(filOTU2nn, taxa_are_rows = TRUE), tax_table(biom), sample_data(s2nn))
FilPhyloseq3nn = phyloseq(otu_table(filOTU3nn, taxa_are_rows = TRUE), tax_table(biom), sample_data(s3nn))

tancontOTU1 =  normOTU1[!(rownames(normOTU1) %in% rownames(filOTU1nn)),]
tancontOTU2 =  normOTU2[!(rownames(normOTU2) %in% rownames(filOTU2nn)),]
tancontOTU3 =  normOTU3[!(rownames(normOTU3) %in% rownames(filOTU3nn)),]
tancontOTU1 = tancontOTU1[rowSums(tancontOTU1)!=0,]
tancontOTU2 = tancontOTU2[,colSums(tancontOTU2)!=0]
tancontOTU3 = tancontOTU3[,colSums(tancontOTU3)!=0]



intersect(intersect(intersect(rownames(tancontOTU1), rownames(tancontOTU2)), rownames(tancontOTU3)) )





Fil1glomphy = tax_glom(FilPhyloseq1nn, taxrank = "Class")
FilglomOTU1 = as.data.frame(Fil1glomphy@otu_table@.Data)
FilabOTU1 = as.data.frame(abundances(FilglomOTU1, transform = 'compositional'))
FilabOTU1 = data.matrix(t(FilabOTU1))
FilabOTU_df1 = as.data.frame(FilabOTU1)
FilabOTU_df1$sample_name = rownames(FilabOTU_df1)
FilabOTU_long1 = pivot_longer(
  data = FilabOTU_df1,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
FilabOTU_long1$Virus_Class = tax[FilabOTU_long1$Virus_Class, 'Class']

ggplot(data = FilabOTU_long1, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
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
  scale_fill_manual(values = c("#cd9500", "#00be68", rep("gray", length(unique(FilabOTU_long3$Virus_Class)) - 2)))



Fil3glomphy = tax_glom(FilPhyloseq3nn, taxrank = "Class")
FilglomOTU3 = as.data.frame(Fil3glomphy@otu_table@.Data)
FilabOTU3 = as.data.frame(abundances(FilglomOTU3, transform = 'compositional'))
FilabOTU3 = data.matrix(t(FilabOTU3))
FilabOTU_df3 = as.data.frame(FilabOTU3)
FilabOTU_df3$sample_name = rownames(FilabOTU_df3)
FilabOTU_long3 = pivot_longer(
  data = FilabOTU_df3,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
FilabOTU_long3$Virus_Class = tax[FilabOTU_long3$Virus_Class, 'Class']

ggplot(data = FilabOTU_long3, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
  geom_bar(stat = "identity",) + 
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
  scale_fill_manual(values = c("#cd9500", "#00be68", rep("gray", length(unique(FilabOTU_long3$Virus_Class)) - 2)))


Fil2glomphy = tax_glom(FilPhyloseq2nn, taxrank = "Class")
FilglomOTU2 = as.data.frame(Fil2glomphy@otu_table@.Data)
FilabOTU2 = as.data.frame(abundances(FilglomOTU2, transform = 'compositional'))
FilabOTU2 = data.matrix(t(FilabOTU2))
FilabOTU_df2 = as.data.frame(FilabOTU2)
FilabOTU_df2$sample_name = rownames(FilabOTU_df2)
FilabOTU_long2 = pivot_longer(
  data = FilabOTU_df2,
  cols = -sample_name, 
  names_to = "Virus_Class", 
  values_to = "relative_abundance" 
)
FilabOTU_long2$Virus_Class = tax[FilabOTU_long2$Virus_Class, 'Class']

ggplot(data = FilabOTU_long2, aes(x = sample_name, y = relative_abundance, fill = Virus_Class)) +
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
