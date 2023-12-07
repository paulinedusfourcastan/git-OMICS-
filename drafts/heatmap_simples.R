# Chargement des librairies 
library(OmicCircos)
library(heatmaply) 


# Chargement des donnees 
data("TCGA.PAM50_genefu_hg18")
data("TCGA.BC.fus")
data("TCGA.BC.cnv.2k.60")
data("TCGA.BC.gene.exp.2k.60")
data("TCGA.BC.sample60")
data("TCGA.BC_Her2_cnv_exp") 


# Arrangement des donnees 
Her2.i = which(TCGA.BC.sample60[,2] == "Her2")
Her2.n = TCGA.BC.sample60[Her2.i,1]
Her2.j = which(colnames(TCGA.BC.cnv.2k.60)%in%Her2.n)
cnv = TCGA.BC.cnv.2k.60[,c(1:3,Her2.j)]

cnv.m = cnv[,c(4:ncol(cnv))]
cnv.m[cnv.m > 2] = 2
cnv.m[cnv.m < -2] = -2
cnv = cbind(cnv[,1:3], cnv.m)

Her2.j = which(colnames(TCGA.BC.gene.exp.2k.60)%in% Her2.n)

gene.exp = TCGA.BC.gene.exp.2k.60[,c(1:3,Her2.j)]



## Affichage des heatmap ##

# Choix des couleurs 
gradient_col <- ggplot2::scale_fill_gradient2(
  low = "blue", high = "red", midpoint = 0, limits = c(-8, 8))




## Chromosome 1 ##
chr1 <- subset(gene.exp, gene.exp$chr == 1 ) # on selectionne le chromosome souhaité dans notre dataframe
mydata_1 <- chr1[,3:18] # on selectionne les colonnes qui nous interessent
rownames(mydata_1) <- mydata_1[,1] # nom des lignes 
mydata_1 <- mydata_1[,-1] # on supprime la première colonne qui ne nous interesse plus

# creation de la heatmap 
heatmaply(mydata_1, scale_fill_gradient_fun = gradient_col, 
          Colv=NA, Rowv=NA, scale='none', limits = c(-7, 8))




## Chromosome 11 ##
chr11 <- subset(gene.exp, gene.exp$chr == 11 )
mydata_11 <- chr11[,3:18] 
rownames(mydata_11) <- mydata_11[,1]
mydata_11 <- mydata_11[,-1]

heatmaply(mydata_11, scale_fill_gradient_fun = gradient_col, 
          Colv=NA, Rowv=NA, scale='none', limits = c(-4, 8))



## Chromosome 17 ##
chr17 <- subset(gene.exp, gene.exp$chr == 17 )
mydata_17 <- chr17[,3:18] 
rownames(mydata_17) <- mydata_17[,1]
mydata_17 <- mydata_17[,-1]

heatmaply(mydata_17, scale_fill_gradient_fun = gradient_col, 
          Colv=NA, Rowv=NA, scale='none', limits = c(-7, 8))



