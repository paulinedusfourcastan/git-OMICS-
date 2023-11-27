library(OmicCircos);
data("TCGA.PAM50_genefu_hg18"); ## vrai jeu de données
# se renseigner sur les données
data("TCGA.BC.fus"); # nombre de fusion des protéines
data("TCGA.BC.cnv.2k.60"); # quand y a un cancer, le nombre de genes copies
data("TCGA.BC.gene.exp.2k.60");
data("TCGA.BC.sample60");
data("TCGA.BC_Her2_cnv_exp"); # etude de la protéine her2 

# on redimensionne les p-values 
pvalue = -1 * log10(TCGA.BC_Her2_cnv_exp[,5]);

# on crée un tableau avec les nouvelles p-value et les 3 premières colonnes du dataframe
# chr - po - gene - p-value
pvalue = cbind(TCGA.BC_Her2_cnv_exp[,c(1:3)], pvalue);

# numéro du tableau de type de cancer Her2
Her2.i = which(TCGA.BC.sample60[,2] == "Her2");

# appellation des cancers Her2
Her2.n = TCGA.BC.sample60[Her2.i,1];

# 
Her2.j = which(colnames(TCGA.BC.cnv.2k.60)%in%Her2.n)


# data avec juste les Her.2
cnv = TCGA.BC.cnv.2k.60[,c(1:3,Her2.j)];

cnv.m = cnv[,c(4:ncol(cnv))];
cnv.m[cnv.m > 2] = 2;
cnv.m[cnv.m < -2] = -2;
cnv = cbind(cnv[,1:3], cnv.m)

Her2.j = which(colnames(TCGA.BC.gene.exp.2k.60)%in% Her2.n);

gene.exp = TCGA.BC.gene.exp.2k.60[,c(1:3,Her2.j)]


#circos(R = 300, cir = "hg18", W = 100, mapping = gene.exp, 
       #col.v = 4, type = "heatmap2", cluster = FALSE, col.bar = TRUE, 
       #lwd = 0.01, zoom = zoom);


library(pheatmap)
library(plotly)
library(heatmaply)
library(ggplot2)

gradient_col <- ggplot2::scale_fill_gradient2(
  low = "blue", high = "red", 
  midpoint = 0, limits = c(-7, 8)
)

# Chromosome 1 
chr1 <- subset(gene.exp, gene.exp$chr == 1 )
mydata_1 <- chr1[,3:18] 
rownames(mydata_1) <- mydata_1[,1]
mydata_1 <- mydata_1[,-1]

heatmaply(mydata_1, scale_fill_gradient_fun = gradient_col, 
          Colv=NA, Rowv=NA, scale='none', limits = c(-7, 8))


# Chromose 11
chr11 <- subset(gene.exp, gene.exp$chr == 11 )
mydata_11 <- chr11[,3:18] 
rownames(mydata_11) <- mydata_11[,1]
mydata_11 <- mydata_11[,-1]

heatmaply(mydata_11, scale_fill_gradient_fun = gradient_col, 
          Colv=NA, Rowv=NA, scale='none', limits = c(-4, 8))

#Chromose 17 
chr17 <- subset(gene.exp, gene.exp$chr == 17 )
mydata_17 <- chr17[,3:18] 
rownames(mydata_17) <- mydata_17[,1]
mydata_17 <- mydata_17[,-1]

heatmaply(mydata_17, scale_fill_gradient_fun = gradient_col, 
          Colv=NA, Rowv=NA, scale='none', limits = c(-7, 8))
