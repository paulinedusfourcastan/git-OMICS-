# Chargement des librairies 
library(OmicCircos)
library(circlize) 
library(ComplexHeatmap)


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

gene.exp <- gene.exp[, -2] # on supprime la colonne position qui est inutile dans la suite



### HEATMAP SUR TOUS LES CHROMOSOMES ###

# Mettre les donnees dans le bon format :

mat <- as.matrix(gene.exp, labels=TRUE) # on transforme le dataframe en matrice 

# On définit les noms des colonnes et des lignes 
row_names <- as.vector(gene.exp$NAME)
col_names <- c("chr", "NAME", "TCGA.A2.A04W.01A", "TCGA.A2.A0D1.01A",
               "TCGA.A2.A0EQ.01A", "TCGA.A8.A081.01A", "TCGA.A8.A08J.01A",
               "TCGA.A8.A09X.01A", "TCGA.AN.A0FV.01A", "TCGA.AO.A0J2.01A",
               "TCGA.AO.A12D.01A", "TCGA.B6.A0RS.01A", "TCGA.BH.A0EE.01A", 
               "TCGA.C8.A12L.01A", "TCGA.C8.A12P.01A", "TCGA.D8.A13Z.01A", 
               "TCGA.E2.A1B0.01A")
rownames_factor <- factor(row_names)

# On met notre matrice sous format "numeric"
mat_numeric <- matrix(as.numeric(mat), ncol = ncol(mat))
rownames(mat_numeric) = row_names
colnames(mat_numeric) = col_names

# Vecteur des chromosomes :
split <- as.vector(gene.exp$chr)
split <- factor(split)



## Affichage des heatmaps
circos.clear()
circos.par(start.degree = 90, gap.degree = 3)
col_fun1 = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))

# heatmap simple 
circos.heatmap(mat_numeric[, !(colnames(mat_numeric) %in% c("chr", "NAME"))], 
               split = split, 
               col = col_fun1)


# heatmap avec les chromosomes associes
circos.clear()
circos.heatmap(mat_numeric[, !(colnames(mat_numeric) %in% c("chr", "NAME"))], 
               split = split, col = col_fun1, track.height = 0.4, 
               bg.border = "green", bg.lwd = 2, bg.lty = 2, 
               show.sector.labels = TRUE)


# heatmap avec les chromosomes et les noms des genes
circos.clear()
circos.heatmap(mat_numeric[, !(colnames(mat_numeric) %in% c("chr", "NAME"))], 
               split = split, col = col_fun1, track.height = 0.4, 
               bg.border = "green", bg.lwd = 2, bg.lty = 2, 
               show.sector.labels = TRUE, rownames.side = "outside")




# heatmap avec les dendrogrammes 
circos.clear()
circos.heatmap(mat_numeric[, !(colnames(mat_numeric) %in% c("chr", "NAME"))], 
               split = split, col = col_fun1, track.height = 0.4, 
               bg.border = "green", bg.lwd = 2, bg.lty = 2, 
               show.sector.labels = TRUE, dend.side = "outside")
lgd = Legend(title = "mat", col_fun = col_fun1)
grid.draw(lgd)






### HEATMAP SUR UNE SELECTION DE CHROMOSOMES ###

x <- c(1, 3, 6, 11, 17, 19) # vecteur des chromosomes que l'on souhaite selectionner
data_chr <- gene.exp[gene.exp$chr %in% x, ] # dataframe avec nos chromosomes selectionnes


# Convertion en matrice 
mat_chr <- as.matrix(data_chr, labels=TRUE)
mat_chr_numeric <- matrix(as.numeric(mat_chr), ncol = ncol(mat_chr))

# On définit les noms des colonnes et des lignes 
row_names_chr <- as.vector(data_chr$NAME)
col_names_chr <- c("chr", "NAME", "TCGA.A2.A04W.01A", "TCGA.A2.A0D1.01A",
               "TCGA.A2.A0EQ.01A", "TCGA.A8.A081.01A", "TCGA.A8.A08J.01A",
               "TCGA.A8.A09X.01A", "TCGA.AN.A0FV.01A", "TCGA.AO.A0J2.01A",
               "TCGA.AO.A12D.01A", "TCGA.B6.A0RS.01A", "TCGA.BH.A0EE.01A", 
               "TCGA.C8.A12L.01A", "TCGA.C8.A12P.01A", "TCGA.D8.A13Z.01A", 
               "TCGA.E2.A1B0.01A")


rownames(mat_chr_numeric) = row_names_chr
colnames(mat_chr_numeric) = col_names_chr

# Vecteur des chromosomes :
split2 <- as.vector(data_chr$chr)
split2 <- factor(split2)




# heatmap simple 
circos.clear()
circos.heatmap(mat_chr_numeric[, !(colnames(mat_chr_numeric) %in% c("chr", "NAME"))], 
               split = split2, 
               col = col_fun1)


# heatmap avec les chromosomes associes
circos.clear()
circos.heatmap(mat_chr_numeric[, !(colnames(mat_chr_numeric) %in% c("chr", "NAME"))], 
               split = split2, col = col_fun1, track.height = 0.4, 
               bg.border = "green", bg.lwd = 2, bg.lty = 2, 
               show.sector.labels = TRUE)


# heatmap avec les chromosomes et les noms des genes
circos.clear()
circos.heatmap(mat_chr_numeric[, !(colnames(mat_chr_numeric) %in% c("chr", "NAME"))], 
               split = split2, col = col_fun1, track.height = 0.4, 
               bg.border = "green", bg.lwd = 2, bg.lty = 2, 
               show.sector.labels = TRUE, rownames.side = "outside")




# heatmap avec les dendrogrammes 
circos.clear()
circos.heatmap(mat_chr_numeric[, !(colnames(mat_chr_numeric) %in% c("chr", "NAME"))], 
               split = split2, col = col_fun1, track.height = 0.4, 
               bg.border = "green", bg.lwd = 2, bg.lty = 2, 
               show.sector.labels = TRUE, dend.side = "outside")
lgd = Legend(title = "mat", col_fun = col_fun1)
grid.draw(lgd)





