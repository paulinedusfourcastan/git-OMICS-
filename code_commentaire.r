rm(list = ls())
options(stringsAsFactors = FALSE)
library(OmicCircos)

data("TCGA.PAM50_genefu_hg18")
data("TCGA.BC.fus")
data("TCGA.BC.cnv.2k.60")
data("TCGA.BC.gene.exp.2k.60")
data("TCGA.BC.sample60")
data("TCGA.BC_Her2_cnv_exp")

pvalue <- -1 * log10(TCGA.BC_Her2_cnv_exp[,5])
pvalue <- cbind(TCGA.BC_Her2_cnv_exp[,c(1:3)], pvalue)

Her2.i <- which(TCGA.BC.sample60[,2] == "Her2")
Her2.n <- TCGA.BC.sample60[Her2.i,1]

Her2.j <- which(colnames(TCGA.BC.cnv.2k.60) %in% Her2.n)
cnv <- TCGA.BC.cnv.2k.60[,c(1:3,Her2.j)]

cnv.m <- cnv[,c(4:ncol(cnv))]
cnv.m[cnv.m > 2] <- 2
cnv.m[cnv.m < -2] <- -2
cnv <- cbind(cnv[,1:3], cnv.m)

Her2.j <- which(colnames(TCGA.BC.gene.exp.2k.60) %in% Her2.n)
gene.exp <- TCGA.BC.gene.exp.2k.60[,c(1:3,Her2.j)]

colors <- rainbow(10, alpha = 0.5)

pdf("OmicCircos4vignette10_paper1.pdf", 8, 8)
par(mar = c(2, 2, 2, 2))
plot(c(1, 800), c(1, 800), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")

zoom <- c(1, 22, 939245.5, 154143883, 0, 180)
#  Mise  des chromosomes 
circos(R = 400, cir = "hg18", W = 4, type = "chr", print.chr.lab = TRUE, scale = TRUE, zoom = zoom)
# génes 
circos(R = 300, cir = "hg18", W = 100, mapping = gene.exp, col.v = 4, type = "heatmap2", cluster = FALSE, col.bar = TRUE, lwd = 0.01, zoom = zoom)

# 
circos(R = 220, cir = "hg18", W = 80, mapping = cnv, col.v = 4, type = "ml3", B = FALSE, lwd = 1, cutoff = 0, zoom = zoom)

circos(R = 140, cir = "hg18", W = 80, mapping = pvalue, col.v = 4, type = "l", B = TRUE, lwd = 1, col = colors[1], zoom = zoom)



circos(R = 130, cir = "hg18", W = 10, mapping = TCGA.BC.fus, type = "link", lwd = 2, zoom = zoom)

# Zoom in links using the highlight functions
the.col1 <- rainbow(10, alpha = 0.5)[1]
highlight <- c(140, 400, 11, 282412.5, 11, 133770314.5, the.col1, the.col1)
circos(R = 110, cir = "hg18", W = 40, mapping = highlight, type = "hl", lwd = 2, zoom = zoom)

the.col2 <- rainbow(10, alpha = 0.5)[6]
highlight <- c(140, 400, 17, 739525, 17, 78385909, the.col2, the.col2)
circos(R = 110, cir = "hg18", W = 40, mapping = highlight, type = "hl", lwd = 2, zoom = zoom)

# Highlight link
highlight.link1 <- c(400, 400, 140, 376.8544, 384.0021, 450, 540.5)
circos(cir = "hg18", mapping = highlight.link1, type = "highlight.link", col = the.col1, lwd = 1)

highlight.link2 <- c(400, 400, 140, 419.1154, 423.3032, 543, 627)
circos(cir = "hg18", mapping = highlight.link2, type = "highlight.link", col = the.col2, lwd = 1)

# Zoom in chromosome 11
zoom <- c(11, 11, 282412.5, 133770314.5, 180, 270)
circos(R = 400, cir = "hg18", W = 4, type = "chr", print.chr.lab = TRUE, scale = TRUE, zoom = zoom)
circos(R = 300, cir = "hg18", W = 100, mapping = gene.exp, col.v = 4, type = "heatmap2", cluster = FALSE, lwd = 0.01, zoom = zoom)
circos(R = 220, cir = "hg18", W = 80, mapping = cnv, col.v = 4, type = "ml3", B = FALSE, lwd = 1, cutoff = 0, zoom = zoom)
circos(R = 140, cir = "hg18", W = 80, mapping = pvalue, col.v = 4, type = "l", B = TRUE, lwd = 1, col = colors[1], zoom = zoom)

# Gene label
gene.names <- c("ERBB2", "CDC6")
PAM50.17 <- which(TCGA.PAM50_genefu_hg18[,3] == gene.names)
TCGA.PAM50 <- TCGA.PAM50_genefu_hg18[PAM50.17,]
TCGA.PAM50 <- rbind(TCGA.PAM50, TCGA.PAM50[2,])
TCGA.PAM50[3,1] <- 11
TCGA.PAM50[3,2] <- 69165000
TCGA.PAM50[,3] <- as.character(TCGA.PAM50[,3])
TCGA.PAM50[3,3] <- c("CCND1")
circos(R = 410, cir = "hg18", W = 40, mapping = TCGA.PAM50, type = "label", side = "out", col = "blue", zoom = zoom)

# Zoom in chromosome 17
zoom <- c(17, 17, 739525, 78385909, 274, 356)
circos(R = 400, cir = "hg18", W = 4, type = "chr", print.chr.lab = TRUE, scale = TRUE, zoom = zoom)
circos(R = 300, cir = "hg18", W = 100, mapping = gene.exp, col.v = 4, type = "heatmap2", cluster = FALSE, lwd = 0.01, zoom = zoom)
circos(R = 220, cir = "hg18", W = 80, mapping = cnv, col.v = 4, type = "ml3", B = FALSE, lwd = 1, cutoff = 0, zoom = zoom)
circos(R = 140, cir = "hg18", W = 80, mapping = pvalue, col.v = 4, type = "l", B = TRUE, lwd = 1, col = colors[1], zoom = zoom)

# Gene label
gene.names <- c("ERBB2", "CDC6")
PAM50.17 <- which(TCGA.PAM50_genefu_hg18[,3] == gene.names)
TCGA.PAM50 <- TCGA.PAM50_genefu_hg18[PAM50.17,]
TCGA.PAM50 <- rbind(TCGA.PAM50, TCGA.PAM50[2,])
TCGA.PAM50[3,1] <- 11
TCGA.PAM50[3,2] <- 69165000
TCGA.PAM50[,3] <- as.character(TCGA.PAM50[,3])
TCGA.PAM50[3,3] <- c("CCND1")
circos(R = 410, cir = "hg18", W = 40, mapping = TCGA.PAM50, type = "label", side = "out", col = "blue", zoom = zoom)

dev.off()
