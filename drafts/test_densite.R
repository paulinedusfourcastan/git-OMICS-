library(OmicCircos)
library(kdensity)
library(ggplot2)

data("TCGA.PAM50_genefu_hg18")

## Etude densite

# Chromosome 1 

chr1 <- TCGA.PAM50_genefu_hg18[1:7,]

ggplot(chr1, aes(..scaled..)) +
  geom_density(aes(x = Her2, y=..scaled..,colour = "Her2"), show.legend = TRUE) +
  geom_density(aes(x = Basal, y=..scaled.., colour = "Basal"), show.legend = TRUE) +
  geom_density(aes(x = LumA, y=..scaled.., colour = "LumA"), show.legend = TRUE) +
  geom_density(aes(x = LumB, y=..scaled.., colour = "LumB"), show.legend = TRUE) +
  geom_density(aes(x = Normal, y=..scaled..,colour = "Normal"), show.legend = TRUE) +
  labs(color = "Variable") +
  theme_minimal()+
  ggtitle("Densité chromosome 1") +
  scale_y_continuous(labels = scales::percent_format(scale = 1))


# Chromosome 17 

chr17 <- TCGA.PAM50_genefu_hg18[35:41,]

ggplot(chr17, aes(..scaled..)) +
  geom_density(aes(x = Her2, y=..scaled..,colour = "Her2"), show.legend = TRUE) +
  geom_density(aes(x = Basal, y=..scaled.., colour = "Basal"), show.legend = TRUE) +
  geom_density(aes(x = LumA, y=..scaled.., colour = "LumA"), show.legend = TRUE) +
  geom_density(aes(x = LumB, y=..scaled.., colour = "LumB"), show.legend = TRUE) +
  geom_density(aes(x = Normal, y=..scaled..,colour = "Normal"), show.legend = TRUE) +
  labs(color = "Variable") +
  theme_minimal()+
  ggtitle("Densité chromosome 1") +
  scale_y_continuous(labels = scales::percent_format(scale = 1))



# Tous les chromosomes 

ggplot(TCGA.PAM50_genefu_hg18, aes(..scaled..)) +
  geom_density(aes(x = Her2, y=..scaled..,colour = "Her2"), show.legend = TRUE) +
  geom_density(aes(x = Basal, y=..scaled.., colour = "Basal"), show.legend = TRUE) +
  geom_density(aes(x = LumA, y=..scaled.., colour = "LumA"), show.legend = TRUE) +
  geom_density(aes(x = LumB, y=..scaled.., colour = "LumB"), show.legend = TRUE) +
  geom_density(aes(x = Normal, y=..scaled..,colour = "Normal"), show.legend = TRUE) +
  labs(color = "Variable") +
  theme_minimal()+
  ggtitle("Densité tous les chromosomes") +
  scale_y_continuous(labels = scales::percent_format(scale = 1))
