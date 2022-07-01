C <- as.matrix(read.table("pcangsd_covmatrix_noout.cov"))
ids <- read.table("pca_pops_noout.txt") #text file with 20 lines of the single word BOD, then 20 lines of CAP etc in the order they appeared in all_rmdups_jo.txt
e <- eigen(C)
# base R
plot(e$vectors[,1:2],xlab="PC1",ylab="PC2", bg=ids$V1, pch=21)
#ggplot
library(ggplot2)
library(tidyverse)
df <- data.frame(Population = ids$V1, PC1 = e$vectors[,1], PC2 = e$vectors[,2])
df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=Population)) +
  geom_point(size=5, shape=21) +
  theme_bw() +
  labs(y = "PC2 (1.03%)", x = "PC1 (1.4%)") +
  scale_fill_brewer(palette="Accent",breaks=c('FOG', 'CAP', 'KIB','BOD','TER','LOM','SAN'))

#scale_fill_discrete() +