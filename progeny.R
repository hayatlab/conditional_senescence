## We load the required packages
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(readr)


setwd("~/Desktop/progeny/xioxong/DESEQ2.EC_AA.versus.EC_AB/")
counts <- read_csv("~/Desktop/progeny/xioxong/DESEQ2.EC_AA.versus.EC_AB/post_TMM.csv")
colnames(counts)[1] <- "gene"
# Remove NAs and set row names
counts <- counts %>%
  dplyr::mutate_if(~ any(is.na(.x)), ~ if_else(is.na(.x),0,.x)) %>% 
  column_to_rownames(var = "gene") %>% 
  as.matrix()
head(counts)
counts <- counts[,-c(1,2,3,7,8,9)]
counts <- counts[, c(4, 5,6, 1, 2, 3)]

design <- read.csv("~/Desktop/progeny/xioxong/pdgfrb/pdgfrb_targets.csv",sep = ";",row.names = NULL)
design
design <- design[,-c(1,2,3,7,8,9)]
design
## We read the results from the differential analysis. 
ttop_KOvsWT <- read.csv("~/Desktop/progeny/xioxong/pdgfrb/DESEQ2.doxN.pdgfrb.versus.doxP.pdgfrb.csv",sep = ",",row.names = 1)
colnames(ttop_KOvsWT)[1] <- "ID"
colnames(ttop_KOvsWT)[5] <- "t"
# Extract t-values per gene
deg <- ttop_KOvsWT %>%
  dplyr::select(ID, t) %>% 
  filter(!is.na(t)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()
head(deg)
head(deg)
net <- get_progeny(organism = 'human', top = 100)
net
# Run wmean
sample_acts <- run_wmean(mat=counts, net=net, .source='source', .target='target',
                         .mor='weight', times = 100, minsize = 5)
sample_acts
# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()
# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)
sample_acts_mat=t(sample_acts_mat)
sample_acts_mat <- sample_acts_mat[, c(4, 5,6, 1, 2, 3)]
sample_acts_mat <- sample_acts_mat[!grepl("^VEGF", rownames(sample_acts_mat)), ]
# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

jpeg(filename = paste0("progeny_results_pdgfrb.jpeg"),  width=1500,height =1500,quality = 100,res = 300)

# Plot
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks, cluster_cols = F) 
dev.off()
pdf("progeny_results_pdgfrb.pdf", width = 8, height = 8)

# Plot
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks, cluster_cols = F) 
dev.off()
