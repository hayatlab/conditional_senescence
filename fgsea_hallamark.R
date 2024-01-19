library(fgsea)
library(tidyverse)


res <- read.table("DESEQ2.doxN.cd10.versus.doxP.cd10.txt",header = TRUE, sep = "")
head(res)
colnames(res)[1]  <- "SYMBOL"    # change column name for x column


### Further, all you’ll care about later on is the gene symbol and the test statistic.####
##Get just those, and remove the NAs. Finally, if you have multiple test statistics for###
####the same symbol, you’ll want to deal with that in some way. Here I’m just averaging them.####
res2 <- res %>%
  dplyr::select(SYMBOL, stat) %>%
  na.omit() %>%
  filter(!is.na(stat) & is.numeric(stat)) %>%  # Filter out non-numeric values
  distinct() %>%
  group_by(SYMBOL) %>%
  summarize(stat = mean(stat, na.rm = TRUE))
res2
ranks <- deframe(res2)
head(ranks, 20)



# Load the pathways into a named list
pathways.hallmark <- gmtPathways("~/Desktop/projects/fgsea/h.all.v6.2.symbols.gmt")


# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

###### Now, run the fgsea algorithm with 1000 permutations:####
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
######## Representation#####
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
fgseaResTidy$pathway=str_replace(fgseaResTidy$pathway, "HALLMARK_", "")
p=ggplot(fgseaResTidy %>%filter(padj <0.05) ,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
pdf("Hallmark pathways NES.pdf", width=10,height=15)
print(p)
dev.off()

readr::write_csv(fgseaResTidy, "fgsea_hallmark_.csv")

