
# written by Hyojin Kim
# ------------------------- # 

# Strategy 
# ------------------------- # 
# counts (Normalize) -> cpm (Normalize) -> TMM (inter-sample normalize) 


library(ggplot2)
library(reshape)
library(dplyr)
library(cowplot)
library(ggrepel)
library(DESeq2)
library(textshape)
library(HTSFilter)
library(tidyverse)
library(tibble)
library(edgeR)
library(glmpca)
library(readr)
library(data.table)
library(progeny)
library(dorothea)
library(readr)
library(ggrepel)
library(matrixStats)
library(plotly)
library(stats)
library(htmlwidgets)
library(viridis)
library(RColorBrewer)
library(patchwork)

'%ni%' = Negate('%in%')
args <- commandArgs(trailingOnly=TRUE)
set.seed(42)


## ---------------------------------- ##
args[1]        <- type_your_folder_name 
args[2]        <- type_your_species           ### "human" or "mouse" 
outdir_common  <- type_your_output_directory
indir          <- paste0(outdir_common, args[1], "/star_salmon/")
outdir         <- paste0(outdir_common, args[1], "/star_salmon/preprocessed/") 
if(!dir.exists(outdir)) dir.create(outdir);


## ---------------------------------- ##
D <- "salmon.merged.gene_counts.tsv" 
counts_df <- as.data.frame(read_delim(paste0(indir, D), 
				"\t", 
				escape_double = FALSE, 
				trim_ws = TRUE)) 


## ---------------------------------- ##
## filter out rRNA, mtRNA 
## ---------------------------------- ##

if (args[2] == "human") {
	clear_G_db <- as.data.frame(read_delim("/XXX/genecode_db/gencode.v38.annotation.gtf.no.rRNA.tRNA.gtf",
				'gene_name "',
                                col_names = FALSE,
                                escape_double = FALSE,
                                trim_ws = TRUE, 
                                ))   
} else { 
	clear_G_db <- as.data.frame(read_delim("/XXX/gencode.vM27.annotation.gtf.no.rRNA.tRNA.gtf",
				'gene_name "',
				col_names = FALSE,
				escape_double = FALSE,
				trim_ws = TRUE, 
				))  
	}
print ("it'll take some time")


clear_G <- strsplit(clear_G_db$X2, '["\"]' ) %>% lapply(., function(x) x[1]) %>% do.call("rbind", .) %>% unique()
rm(clear_G_db)

counts_df <- counts_df %>%
		dplyr::select(-gene_id) %>%  
		filter(gene_name %in% clear_G)  


counts_df_orig <- counts_df
counts_df_colnames <- colnames(counts_df[2:length(colnames(counts_df))])



## ---------------------------------- ##
## get maximum gene exp for the identical gene id 
## ---------------------------------- ##

order_G <- counts_df$gene_name %>% unique()
max_col <- length(colnames(counts_df)) 
columns <- colnames(counts_df)[2:max_col]
counts_df_max <- c()
for ( column in columns ) {
	counts_df_max[[column]] <- counts_df %>% 
				dplyr::select(gene_name, column) %>% 
				group_by(gene_name) %>% 
				mutate(max = max(eval(parse(text=column)))) %>% 
				dplyr::select(-column) %>%
				unique() %>% 
				as.data.frame() %>%  
				mutate(gene_name = fct_reorder(gene_name, order_G)) 
	
	rownames(counts_df_max[[column]]) <- counts_df_max[[column]]$gene_name
	counts_df_max[[column]] <- counts_df_max[[column]] %>% dplyr::select(-gene_name)
	colnames(counts_df_max[[column]]) <- column

	}


T <- do.call("cbind", counts_df_max) 



## ---------------------------------- ##
## make density plot before normalization
## ---------------------------------- ## 
pdf(paste0(outdir, "density.pre.pdf"),width = 5, height = 5)
d <- density(log2(as.matrix(T)+1)) # returns the density data 
print(plot(d)) # plots the result
dev.off()


## ---------------------------------- ##
## prepare meta 
## ---------------------------------- ##
colnames(T) <- gsub("__", "_", colnames(T))
condition <- gsub( paste0("_REP","[0-9]+"),"", colnames(T))
targets <- data.frame(sample=colnames(T), condition=condition)


## ---------------------------------- ##
## run HTSFilter to remove lowly expressed genes 
## ---------------------------------- ##
if (args[3] == "REP") {
	T <- HTSFilter(as.matrix(T), condition, s.min = 0.1, normalization ="TMM", s.len=25, plot=TRUE)
	T <- T$filteredData
} else {
	cutoff <- quantile(as.matrix(T), probs = c(0.00))[[1]]
	print ("------------------------------")
	print ("print 5% quantile for cutoff")
	print (cutoff)
	print ("------------------------------")
	keep <- rowSums(as.matrix(T)) > cutoff
	T <- T[keep,]
}


## ---------------------------------- ##
## make density plot after normalization
## ---------------------------------- ##
pdf(paste0(outdir, "density.post.pdf"),width = 5, height = 5)
d <- density(log2(as.matrix(T)+1)) # returns the density data
print(plot(d)) # plots the result
dev.off()



## ---------------------------------- ##
## add data to edgeR object
## ---------------------------------- ##
T <- DGEList(T, group=condition)
T <- calcNormFactors(T, method="TMM")

plotMD(T, column=1)
abline(h=0, col="red", lty=2, lwd=2)



## ---------------------------------- ##
## box plot 
## ---------------------------------- ##
pre_TMM <- cpm(T, normalized.lib.sizes=FALSE)
pdf(paste0(outdir, "density.boxplot.wo.w.TMM.pdf"),width = 5, height = 5)
# Basic box plot
p1 <- log2(pre_TMM+1) %>% 
	melt %>% 
	ggplot(aes(x=Var2, y=value)) + 
	geom_boxplot()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=4, hjust=1))+
	xlab("")

post_TMM <- cpm(T, normalized.lib.sizes=TRUE)
# Basic box plot
p2 <- log2(post_TMM+1) %>% 
        melt %>% 
        ggplot(aes(x=Var2, y=value)) +             
        geom_boxplot() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=4, hjust=1))+
	xlab("")

print(plot_grid(p1, p2))
dev.off()

write.table(post_TMM, file=paste0(outdir, "post_TMM.txt"),row.names = TRUE, quote = FALSE, sep="\t")




## ---------------------------------- ##
## Save filered count matrix          ## 
## ---------------------------------- ##
T$counts %>% write.table(file=paste0(outdir, "filtered.count.txt"),
                         row.names = TRUE, 
			 col.names = TRUE,
                         quote = FALSE,
                         sep="\t")





## ---------------------------------- ##
## PCA 1 :: Prcomp	              ##
## ---------------------------------- ##
sample <- targets$sample
X <- t(post_TMM)
pc <- prcomp(X)
explained_vairiance_ratio <- summary(pc)[["importance"]]['Proportion of Variance',]
percentVar <- explained_vairiance_ratio
df <- cbind(pc$x[,c(1,2,3)], condition, sample) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(X)))
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(X)))
df$PC3 <- as.numeric(df$PC3) / (pc$sdev[3] * sqrt(nrow(X)))

## ---------------------------------- ##
pdf(paste0(outdir, "pca.2d.prcomp.condition.detail.pdf"),width = 6, height = 6)
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = sample)) +
        geom_point(size =3) +
        xlab(paste0("PC1: ", percentVar[1]*100, "% variance")) +
        ylab(paste0("PC2: ", percentVar[2]*100, "% variance")) +
        coord_fixed() +
        theme_light() +
        ggtitle("PCA with inter & intra normalized data")
print(p1)
dev.off()

## ---------------------------------- ##
pdf(paste0(outdir, "pca.2d.prcomp.condition.pdf"),width = 6, height = 10)
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = condition)) +
        geom_point(size =3) +
        xlab(paste0("PC1: ", percentVar[1]*100, "% variance")) +
        ylab(paste0("PC2: ", percentVar[2]*100, "% variance")) +
        coord_fixed() +
        theme_light() +
        ggtitle("PCA with inter & intra normalized data")
p2 <- ggplot(df, aes(x = PC1, y = PC3, color = condition)) +
        geom_point(size =3) +
        xlab(paste0("PC1: ", percentVar[1]*100, "% variance")) +
        ylab(paste0("PC3: ", percentVar[3]*100, "% variance")) +
        coord_fixed() +
        theme_light() 
p3 <- ggplot(df, aes(x = PC2, y = PC3, color = condition)) +
        geom_point(size =3) +
        xlab(paste0("PC2: ", percentVar[2]*100, "% variance")) +
        ylab(paste0("PC3: ", percentVar[3]*100, "% variance")) +
        coord_fixed() +
        theme_light() 
p1 + p2 + p3 + plot_layout(ncol = 1)
dev.off()




## ---------------------------------- ##
## PCA 2 :: DESeq2
## ---------------------------------- ##
data = "filtered.count.txt"
dat = read.table(paste0(outdir, "/", data))
design = colnames(dat) %>% strsplit(., "_REP") %>% lapply(., function(x) x[1]) %>% unlist() %>% factor()
sample <- colnames(dat) %>% factor()
meta = data.frame(sample, design)
colnames(meta) = c("sample", "condition")

## ---------------------------------- ##
dds <- DESeqDataSetFromMatrix(countData = round(dat,0), colData = meta, design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
resultsNames(dds)

## ---------------------------------- ##
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c( "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

## ---------------------------------- ##
pdf(paste0(outdir, "pca.2d.plotPCA.condition.pdf"),width = 5, height = 8)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
        geom_point(size =3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed() +
        theme_light() +
        ggtitle("PCA with VST data")
dev.off()

## ---------------------------------- ##
pdf(paste0(outdir, "pca.2d.plotPCA.sample.pdf"),width = 6, height = 8)
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
        geom_point(size =3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed() +
        theme_light() +
        ggtitle("PCA with VST data")
dev.off()






