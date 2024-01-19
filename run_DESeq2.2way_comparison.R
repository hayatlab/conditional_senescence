

library(DESeq2)
library(dplyr)
library(tibble)
library(gplots)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(progeny)
library(dorothea)
args <- commandArgs(trailingOnly=TRUE)

project = type_your_project_name 
dir = paste0("/XXX/", project, "/star_salmon/preprocessed/")
data = "filtered.count.txt"
norm = "post_TMM.txt"
dat = read.table(paste0(dir, "/", data))
nom = read.table(paste0(dir, "/", norm))
design = colnames(dat) %>% strsplit(., "_REP") %>% lapply(., function(x) x[1]) %>% unlist() %>% factor()
sample <- colnames(dat) %>% factor()
meta = data.frame(sample, design)
colnames(meta) = c("sample", "condition")
outdir = dir

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
pdf(paste0(outdir, "pca.2d.plotPCA.sample.pdf"),width = 6, height = 6)
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
        geom_point(size =3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed() +
        theme_light() +
        ggtitle("PCA with VST data")
dev.off()




get_top20DE <- function(E_up, E_up_gene_df, N) {
        E_up$log2FoldChange.avg <- E_up %>% dplyr::select(log2FoldChange.x,  log2FoldChange.y) %>% rowMeans(.)
        E_up_gene <- E_up[ E_up$rowname %in% E_up_gene_df$rowname, ] %>% arrange(., desc(log2FoldChange.avg)) %>% top_n(N) %>% pull(rowname)
        return (E_up_gene)
        }


# ge DE
get_DEgene <- function(dds,A,B,folder,qval_cutoff, log2_cutoff) {
	# ----------------------------------- #
	# condition column name should be "condition"
	# ----------------------------------- #
	A_B <- c("condition", A, B)
	# A vs B
	res_table_A_B <- results(dds, contrast=A_B, name=paste0(A, ".vs.", B), alpha = qval_cutoff)
	res_table_A_B["index"] = paste0(A, ".vs.", B)
	res_table_A_B %>% as.data.frame() -> res_table_A_B
	res_table_A_B[!is.na(res_table_A_B$padj),] -> res_table_A_B
	res_table_A_B_df = res_table_A_B %>% as.data.frame() %>% rownames_to_column(.)
	#
	res = res_table_A_B
	x_max = res$log2FoldChange %>% max()
	x_min = res$log2FoldChange %>% min() %>% abs()
	x_lim_max <- max(x_max, x_min)
	pdf(paste0(outdir, folder, "/",  "DESeq2.volcano.log2FC", log2_cutoff, ".qval.", qval_cutoff, ".pdf"),width=7, height=7)
	p <- EnhancedVolcano(res, lab = rownames(res),
		x = 'log2FoldChange',
		y = 'padj',
		xlim = c(x_lim_max*-1, x_lim_max),
		pointSize = 3.0,
		labSize = 1.0,
		ylab = bquote(~-Log[10]~adjusted~italic(P)),
		#drawConnectors = TRUE,
		#widthConnectors = 0.5,
		#colAlpha = 1,
		#maxoverlapsConnectors = Inf,
		FCcutoff = log2_cutoff,
		pCutoff = qval_cutoff
		)
	print (p)
	dev.off()
	return (res_table_A_B_df)	
	}


folders     = c( "DESEQ2.doxN.cd10.versus.doxP.cd10", "DESEQ2.doxN.pdgfrb.versus.doxP.pdgfrb" )
AS          = c( "doxN.cd10", "doxN.pdgfrb")
BS          = c( "doxP.cd10", "doxP.pdgfrb")
qval_cutoff = 0.001
log2_cutoff = 2.0


N = length(folders)
for ( i in seq(1,N)) {
	new_dir = paste0(outdir, folders[i], "/")
	if(!dir.exists(new_dir)) dir.create(new_dir);	
	A_up <- get_DEgene(dds, AS[i], BS[i], folders[i], qval_cutoff, log2_cutoff)
	A_up %>% write.table(., file=paste0(outdir, folders[i], "/", 
				"DESEQ2.", AS[i], ".versus.", BS[i], ".txt"), 
				sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
	}

nom_list <- list(nom)
names(nom_list) <- c("all")

##PROGENY
run_progeny <- function(nom, nom_name) {
	nom <- as.matrix(nom)
	PathwayActivity_counts <- progeny(nom, scale=TRUE, organism="Human", top = 100)
	Activity_counts <- as.vector(PathwayActivity_counts)
	
	paletteLength <- 100
	myColor <- colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)
	
	progenyBreaks <- c(seq(min(Activity_counts), 0,
                        length.out=ceiling(paletteLength/2) + 1),
                        seq(max(Activity_counts)/paletteLength,
                        max(Activity_counts),
                        length.out=floor(paletteLength/2)))


	pdf(paste0(outdir, "PROGENy.", nom_name, ".pdf"), width = 7, height = 6)
	print(pheatmap(t(PathwayActivity_counts),fontsize=14,
			cluster_rows=FALSE, cluster_cols=FALSE,
                        fontsize_row = 10, fontsize_col = 10,
                        color=myColor, breaks = progenyBreaks,
                        main = "PROGENy (100)", angle_col = 90,
                        treeheight_col = 0,  border_color = NA)
			)
	dev.off()
	}


run_progeny(nom_list[[1]], "all")

## DOROTHEA
run_dorothea <- function(nom, nom_name) {
	data(dorothea_hs, package = "dorothea")
	regulons <- dorothea_hs %>%
                        dplyr::filter(confidence %in% c("A", "B","C"))
	
	tf_activities_counts <-
                dorothea::run_viper(nom, regulons,
                options =  list(minsize = 5, eset.filter = FALSE,
                cores = 1, verbose = FALSE, method = c("scale")))

	var_df <- rowVars(tf_activities_counts)
	tf_activities_counts <- as.data.frame(tf_activities_counts)
	tf_activities_counts$var <- var_df
	
	sample_num <- nom %>% dim() %>% .[2]
	tf_activities_counts_filter <- tf_activities_counts %>%
                top_n(sample_num*2, var) %>%
                dplyr::select(-var) %>%
                as.matrix()

	tf_activities_vector <- as.vector(tf_activities_counts_filter)


	paletteLength <- 100
	myColor <- colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

	dorotheaBreaks <- c(seq(min(tf_activities_vector), 0,
	    length.out=ceiling(paletteLength/2) + 1),
	    seq(max(tf_activities_vector)/paletteLength,
	    max(tf_activities_vector),
	    length.out=floor(paletteLength/2)))
	
	pdf(paste0(outdir, "DOROTHEA.", nom_name, ".pdf"), width = 6, height = 6)
	pheatmap(tf_activities_counts_filter,
	    cluster_rows=FALSE, cluster_cols=FALSE,
	    fontsize=14, fontsize_row = 8, fontsize_col = 8,
	    color=myColor, breaks = dorotheaBreaks,
	    main = "Dorothea ABC", angle_col = 90,
	    treeheight_col = 0,  border_color = NA)
	dev.off()

	}


run_dorothea(nom_list[[1]], "all")

