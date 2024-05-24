# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(limma)
# library(AnnotationDbi) 
# library(AnnotationHub)
library(ComplexHeatmap)
library(ashr)

# build the object for the annotation  ------------------------------------
# ah <- AnnotationHub()
# ah
# 
# #query the database to identify the version of interest
# query(ah, pattern = c("Homo Sapiens", "EnsDb"))
# 
# # copy the specific code ofr the database of interest
# edb <- ah[["AH109606"]]
# columns(edb)

# read in the biotype annotation
LUT_gene <- read_tsv("../../out/table/df_LUT_gtf_genes.tsv") %>% 
  mutate(transcript_id=gene_ids) %>% 
  separate(transcript_id,into = c("gene_id","version"),remove = F,sep = "\\.") %>% 
  dplyr::select(-c(version,gene_ids))

# read in the data --------------------------------------------------------
ddsHTSeq_filter <- readRDS("../../out/object/dds_all_filter_GMPvsResearch.rds")
design <- readRDS("../../out/object/design_all_dds_GMPvsResearch.rds")
vds_filter <- readRDS("../../out/object/vds_all_filter_GMP vs RESEARCH ECFCs.rds")

LUT_sample <- colData(ddsHTSeq_filter) %>% 
  data.frame()

# differential expression analyisis ---------------------------------------
ddsHTSeq_filter <- DESeq(ddsHTSeq_filter)
# if needed is possible to check the distributions of the counts before and after the normalizatoin
# boxplot(log(counts(ddsHTSeq_structure,normalized = T)))
# boxplot(log(counts(ddsHTSeq_structure,normalized = F)))

# save the filtered object
saveRDS(ddsHTSeq_filter,"../../out/object/ddsHTSeq_filter_DESeq_GMPvsResearch.rds")

# print the contrast
resultsNames(ddsHTSeq_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(GMPvsRESEARCH = conditionGMP,
                          levels = design)

res <- results(ddsHTSeq_filter, contrast=contrast[,"GMPvsRESEARCH"],alpha = 0.05)

summary(res)

# add the gene symbols
list_df <- 
  list("GMPvsRESEARCH" = res) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue) %>%
      left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
    })

# check the expression of Procr
list_df$GMPvsRESEARCH %>% 
  filter(str_detect(symbol,pattern = "PTX3"))

# save the tables
pmap(list(list_df,names(list_df)),function(x,y){
  name <- paste0("res_",y,".txt")
  # name
  write_tsv(x,file = paste0("../../out/table/",name))
})

# shrink ------------------------------------------------------------------  
res_shr <- lfcShrink(ddsHTSeq_filter, res = res, type = "ashr")
summary(res_shr)

list_df_shr <- 
  list("GMPvsRESEARCH_shr" = res_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue) %>%
      left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  })

# save the tables
pmap(list(list_df_shr,names(list_df_shr)),function(x,y){
  name <- paste0("res_",y,".txt")
  # name
  write_tsv(x,file = paste0("../../out/table/",name))
})

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
list_df$GMPvsRESEARCH %>%
  data.frame()%>%
  dplyr::filter(baseMean>1)%>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  theme_bw()
ggsave("../../out/image/histogram_pvalue_GMPvsRESEARCH.pdf",width = 4,height = 3)

# PLOTTING RESULTS --------------------------------------------------------
# add the info of the genename
test_plot <- list_df$GMPvsRESEARCH %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0))

test_plot %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = test_plot[test_plot$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = test_plot[test_plot$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = test_plot[test_plot$col==1,][1:1000,],
    aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_text_GMPvsRESEARCH.pdf",width = 12,height = 12)
#
# test_plot_fake <- list_df$HYPvsNORM %>%
#   data.frame()%>%
#   # add a clor variable in case significant
#   mutate(col=ifelse(((pvalue<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0))
# 
# test_plot_fake %>%
#   ggplot(aes(x=log2FoldChange,y=-log(pvalue)))+
#   # geom_point()
#   geom_point(data = test_plot_fake[test_plot_fake$col==0,],aes(x=log2FoldChange,y=-log(pvalue),col=factor(col)),alpha=0.05)+
#   geom_point(data = test_plot_fake[test_plot_fake$col==1,],aes(x=log2FoldChange,y=-log(pvalue),col=factor(col)),alpha=0.5)+
#   geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
#   geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
#   scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
#   ggrepel::geom_text_repel(
#     data = test_plot_fake[test_plot_fake$col==1,],
#     aes(label = symbol),
#     size = 5,
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines")) +
#   theme_bw() +
#   theme(legend.position = "none")
# ggsave("../../out/image/vulcano_plot_text_DIAvsCTRL_fake.pdf",width = 12,height = 12)

#
test_plot_shr <- list_df_shr$GMPvsRESEARCH_shr%>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0))

test_plot_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = test_plot_shr[test_plot_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = test_plot_shr[test_plot_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = test_plot_shr[test_plot_shr$col==1,],
    aes(label = symbol),max.overlaps = 3,segment.alpha=0.4,
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_text_GMPvsRESEARCH_shr.pdf",width = 12,height = 12)

# plotMA(res, ylim = c(-5, 5))
list_df$GMPvsRESEARCH %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>% 
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) + 
  scale_x_log10() + scale_color_manual(values = c("gray","red")) + theme_bw() + 
  theme(legend.position = "none")
ggsave("../../out/image/MA_plot_GMPvsRESEARCH.pdf",width = 4,height = 3)

# using the shrinked values
list_df_shr$GMPvsRESEARCH_shr %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>% 
  ggplot(aes(baseMean,y = log2FoldChange)) + 
  geom_point(aes(col=color),alpha=0.2) + 
  ggrepel::geom_text_repel(
    data = test_plot_shr[test_plot_shr$col==1,],
    aes(label = symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  scale_x_log10() + 
  scale_color_manual(values = c("gray","red")) + 
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/image/MA_plot_shr_GMPvsRESEARCH.pdf",width = 12,height = 12)

# # the DEGs plot stringent
# DEG_1 <- list_df$ %>%
#   data.frame()%>%
#   # add a clor variable in case significant
#   mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
#   filter(col==1) %>%
#   pull(symbol)

mat_filter <- assay(vds_filter) %>%
  data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# mat <- mat_filter[rownames(vds_filter) %in% DEG_1, ]
# mat2 <- (mat - rowMeans(mat))/rowSds(mat)
# #
# 
# sample_ordered <- str_extract(colnames(mat2),pattern = "BMP9|mock")
# column_ha <- HeatmapAnnotation(treat = sample_ordered,  
#                                col = list(treat = c("mock" = "green", "BMP9" = "gray"))) 
# 
# ht2 <- Heatmap(mat2, 
#                name = "exp", 
#                column_title = "BMP9",
#                row_names_gp = gpar(fontsize = 3),
#                top_annotation = column_ha, 
#                # cluster_rows = F, 
#                # right_annotation = row_ha, 
#                # row_split = rep(c(1,2,3,4),c(2,3,4,7))
# ) 
# pdf("out/image/heatmap_BMP9_DEG.pdf",width = 4,height = 15) 
# draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()

# the DEGs plot
DEG_2 <- list_df_shr$GMPvsRESEARCH_shr %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(symbol)

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

mat_shr <- mat_filter[rownames(vds_filter) %in% DEG_2, ]
mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr,useNames = TRUE)
#
meta_sample <- data.frame(colname = colnames(mat2_shr)) %>% 
  left_join(LUT_sample,by=c("colname"="sample"))

# make the column of the matrix more readable
colnames(mat2_shr) <- meta_sample$Sample.name

column_ha_shr <- HeatmapAnnotation(treat = meta_sample$Condition,  
                                   col = list(treat = c("GMP" = "green", "Research" = "black"))) 

ht2_shr <- Heatmap(mat2_shr, show_column_names = T,
                   name = "exp", 
                   column_title = "GMPvsRESEARCHshr",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr,show_row_names = F,
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                  
) 
pdf("../../out/image/heatmap_DEG_shr_GMPvsRESEARCH.pdf",width = 7,height = 8) 
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# PLOT DISPERSION ---------------------------------------------------------
pdf("../../out/image/ddsHTSeq_filter_dispersion_GMPvsRESEARCH.pdf",width = 5,height = 5) 
plotDispEsts(ddsHTSeq_filter)
dev.off()