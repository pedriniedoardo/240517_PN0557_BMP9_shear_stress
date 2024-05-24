# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(DESeq2)
# library(ComplexHeatmap)
# library(UpSetR)
# library(gplots)
# # library(msigdbr)
# 
# # plot genes with heatmap -------------------------------------------------
# vds_filter <- readRDS(file = "../../out/object/vds_filter.rds")
# 
# # load the table of degs
# df_deg <- read_tsv("../../out/table/res_PROCRvsSCR_shr.txt")
# 
# # read the whole signatures
# # kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
# 
# # load the genes of interest
# df_signatures <- read_tsv("../../out/table/df_tables_GSEA_res_PROCRvsSCR_shr_all_KEGG.tsv")
# 
# id_signature <- c("KEGG_FOCAL_ADHESION",
#                   "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON",
#                   "KEGG_MAPK_SIGNALING_PATHWAY",
#                   "KEGG_ADHERENS_JUNCTION")
# 
# list_sigantures <- lapply(id_signature,function(x){
#            df_signatures %>% 
#              filter(pathway %in% x) %>% 
#              pull(leadingEdge) %>%
#              str_split(pattern = "\\|") %>% 
#              unlist() %>% 
#              data.frame(gene = .) %>% 
#              mutate(siganture = x) %>% 
#              left_join(df_deg,by = c("gene"="symbol"))
#          }) %>% 
#   setNames(id_signature)
# 
# # regular heatmap ---------------------------------------------------------
# # x<-list_sigantures$KEGG_FOCAL_ADHESION
# # name_sig <- "KEGG_FOCAL_ADHESION"
# 
# pmap(list(list_sigantures,names(list_sigantures)),function(x,name_sig){
#   gene_id <- rownames(assay(vds_filter)) %in% x$gene
#   gene_id
#   
#   # define the matrices
#   mat <- assay(vds_filter)[gene_id, ]
#   mat2 <- (mat - rowMeans(mat))/rowSds(mat)
#   
#   mat2
#   
#   # build the annotation object  
#   meta <- colData(vds_filter) %>% 
#     data.frame()
#   
#   LUT_sample_matrix <- data.frame(sample_matrix = colnames(mat2)) %>%
#     left_join(meta,by = c("sample_matrix"="sample"))
#   
#   sample_ordered <- LUT_sample_matrix$treat
#   
#   column_ha <- HeatmapAnnotation(treat = sample_ordered
#                                  # col = list(treat = c("A" = "green", "B" = "gray","D"="black"))
#   ) 
#   
#   # change the name of the column in the matrix
#   colnames(mat2) <- LUT_sample_matrix$Description
#   
#   # label the significant genes
#   gene_significant <- df_deg %>% 
#     filter(abs(log2FoldChange)>1&padj<0.05) %>% 
#     pull(symbol)
#   
#   rowname_color <- case_when(rownames(mat2)%in%gene_significant~"red",
#                              T~"black")
#   # rowname_color
#   
#   ht2 <- Heatmap(mat2, 
#                  name = "NES",
#                  top_annotation = column_ha, 
#                  row_names_gp = gpar(col = rowname_color),
#                  # cluster_rows = F, 
#                  # col = colorRamp2(c(-2, 0, 1), c("green", "white", "red")),
#                  # right_annotation = row_ha, 
#                  # row_split = rep(c(1,2,3,4),c(2,3,4,7)),
#                  column_title = name_sig) 
#   
#   pdf(paste0("../../out/image/heatmap_GOI_LeadingEdges_",name_sig,".pdf"),width = 10,height = 15) 
#   draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
#   dev.off()
# })
# 
# # # rank heatmap ------------------------------------------------------------
# # 
# # df_rank <- df_deg %>% 
# #   arrange(desc(log2FoldChange)) %>% 
# #   mutate(rank = nrow(.):1) %>% 
# #   mutate(rank2 = 1:nrow(.))
# # 
# # df_rank 
# # 
# # id_UP <- df_rank %>% 
# #   filter(symbol %in% df_signatures)
# # 
# # id_UP 
# # # A tibble: 31 x 10 
# # symbol  baseMean log2FoldChange  lfcSE   pvalue     padj ensembl         entrez  rank rank2 
# # <chr>      <dbl>          <dbl>  <dbl>    <dbl>    <dbl> <chr>            <dbl> <int> <int> 
# #   1 LPL         364.         1.21   0.226  2.26e- 9 1.03e- 7 ENSG00000175445   4023 22230    87 
# # 2 VCAN       1033.         1.18   0.525  7.53e- 5 1.31e- 3 ENSG00000038427   1462 22225    92 
# # 3 VEGFA      3153.         1.16   0.179  2.34e-12 1.56e-10 ENSG00000112715   7422 22220    97 
# # 4 JAG2       2246.         0.737  0.104  2.12e-14 1.85e-12 ENSG00000184916   3714 22122   195 
# # 5 JAG1       6684.         0.671  0.223  3.06e- 5 5.92e- 4 ENSG00000101384    182 22097   220 
# # 6 SLCO2A1   10520.         0.536  0.524  2.20e- 3 2.44e- 2 ENSG00000174640   6578 22032   285 
# # 7 TIMP1      3183.         0.395  0.0718 4.81e- 9 2.05e- 7 ENSG00000102265   7076 21894   423 
# # 8 POSTN       987.         0.365  0.527  6.30e- 3 5.82e- 2 ENSG00000133110  10631 21850   467 
# # 9 APP       89813.         0.346  0.0504 8.07e-13 5.68e-11 ENSG00000142192    351 21824   493 
# # 10 VAV2       1504.         0.0854 0.0761 2.39e- 2 1.65e- 1 ENSG00000160293   7410 20873  1444 
# # # ... with 21 more rows
# # ```
# # produce the heatmap
# # ```
# # m = matrix(df_rank$rank,ncol = 1) 
# # ha_up = rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$symbol)) 
# # hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up)
# # 
# # draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left")
# # ```
# # 
# # 
# # 
# # 
