# # library -----------------------------------------------------------------
# library(tidyverse)
# library(fgsea)
# library(msigdbr)
# library(GSEABase)
# library(patchwork)
# library(ComplexHeatmap)
# library(circlize)
# 
# # -------------------------------------------------------------------------
# ddsHTSeq_filter <- readRDS("../../out/object/dds_filter.rds")
# 
# LUT_sample <- colData(ddsHTSeq_filter) %>% 
#   data.frame()
# 
# df_tables_GSEA_all_non_redundant <- read_tsv(file = "../../out/table/df_table_GSEA_res_DIAvsCTRL_shr_nonredundant_KEGG.tsv")
# 
# vds_filter <- readRDS(file = "../../out/object/vds_filter.rds")
# 
# # -------------------------------------------------------------------------
# # plot the heatmap of the leading edges for endo cells
# df_test <- df_tables_GSEA_all_non_redundant %>%
#   filter(dataset == "res_DIAvsCTRL_shr")
# 
# df_leading_edges <- df_test %>% 
#   pull(leadingEdge) %>%
#   str_split(pattern = "\\|") %>%
#   setNames(df_test$pathway)
# 
# # pull the fold change data, I could pull the level of expression for the pseudobulk but I would have just one sample, therefore the row normalization would just make the plot useless
# # filter the one in the leading edges that are significant
# GOI <- results$res_DIAvsCTRL_shr %>%
#   dplyr::filter(symbol %in% df_leading_edges$KEGG_STARCH_AND_SUCROSE_METABOLISM) %>%
#   pull(symbol)
# subset_genes <- data.frame(symbol = pathways$KEGG_STARCH_AND_SUCROSE_METABOLISM) %>% 
#   mutate(col = case_when(symbol %in%  GOI ~ 1,
#                    T~0))
# 
# # gene_id <- rownames(assay(vds_filter)) %in% df_leading_edges$signature_tip_goveia_2020
# gene_id <- subset_genes %>% filter(col==1) %>% pull(symbol)
# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()
# 
# mat_shr <- mat_filter[rownames(vds_filter) %in% gene_id, ]
# mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr)
# #
# meta_sample <- data.frame(colname = colnames(mat2_shr)) %>% 
#   left_join(LUT_sample,by=c("colname"="sample"))
# 
# # make the column of the matrix more readable
# colnames(mat2_shr) <- meta_sample$Description
# 
# sample_ordered_shr <- meta_sample$condition
# column_ha_shr <- HeatmapAnnotation(treat = sample_ordered_shr,  
#                                    col = list(treat = c("ctrl" = "green", "dia" = "black"))) 
# 
# ht2_shr <- Heatmap(mat2_shr, show_column_names = T,
#                    name = "exp", 
#                    column_title = "DIAshr",
#                    # row_names_gp = gpar(fontsize = 3),
#                    top_annotation = column_ha_shr
#                    # cluster_rows = F, 
#                    # right_annotation = row_ha, 
#                    # row_split = rep(c(1,2,3,4),c(2,3,4,7))
#                    
# ) 
# pdf("../../out/image/heatmap_DEG_KEGG_STARCH_AND_SUCROSE_METABOLISM.pdf",width = 5,height = 3) 
# draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()
# 
# 
# # plot heatmap with rank of genes -----------------------------------------
# # get the top and the bottom genes of the signatures in the rank
# df_rank <- results$res_DIAvsCTRL_shr %>%
#   arrange(desc(log2FoldChange)) %>%
#   mutate(rank = nrow(.):1) %>%
#   mutate(rank2 = 1:nrow(.))
# 
# # filter the ranks of the gene in the signature
# id_UP <- df_rank %>%
#   filter(symbol %in% pathways$KEGG_STARCH_AND_SUCROSE_METABOLISM) %>% 
#   mutate(color = case_when(symbol %in% df_leading_edges$KEGG_STARCH_AND_SUCROSE_METABOLISM~"red",T~"black"))
# 
# head(id_UP)
# 
# # library(ComplexHeatmap)
# m = matrix(df_rank$rank,ncol = 1)
# ha_up = rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$symbol,
#                                       labels_gp = gpar(col = id_UP$color,
#                                                        fontsize = 10)))
# 
# hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up)
# 
# pdf("../../out/image/heatmap_id_DEG_KEGG_STARCH_AND_SUCROSE_METABOLISM.pdf",width = 3,height = 10) 
# draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()
# 
# # -------------------------------------------------------------------------
# # pull the fold change data, I could pull the level of expression for the pseudobulk but I would have just one sample, therefore the row normalization would just make the plot useless
# # filter the one in the leading edges that are significant
# GOI <- results$res_DIAvsCTRL_shr %>%
#   dplyr::filter(symbol %in% df_leading_edges$KEGG_TYPE_I_DIABETES_MELLITUS) %>%
#   pull(symbol)
# subset_genes <- data.frame(symbol = pathways$KEGG_TYPE_I_DIABETES_MELLITUS) %>% 
#   mutate(col = case_when(symbol %in%  GOI ~ 1,
#                          T~0))
# 
# # gene_id <- rownames(assay(vds_filter)) %in% df_leading_edges$signature_tip_goveia_2020
# gene_id <- subset_genes %>% filter(col==1) %>% pull(symbol)
# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()
# 
# mat_shr <- mat_filter[rownames(vds_filter) %in% gene_id, ]
# mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr)
# #
# meta_sample <- data.frame(colname = colnames(mat2_shr)) %>% 
#   left_join(LUT_sample,by=c("colname"="sample"))
# 
# # make the column of the matrix more readable
# colnames(mat2_shr) <- meta_sample$Description
# 
# sample_ordered_shr <- meta_sample$condition
# column_ha_shr <- HeatmapAnnotation(treat = sample_ordered_shr,  
#                                    col = list(treat = c("ctrl" = "green", "dia" = "black"))) 
# 
# ht2_shr <- Heatmap(mat2_shr, show_column_names = T,
#                    name = "exp", 
#                    column_title = "DIAshr",
#                    # row_names_gp = gpar(fontsize = 3),
#                    top_annotation = column_ha_shr
#                    # cluster_rows = F, 
#                    # right_annotation = row_ha, 
#                    # row_split = rep(c(1,2,3,4),c(2,3,4,7))
#                    
# ) 
# pdf("../../out/image/heatmap_DEG_KEGG_TYPE_I_DIABETES_MELLITUS.pdf",width = 5,height = 3) 
# draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()
# 
# 
# # plot heatmap with rank of genes -----------------------------------------
# # get the top and the bottom genes of the signatures in the rank
# df_rank <- results$res_DIAvsCTRL_shr %>%
#   arrange(desc(log2FoldChange)) %>%
#   mutate(rank = nrow(.):1) %>%
#   mutate(rank2 = 1:nrow(.))
# 
# # filter the ranks of the gene in the signature
# id_UP <- df_rank %>%
#   filter(symbol %in% pathways$KEGG_TYPE_I_DIABETES_MELLITUS) %>% 
#   mutate(color = case_when(symbol %in% df_leading_edges$KEGG_TYPE_I_DIABETES_MELLITUS~"red",T~"black"))
# 
# head(id_UP)
# 
# # library(ComplexHeatmap)
# m = matrix(df_rank$rank,ncol = 1)
# ha_up = rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$symbol,
#                                       labels_gp = gpar(col = id_UP$color,
#                                                        fontsize = 10)))
# 
# hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up)
# 
# pdf("../../out/image/heatmap_id_DEG_KEGG_TYPE_I_DIABETES_MELLITUS.pdf",width = 3,height = 10) 
# draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()
