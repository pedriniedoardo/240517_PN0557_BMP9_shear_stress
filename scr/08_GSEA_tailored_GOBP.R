# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)

# pull the expression values ----------------------------------------------
ddsHTSeq_filter <- readRDS("../../out/object/dds_all_filter_GMPvsResearch.rds")

LUT_sample <- colData(ddsHTSeq_filter) %>%
  data.frame()

df_tables_GSEA_all_non_redundant <- read_tsv(file = "../../out/table/df_tables_GSEA_res_GMPvsRESEARCH_shr_all_GOBP.tsv")

vds_filter <- readRDS("../../out/object/vds_all_filter_GMP vs RESEARCH ECFCs.rds")

# pull the DE results -----------------------------------------------------
file <- dir("../../out/table/") %>%
  str_subset(pattern = "res_") %>%
  str_subset(pattern = ".txt") %>%
  str_subset(pattern = "shr",negate = F) %>%
  str_subset(pattern = "GMPvsRESEARCH",negate = F)
file 

# load the results 
results <- lapply(paste0("../../out/table/",file),function(x){
  read_tsv(x) 
}) %>%
  setNames(str_remove_all(file,pattern = ".txt"))

# generate the list of ranks
list_ranks <- lapply(results, function(x){
  
  x <- dplyr::filter(x,!is.na(symbol)) %>%
    # average logFC in case of duplicated genenames
    group_by(symbol) %>%
    summarise(logFC = mean(log2FoldChange))
  
  ranks <- setNames(x$logFC, x$symbol)
  ranks
}) 
glimpse(list_ranks)

# pull the pathways reference ---------------------------------------------
GOBP_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
head(GOBP_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = GOBP_gene_sets$gene_symbol, f = GOBP_gene_sets$gs_name)
# head(pathways)

# wrangling ---------------------------------------------------------------
# plot the heatmap of the leading edges for endo cells
df_test <- df_tables_GSEA_all_non_redundant %>%
  filter(dataset == "res_GMPvsRESEARCH_shr")

"HALLMARK_GLYCOLYSIS"

# identify is the terms of interest are in the dataset
df_test |> 
  filter(str_detect(pathway,pattern = c("GOBP_ENDOTHELIAL_CELL_MATRIX_ADHESION|
  GOBP_POSITIVE_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS|GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS|GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_MIGRATION|GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_CHEMOTAXIS|GOBP_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE|GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_PROLIFERATION|GOBP_REGULATION_OF_ENDOTHELIAL_TUBE_MORPHOGENESIS|GOBP_ENDOTHELIAL_CELL_MORPHOGENESIS|GOBP_VASCULOGENESIS|GOBP_ESTABLISHMENT_OF_ENDOTHELIAL_BARRIER|GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_DIFFERENTIATION")))

# plot the GSEA profile ---------------------------------------------------
# plotEnrichment(pathways$KEGG_DNA_REPAIR, list_ranks$res_HYPvsNORM_shr) + labs(title = "KEGG_DNA_REPAIR") 
# ggsave(filename = "../../out/image/profile_GSEA_KEGG_DNA_REPAIR.pdf",width = 6,height = 4)

# plot volcano GSEA highlighlght ------------------------------------------
df_plot <- df_tables_GSEA_all_non_redundant %>%
  # add color to the text
  # mutate(color = case_when(pathway %in% c("GOBP_RECEPTOR_SIGNALING_PATHWAY_VIA_STAT")~"red",
  #                          T~"black")) %>%
  mutate(color = case_when(pathway %in% c("GOBP_ENDOTHELIAL_CELL_MATRIX_ADHESION",
                                          "GOBP_POSITIVE_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
                                          "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
                                          "GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_MIGRATION",
                                          "GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_CHEMOTAXIS",
                                          "GOBP_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE",
                                          "GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_PROLIFERATION",
                                          "GOBP_REGULATION_OF_ENDOTHELIAL_TUBE_MORPHOGENESIS",
                                          "GOBP_ENDOTHELIAL_CELL_MORPHOGENESIS",
                                          "GOBP_VASCULOGENESIS",
                                          "GOBP_ESTABLISHMENT_OF_ENDOTHELIAL_BARRIER",
                                          "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_DIFFERENTIATION")~"red",
                           T~"black")) %>%
  mutate(color = factor(color))

df_plot2 <- df_plot %>%   
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_|KEGG_|GOBP_") %>%
           str_sub(start = 1,end = 35))

df_plot2 %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES)) + 
  geom_point(data = df_plot2 %>% filter(color == "black"),aes(size = size),alpha = 0.2,color="black") +
  geom_point(data = df_plot2 %>% filter(color == "red"),aes(size = size),alpha = 0.2,color="red") +
  facet_wrap(~dataset) +
  # theme_bw(base_rect_size = 2)+
  theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA),
  #       axis.ticks = element_line(colour = "black"),
  #       #axis.ticks.length = unit(.25, "cm")
  #       legend.position = "none"
  # )+
  theme(strip.background = element_blank(),legend.position = "none")+
  geom_text_repel(data = df_plot2 %>% filter(color == "red"),aes(y = -log10(padj),x = NES,label = pathway2,col=color),
                  size = 2,
                  box.padding = 0.5,
                  segment.alpha = 0.5,
                  max.overlaps = 10,min.segment.length = 0,nudge_x = 0.5,nudge_y = 0.5)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray",alpha=0.8)+
  # scale_color_manual(values = levels(df_plot$color))
  scale_color_manual(values = "red")
ggsave("../../out/image/volcano_GSEA_GOBP_GMPvsRESEARCH_tailored.pdf",width = 8,height = 5)

# plot the leading edges --------------------------------------------------
df_leading_edges <- df_test %>%
  pull(leadingEdge) %>%
  str_split(pattern = "\\|") %>%
  setNames(df_test$pathway)

# pp <- "KEGG_PRIMARY_IMMUNODEFICIENCY"
list_hm <- lapply(c("GOBP_ENDOTHELIAL_CELL_MATRIX_ADHESION",
                    "GOBP_POSITIVE_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
                    "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
                    "GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_MIGRATION",
                    "GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_CHEMOTAXIS",
                    "GOBP_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE",
                    "GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_PROLIFERATION",
                    "GOBP_REGULATION_OF_ENDOTHELIAL_TUBE_MORPHOGENESIS",
                    "GOBP_ENDOTHELIAL_CELL_MORPHOGENESIS",
                    "GOBP_VASCULOGENESIS",
                    "GOBP_ESTABLISHMENT_OF_ENDOTHELIAL_BARRIER",
                    "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_DIFFERENTIATION"),
  function(pp){
    
    # pull the fold change data, I could pull the level of expression for the pseudobulk but I would have just one sample, therefore the row normalization would just make the plot useless
    # filter the one in the leading edges that are significant
    GOI <- results$res_GMPvsRESEARCH_shr %>%
      dplyr::filter(symbol %in% df_leading_edges[[pp]]) %>%
      pull(symbol)
    subset_genes <- data.frame(symbol = pathways[[pp]]) %>%
      mutate(col = case_when(symbol %in%  GOI ~ 1,
                             T~0))
    
    # gene_id <- rownames(assay(vds_filter)) %in% df_leading_edges$signature_tip_goveia_2020
    gene_id <- subset_genes %>% filter(col==1) %>% pull(symbol)
    mat_filter <- assay(vds_filter) %>%
      data.frame() %>%
      # dplyr::select(contains(c("_0_","_6_"))) %>%
      as.matrix()
    
    mat_shr <- mat_filter[rownames(vds_filter) %in% gene_id, ]
    mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr,useNames = T)
    #
    meta_sample <- data.frame(colname = colnames(mat2_shr)) %>%
      left_join(LUT_sample,by=c("colname"="sample"))
    
    # make the column of the matrix more readable
    colnames(mat2_shr) <- meta_sample$Sample.name
    
    sample_ordered_shr <- meta_sample$Condition
    column_ha_shr <- HeatmapAnnotation(treat = sample_ordered_shr,
                                       col = list(treat = c("GMP" = "green",
                                                            "Research" = "black")))
    
    ht2_shr <- Heatmap(mat2_shr, show_column_names = T,
                       name = "exp",
                       column_title = pp,
                       # row_names_gp = gpar(fontsize = 3),
                       top_annotation = column_ha_shr
                       # cluster_rows = F,
                       # right_annotation = row_ha,
                       # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                       
    )
    return(draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left"))
    
  })

pdf("../../out/image/heatmap_DEG_GOBP_GMPvsRESEARCH_leadingEdges.pdf",width = 7,height = 8)
lapply(list_hm,function(x){
  x
})
dev.off()

# # plot heatmap with rank of genes -----------------------------------------
# # get the top and the bottom genes of the signatures in the rank
# df_rank <- results$res_HYPvsNORM_shr %>%
#   arrange(desc(log2FoldChange)) %>%
#   mutate(rank = nrow(.):1) %>%
#   mutate(rank2 = 1:nrow(.))
# 
# # filter the ranks of the gene in the signature
# id_UP <- df_rank %>%
#   filter(symbol %in% pathways$KEGG_DNA_REPAIR) %>%
#   mutate(color = case_when(symbol %in% df_leading_edges$KEGG_DNA_REPAIR~"red",T~"black")) |> 
#   filter(color != "black")
# 
# head(id_UP)
# 
# # library(ComplexHeatmap)
# m <- matrix(df_rank$rank,ncol = 1)
# ha_up <- rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$symbol,
#                                        labels_gp = gpar(col = id_UP$color,
#                                                         fontsize = 10),
#                                        link_width = unit(20, "mm"),
#                                        extend = unit(100, "mm")))
# 
# hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up)
# 
# pdf("../../out/image/heatmap_id_DEG_KEGG_DNA_REPAIR.pdf",width = 3,height = 30)
# draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(100, 2, 2, 2), "mm"))
# dev.off()
# 
# 
# ha_up <- rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$symbol,
#                                        labels_gp = gpar(col = id_UP$color,
#                                                         fontsize = 10)))
# 
# hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up)
# 
# pdf("../../out/image/heatmap_id_DEG_KEGG_DNA_REPAIR.pdf",width = 3,height = 30)
# draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left")
# dev.off()