# libraries ---------------------------------------------------------------
library(tidyverse)
library(GSVA)
library(limma)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(msigdbr)
library(ggrepel)

# read in the data --------------------------------------------------------
# dds_filter <- readRDS("../../out/object/ddsHTSeq_filter_DESeq_GMPvsResearch.rds")
# 
# # fix the column name to make it more readable
# # pull the order form the count table
# id_column <- counts(dds_filter,normalized = T) %>%
#   data.frame() |>
#   colnames()
# # pull the new name from the ordered meta
# new_column_name <- colData(dds_filter)%>%
#   data.frame() %>%
#   arrange(factor(sample,levels=id_column)) %>%
#   pull("Sample.name")
# 
# test_01 <- counts(dds_filter,normalized = T) %>%
#   data.frame()
# colnames(test_01) <- new_column_name
# 
# test_01 %>%
#   rownames_to_column(var = "symbol") %>%
#   saveRDS("../../out/object/ddsHTSeq_filter_DESeq_GMPvsResearch_counts_norm.rds")
# 
# test_02 <- counts(dds_filter,normalized = F) %>%
#   data.frame()
# colnames(test_02) <- new_column_name
# 
# test_02 %>%
#   rownames_to_column(var = "symbol") %>%
#   saveRDS("../../out/object/ddsHTSeq_filter_DESeq_GMPvsResearch_counts_raw.rds")

# read in the normalzied expression table geneXsample
exp <- readRDS("../../out/object/ddsHTSeq_filter_DESeq_GMPvsResearch_counts_norm.rds") %>%
  # dplyr::select("symbol") %>%
  column_to_rownames("symbol") %>%
  as.matrix()

# colnames(exp) <- str_sub(colnames(exp),start = 18,end = -1)

# read in the signatures
# pathways_all <- readRDS("data/list_costume_signatures_pathways.rds")
# merge all the pathways in a single list
# pathways <- pathways_all
reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
  # only use a subset of the dataset
  # filter(gs_name %in% c("KEGG_PRIMARY_IMMUNODEFICIENCY",
  #   "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  #   "KEGG_VEGF_SIGNALING_PATHWAY",
  #   "KEGG_CELL_CYCLE",
  #   "KEGG_CELL_ADHESION_MOLECULES_CAMS",
  #   "KEGG_ECM_RECEPTOR_INTERACTION"))
head(reactome_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)

# -------------------------------------------------------------------------
# perform the GSEA based on the normalized counts
es <- gsva(exp,
           pathways,
           min.sz=5,
           max.sz=500,
           kcdf="Poisson",
           mx.diff=TRUE,
           verbose=FALSE,
           parallel.sz=1)

es_log <- gsva(log(exp+1),
               pathways,
               min.sz=5,
               max.sz=500,
               kcdf="Poisson",
               mx.diff=TRUE,
               verbose=FALSE,
               parallel.sz=1)

# show correlation between the two estimates ------------------------------
left_join(
  es %>%
    data.frame() %>%
    rownames_to_column("pathway") %>%
    pivot_longer(names_to = "sample",values_to = "estimate",-pathway),
  es_log %>%
    data.frame() %>%
    rownames_to_column("pathway") %>%
    pivot_longer(names_to = "sample",values_to = "estimate_log",-pathway),by = c("pathway","sample")) %>%
  ggplot(aes(x=estimate_log,y=estimate))+geom_point()

# STATISTICAL TESTING -----------------------------------------------------
# use limma for testing significance
# library(limma)
# adjPvalueCutoff <- 0.001
# logFCcutoff <- log2(2)
# logFCcutoff

# define the factor for the treatmentnt basesd on the colnames
lut <- data.frame(colname = colnames(es)) %>%
  mutate(treat = unlist(str_extract_all(colname,pattern = "GMP|Research"))) %>%
  mutate(treat = factor(treat,levels = c("Research","GMP")))

design <- model.matrix(~ lut$treat)
colnames(design) <- c("intercept", "GMPVsResearch")
design

fit <- lmFit(es, design)
fit_log <- lmFit(es_log, design)
fit <- eBayes(fit)
fit_log <- eBayes(fit_log)

allGenesets <- topTable(fit, coef="GMPVsResearch", number=Inf)
allGenesets_log <- topTable(fit_log, coef="GMPVsResearch", number=Inf)

# volcano GSVA ------------------------------------------------------------
df_plot <- allGenesets %>%
  rownames_to_column("pathway") %>%
  # add color to the text
  # mutate(color = case_when(pathway %in% c("GOBP_RECEPTOR_SIGNALING_PATHWAY_VIA_STAT")~"red",
  #                          T~"black")) %>%
  mutate(color = case_when(pathway %in% c("KEGG_PRIMARY_IMMUNODEFICIENCY",
                                          "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
                                          "KEGG_VEGF_SIGNALING_PATHWAY",
                                          "KEGG_CELL_CYCLE",
                                          "KEGG_CELL_ADHESION_MOLECULES_CAMS",
                                          "KEGG_ECM_RECEPTOR_INTERACTION")~"red",
                           T~"black")) %>%
  mutate(color = factor(color))

df_plot2 <- df_plot %>%   
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_|KEGG_|GOBP_") %>%
           str_sub(start = 1,end = 35))
  # mutate(min_log10_padj = -log10(padj)) %>%
df_plot2 %>%
  ggplot(aes(y = -log10(adj.P.Val),x = logFC,label = pathway2,col=color)) +
  geom_point(alpha = 0.2) +
  # geom_point(aes(size = size),alpha = 0.2) +
  # facet_wrap(~dataset) +
  # theme_bw(base_rect_size = 2)+
  theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA),
  #       axis.ticks = element_line(colour = "black"),
  #       #axis.ticks.length = unit(.25, "cm")
  #       legend.position = "none"
  # )+
  theme(strip.background = element_blank(),legend.position = "none")+
  geom_text_repel(data = df_plot2 %>% filter(color == "black"),aes(-log10(adj.P.Val),x = logFC,label = pathway2,col=color),
                  size = 2,
                  box.padding = 0.5,
                  segment.alpha = 0.5,
                  max.overlaps = 10)+
  geom_text_repel(data = df_plot2 %>% filter(color == "red"),aes(-log10(adj.P.Val),x = logFC,label = pathway2,col=color),
                  size = 2,
                  box.padding = 0.5,
                  segment.alpha = 0.5,
                  max.overlaps = 10,min.segment.length = 0,nudge_x = 0.2,nudge_y = 0.2)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray",alpha=0.8)+
  scale_color_manual(values = levels(df_plot$color))
ggsave("../../out/image/volcano_GSVA_KEGG_GMPVsResearch.pdf",width = 12,height = 6)

# plot heatmap ------------------------------------------------------------
# library(ComplexHeatmap)
# define only the significnat terms
DEgeneSets <- allGenesets %>%
  rownames_to_column("pathway") %>%
  filter(adj.P.Val < 0.05) %>%
  pull(pathway) %>%
  # force in the one of interest
  c("KEGG_PRIMARY_IMMUNODEFICIENCY",
      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
      "KEGG_VEGF_SIGNALING_PATHWAY",
      "KEGG_CELL_CYCLE",
      "KEGG_CELL_ADHESION_MOLECULES_CAMS",
      "KEGG_ECM_RECEPTOR_INTERACTION")

mat <- es
mat_norm <- es %>%
  data.frame() %>%
  rownames_to_column() %>%
  # plot only the significant terms
  filter(rowname %in% DEgeneSets) %>%
  # scale the values rowwise
  gather(key = sample,value = exp,-rowname) %>%
  group_by(rowname) %>%
  mutate(norm = (exp - mean(exp))/sd(exp)) %>%
  dplyr::select(-exp) %>%
  spread(key = sample,value = norm) %>%
  column_to_rownames()

sample_ordered <- str_extract(colnames(mat_norm),pattern = "GMP|Research")
sample_ordered

# build the annotation object
column_ha <- HeatmapAnnotation(treat = sample_ordered,  
                               col = list(treat = c("GMP" = "green", "Research"="black"))) 

# highlight some terms
rowname_color <- case_when(rownames(mat_norm)%in%c("KEGG_PRIMARY_IMMUNODEFICIENCY",
                                                   "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
                                                   "KEGG_VEGF_SIGNALING_PATHWAY",
                                                   "KEGG_CELL_CYCLE",
                                                   "KEGG_CELL_ADHESION_MOLECULES_CAMS",
                                                   "KEGG_ECM_RECEPTOR_INTERACTION")~"red",
                           T~"black")

hm <- Heatmap(mat_norm,
              # add annotation for the columns
              # hide columns labels
              # show_column_names = F,
              # fix width of the lables
              top_annotation = column_ha,
              row_names_gp = gpar(col = rowname_color),
              row_names_max_width = max_text_width(
                rownames(mat),
                gp = gpar(fontsize = 12)
              ))

pdf(file = "../../out/image/heatmap_GSVA_KEGG_GMPVsResearch_tailored_es_nonLog_ZScore.pdf", width = 16, height = 5)
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 100), "mm"))
dev.off()