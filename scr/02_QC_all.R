# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(vsn)
library(hexbin)
library(viridis)
library(pheatmap)
library(PoiClaClu)
library(AnnotationDbi) 
library(AnnotationHub)
library(GGally)
library(RNAseqQC)

# read in the data --------------------------------------------------------
ddsHTSeq <- readRDS("../../out/object/dds_all_raw.rds")

# remove low epressed genes -----------------------------------------------
colSums(counts(ddsHTSeq)) %>%
  data.frame(tot_counts=.) %>%
  rownames_to_column("sample") %>% 
  # add the metadata
  left_join(colData(ddsHTSeq) %>% 
              data.frame(),by = "sample") %>% 
  ggplot(aes(x=Sample.name,y = tot_counts))+geom_col()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = margin(0.5, 0.5, 2, 2, "cm"))
ggsave("../../out/image/barplot_tot_count_all.pdf",width = 7,height = 5)

nrow(ddsHTSeq)
# remove potential non infirmative genes
ddsHTSeq_filter <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 100, ]

nrow(ddsHTSeq_filter)

# save the filtered object
saveRDS(ddsHTSeq_filter,file = "../../out/object/dds_all_filter.rds")

# scaling transformation of the data --------------------------------------
# vsd
vds_filter <- vst(ddsHTSeq_filter, blind = F)
vsd_blind <- vst(ddsHTSeq_filter, blind = T)
head(assay(vds_filter), 3)
head(assay(vsd_blind), 3)

# rlog
rld_filter <- rlog(ddsHTSeq_filter, blind = F)
rld_blind <- rlog(ddsHTSeq_filter, blind = T)
head(assay(rld_filter), 3)
head(assay(rld_blind), 3)

meanSdPlot_vsd <- meanSdPlot(assay(vds_filter))
ggsave(plot = meanSdPlot_vsd$gg+theme_bw(),filename = "../../out/image/meanSdPlot_vsd.pdf",width = 4,height = 4)

meanSdPlot_rlog <- meanSdPlot(assay(rld_filter))
ggsave(plot = meanSdPlot_rlog$gg+theme_bw(),filename = "../../out/image/meanSdPlot_rlog.pdf",width = 4,height = 4)

ddsHTSeq_filter <- estimateSizeFactors(ddsHTSeq_filter)

# sample distance ---------------------------------------------------------
sampleDists_vsd <- dist(t(assay(vds_filter)))
sampleDists_vsd

head(assay(vds_filter))

sampleDistMatrix_vsd <- as.matrix(sampleDists_vsd)

rownames(sampleDistMatrix_vsd) <- paste(vds_filter$`Sample name`)
colnames(sampleDistMatrix_vsd) <- NULL

map_colors<-colorRampPalette(viridis(12))(255)

hm_1 <- pheatmap(sampleDistMatrix_vsd,
                 clustering_distance_rows = sampleDists_vsd,
                 clustering_distance_cols = sampleDists_vsd,
                 col = map_colors)

pdf("../../out/image/heatmap_vsd.pdf",width = 7,height = 5)
hm_1
dev.off()

sampleDists_rld <- dist(t(assay(rld_filter)))
sampleDists_rld

head(assay(rld_filter))
head(assay(rld_blind))

sampleDistMatrix_rld <- as.matrix(sampleDists_rld)

rownames(sampleDistMatrix_rld) <- paste(rld_filter$`Sample name`)
colnames(sampleDistMatrix_rld) <- NULL

hm_1_2 <- pheatmap(sampleDistMatrix_rld,
                   clustering_distance_rows = sampleDists_rld,
                   clustering_distance_cols = sampleDists_rld,
                   col = map_colors)

pdf("../../out/image/heatmap_rld.pdf",width = 7,height = 5)
hm_1_2
dev.off()

poisd <- PoissonDistance(t(counts(ddsHTSeq_filter,normalized = F)))

samplePoisDistMatrix <- as.matrix(poisd$dd)

rownames(samplePoisDistMatrix) <- paste(ddsHTSeq_filter$`Sample name`)
colnames(samplePoisDistMatrix) <- NULL

hm_p <- pheatmap(samplePoisDistMatrix,
                 clustering_distance_rows = poisd$dd,
                 clustering_distance_cols = poisd$dd,
                 col = map_colors)

pdf("../../out/image/heatmap_poisd.pdf",width = 7,height = 5)
hm_p
dev.off()

# plot cluster alternative ------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/image/heatmap_cluster_all.pdf",width = 7,height = 5)
set.seed(1)
plot_sample_clustering(vds_filter,
                       anno_vars = c("Condition","Gender","Passage"),
                       distance = "euclidean")
dev.off()

# PCA plot ----------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter,
                    intgroup = c("sample","Description","Sample name","Donor","Cell Type","Condition","Gender","Passage","Maternal Ethnicity","Paternal Ethnicity","Donor Ethnicity","Source","Approx Time between collection and processing","sizeFactor")) +
  theme_bw()

plot_vsd$data %>%
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=Condition,shape=Gender)) +
  geom_point(size =3) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/image/PCA_vsd.pdf",width = 5,height = 4)

# pull more PC
rv <- rowVars(assay(vds_filter),useNames = TRUE)
# select the ntop genes by variance
select_var <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
test <- prcomp(t(assay(vds_filter)[select_var,]))$x %>% 
  data.frame() %>% 
  rownames_to_column("sample")

left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample")) %>% 
  # ggpairs(columns = 5:14,ggplot2::aes(colour=condition),upper = "blank")+
  ggpairs(columns = 5:12,ggplot2::aes(colour=Condition),upper = "blank") +
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/image/PCA_panel_vsd.pdf",width = 20,height = 20)

test <- left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample"))
test_df1 <- test |> 
  dplyr::select(sample,Description,Sample.name,Cell.Type:Passage,Approx.Time.between.collection.and.processing) |> 
  dplyr::rename(approx.time = Approx.Time.between.collection.and.processing) |> 
  pivot_longer(names_to = "var_1",values_to = "value_1",-c(sample,Description,Sample.name))

test_df2 <- test |> 
  dplyr::select(sample,Description,Sample.name,PC1:PC9) |>
  pivot_longer(names_to = "var_2",values_to = "value_2",-c(sample,Description,Sample.name))

left_join(test_df1,test_df2,by=c("sample","Description","Sample.name")) |> 
  mutate(comparison = paste0(var_1,"_vs_",var_2)) |> 
  ggplot(aes(x=value_1,y=value_2)) +
  facet_wrap(~comparison,scales = "free",nrow=5) +
  # facet_grid(var_1~var_2,scales = "free_x",drop = T) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/image/panel_metadata_PC.pdf",width = 20,height = 10)


plot_rld <- plotPCA(rld_filter,
                    intgroup = c("sample","Description","Sample name","Donor","Cell Type","Condition","Gender","Passage","Maternal Ethnicity","Paternal Ethnicity","Donor Ethnicity","Source","Approx Time between collection and processing","sizeFactor")) +
  theme_bw()

plot_rld$data %>%
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=Condition)) +
  geom_point(size =3) +
  # ggrepel::geom_text_repel(show.legend = F)+
  theme_bw() + ylab(plot_rld$labels[1]) + xlab(plot_rld$labels[2])
ggsave("../../out/image/PCA_rld.pdf",width = 5,height = 4)

# MSD plot ----------------------------------------------------------------
mds_vsd <- as.data.frame(colData(vds_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_vsd))

ggplot(mds_vsd, aes(x = `1`, y = `2`, color = Condition)) +
  geom_point(size = 3)+ theme_bw()
ggsave("../../out/image/MDS_vsd.pdf",width = 5,height = 4)

mds_rld <- as.data.frame(colData(rld_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_rld))

ggplot(mds_rld, aes(x = `1`, y = `2`, color = Condition)) +
  geom_point(size = 3) + theme_bw()
ggsave("../../out/image/MDS_rld.pdf",width = 5,height = 4)

mdsPois <- as.data.frame(colData(ddsHTSeq_filter)) %>%
  cbind(cmdscale(samplePoisDistMatrix))

ggplot(mdsPois, aes(x = `1`, y = `2`, color = Condition)) +
  geom_point(size = 3) + theme_bw()
ggsave("../../out/image/PoissonDistance_scatter.pdf",width = 5,height = 4)

# save the object of interest ---------------------------------------------
saveRDS(vds_filter,file = "../../out/object/vds_all_filter.rds")
saveRDS(vsd_blind,file = "../../out/object/vsd_all_blind.rds")
saveRDS(rld_filter,file = "../../out/object/rld_all_filter.rds")
saveRDS(rld_blind,file = "../../out/object/rld_all_blind.rds")

# # define the gender of the samples ----------------------------------------
# gene_id <- rownames(assay(vds_filter)) %in% c("XIST","DDX3Y","RPS4Y1","USP9Y")
# 
# mat <- assay(vds_filter)[gene_id, ]
# mat2 <- (mat - rowMeans(mat))/rowSds(mat)
# 
# anno_structure <- as.data.frame(colData(vds_filter)[, c("sample_id","condition")])
# 
# hm_var <- pheatmap(mat2, annotation_col = anno_structure)
# 
# pdf(file = "out/image/heatmap_gender_genes.pdf", width = 6, height = 3.5)
# hm_var
# dev.off()
# define the composition of the 

# biotype of the dataset --------------------------------------------------
# read in the gtf file
gtf_test <- read.table("/home/edo/Documents/reference/gencode.v31.basic.annotation.gtf", header = FALSE, sep = '\t')
# gtf_test <- read.table("data/GCF_000001405.40_GRCh38.p14_genomic.gtf", header = FALSE, sep = '\t')
head(gtf_test)

# from V9 pull the gene_id and the filter only at the gene level
test <- gtf_test %>%
  dplyr::filter(V3 == "gene")

dim(test)
head(test)
# from the colum 9 pull the biotype and gene name
# test2 <- test %>%
#   separate(V9,into = c("gene_id","transcript_id","db_xref", "db_xref2", "description", "gbkey", "gene", "gene_biotype","pseudo"),sep = ";")

gene_ids <- str_extract(test$V9,pattern = "gene_id .*; gene_type") %>%
  str_remove_all("gene_id |; gene_type")

sum(is.na(gene_ids))
head(gene_ids)

gene_name <- str_extract(test$V9,pattern = "gene_name .*; level") %>%
  str_remove_all("gene_name |; level")

sum(is.na(gene_name))
head(gene_name)

biotype <- str_extract(test$V9,pattern = "gene_type \\w+;") %>%
  str_remove_all("gene_type |;")
sum(is.na(biotype))
head(biotype)


LUT_gtf_file <- data.frame(gene_ids,gene_name,biotype)
sum(is.na(LUT_gtf_file))
write_tsv(LUT_gtf_file,"../../out/table/df_LUT_gtf_genes.tsv")
# 
# sort(table(LUT_gtf_file$biotype))
# LUT_gtf_file %>% 
#   dplyr::filter(biotype=="rRNA")
LUT_gtf_file <- read_tsv("../../out/table/df_LUT_gtf_genes.tsv") %>% 
  group_by(gene_name,biotype) %>% 
  summarise() %>% 
  ungroup()

# confirm there are no NAs
sum(is.na(LUT_gtf_file))
# some gene names have multiple annotation
LUT_gtf_file %>% 
  group_by(gene_name) %>% 
  mutate(n = n()) %>% 
  arrange(desc(n))

# add the annotation to the dataset
df_biotype <- data.frame(gene_name = rownames(ddsHTSeq_filter)) %>% 
  left_join(LUT_gtf_file)

df_biotype %>%
  group_by(biotype) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n/sum(n)) %>% 
  # dplyr::filter(!is.na(GENEBIOTYPE)) %>% 
  mutate(biotype = fct_reorder(biotype,-prop)) %>% 
  ggplot(aes(x=biotype,y=prop))+geom_col(position = "dodge")+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90),strip.background = element_blank())
# facet_wrap(~condition)
ggsave("../../out/image/barplot_biotype.pdf",width = 6,height = 4)
