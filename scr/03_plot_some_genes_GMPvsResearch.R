# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(UpSetR)
library(gplots)

# plot genes dotplot ------------------------------------------------------
data <- readRDS("../../out/object/ddsHTSeq_filter_DESeq_GMPvsResearch.rds")

lut <- colData(data) %>%
  data.frame()

# GOI <- read_tsv("../../out/table/res_DIAvsCTRL.txt") %>% 
#   # add a clor variable in case significant
#   mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
#   dplyr::filter(col==1) %>%
#   pull(symbol)

GOI <- c("BMP4","VWF","KDR")


MR <- counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  pivot_longer(names_to = "sample",values_to = "exp",-symbol)%>%
  group_by(sample)%>%
  summarise(MR = sum(exp)/10^6)

# plot the data following the methods implemented in the plotCounts funciton from DESeq2
# Normalized counts plus a pseudocount of 0.5 are shown by default.
counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=Condition,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6)+facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/boxplot_GOI_GMPvsResearch.pdf",width = 9,height = 3)


counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=Condition,y = count_norm_adj,label=Sample.name))+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6)+facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  geom_text_repel()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../../out/image/boxplot_GOI_label.pdf",width = 6,height = 6)
