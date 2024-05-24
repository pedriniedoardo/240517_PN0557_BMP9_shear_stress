# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)

# read in the data --------------------------------------------------------
# dds values
ddsHTSeq_filter <- readRDS("../../out/object/dds_all_filter.rds") %>% 
  DESeq(.)
# scaled values
vds_filter <- readRDS("../../out/object/vds_all_filter.rds")

# extract expression values
df_norm <- counts(ddsHTSeq_filter,normalized = T) |> 
  data.frame() |> 
  rownames_to_column("gene")
df_raw <- counts(ddsHTSeq_filter,normalized = F) |> 
  data.frame() |> 
  rownames_to_column("gene")
df_scaled <- assay(vds_filter) |> 
  data.frame() |> 
  rownames_to_column("gene")

# save tables -------------------------------------------------------------
write_tsv(df_norm,"../../out/table/read_norm_all.tsv")
write_tsv(df_raw,"../../out/table/read_raw_all.tsv")
write_tsv(df_scaled,"../../out/table/read_scaled_all.tsv")
