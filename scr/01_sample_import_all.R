# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")

# read in the data --------------------------------------------------------
# # read in all the files
# folder <- "data/"
# file <- dir(folder) %>% 
#   str_subset(pattern = "all.counts")
# 
# sample_id <- file %>% 
#   str_extract_all(pattern = "PN0401_\\d+") %>% 
#   unlist()

# sample read one file
Cts <-"../../data/all.counts"
test <- read.delim(Cts,header = T)
head(test)
# colnames(test) <- c("Geneid","Chr","Start","End","Strand","Length","Count")

# # in case of fractional count I should round them to the closest integer https://support.bioconductor.org/p/75848/
# test_fix <- test %>% 
#   mutate(Count_rounded = round(Count))

# df_out <- test_fix %>% 
#   dplyr::select(Geneid,Count_rounded)

# are there any replicates in the geneid ?
# df_out %>% 
#   group_by(Geneid) %>% 
#   summarise(n = n()) %>% 
#   arrange(desc(n))

# test_fix[1:10,]

# # read in all the files
# list_reads <- lapply(sample_id,function(x){
#   Cts <-paste0("../data/05_featureCounts_out_50/",x,"_featureCounts_Multimapping_fraction_output.txt")
#   test <- read.delim(Cts,skip = 1)
#   colnames(test) <- c("Geneid","Chr","Start","End","Strand","Length","Count")
#   # in case of fractional count I should round them to the closest integer https://support.bioconductor.org/p/75848/
#   test_fix <- test %>% 
#     mutate(Count_rounded = round(Count))
#   
#   df_out <- test_fix %>% 
#     dplyr::select(Geneid,Count_rounded)
#   
#   return(df_out)
# }) %>% 
#   setNames(sample_id)

# confirm the dimension of all the files is the same
# lapply(list_reads,function(x){dim(x)})
dim(test)

# confirm no dataset have replicated gene names
# lapply(list_reads,function(x){x %>% group_by(Geneid) %>% 
#     summarise(n = n()) %>% 
#     arrange(desc(n))
# })
table(test$Geneid) %>% 
  sort(decreasing = T) %>% 
  head()

# # make the list as a single matrix
# df_exp <- pmap(list(list_reads,names(list_reads)),function(x,name){
#   x %>%
#     dplyr::rename(!!name := "Count_rounded")
# }) %>% 
#   purrr::reduce(left_join,by=c("Geneid"))

# # save the table
# df_exp %>% 
#   write_tsv("out/table/featureCounts_Multimapping_fraction_output_ALL_raw.tsv")
# 
# S <- df_exp %>% 
#   column_to_rownames("Geneid") %>% 
#   as.matrix()

# build the annotation besed on the sample metadata
LUT_samples <- read_csv("../../data/lut_sample.csv")
  # remove sample c47 from the analysis becousa of same issue with the gender assignament
  # remove a sample that has some issue with the gender assignament
  # dplyr::filter(!(clone %in% c("c03")))
# separate(Description,into = c("clone","treat","sirna"),sep = " ",remove = F)

# notice there are three samples more (PN0433) for which I don't know what they are
# extract only the count information
mat_exp <- test %>% 
  dplyr::select(-c("Chr","Start","End","Strand","Length")) %>% 
  column_to_rownames("Geneid") %>% 
  dplyr::select(contains("PN0453")) %>% 
  # keep all the samples of interest
  dplyr::select(LUT_samples$Sample_ID) %>% 
  as.matrix()

# match the order of the sample in the matrix with the sample in the sample sheet
# build the metadata for the deseq object
coldata <- data.frame(sample = colnames(mat_exp)) %>% 
  left_join(LUT_samples,by = c("sample"="Sample_ID")) %>% 
  mutate(rowname = sample) %>% 
  column_to_rownames("rowname")

# save the table
write.csv(coldata,file = "../../data/LUT_samples_final.csv",row.names = T)

# It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.
# in this case the coldata has been build from the expression matrtix, therefore they have to be in order
all(rownames(coldata) %in% colnames(mat_exp))
all(rownames(coldata) == colnames(mat_exp))

# define the model --------------------------------------------------------
# cannot block by clone as the treatment is confounded by the clone
# clone <- factor(coldata$clone)
# block by gender. use the new gender variable
# clone <- factor(coldata$Donor)
condition <- factor(coldata$Condition,levels = c("HUVEC","Research","GMP"))

# build the design
# design <- model.matrix(~ clone+condition)
design <- model.matrix(~ condition)
# colnames(design) <- c("intercept","clonec02","clonec03","clonec04","clonec05","conditionHypoxia")
colnames(design)[1] <- c("intercept")

# I am not assuming the rep variable has to be blocked
dds <- DESeqDataSetFromMatrix(countData = mat_exp,
                              colData = coldata,
                              design = design)

# save the raw object
saveRDS(object = dds,file = "../../out/object/dds_all_raw.rds")
saveRDS(object = design,file = "../../out/object/design_all_dds.rds")
