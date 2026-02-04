#BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

#BiocManager::install(c("DSS"))
library(DSS)
library(bsseq)
library(GenomicRanges)
library(dplyr)
library(openxlsx)


#Use file.list 
Files <- lapply(file.list, function(f)
  read.table(gzfile(f), header = FALSE)
)


sample07 <- read.table(gzfile("C:/Users/Jessi/Documents/Methylation_data/Extract/07/AAALMCTHV-6434-07-0-1_S7_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"), header = FALSE)

sample08 <- read.table(gzfile("C:/Users/Jessi/Documents/Methylation_data/Extract/08/AAALMCTHV-6434-08-0-1_S8_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"), header = FALSE)

sample09 <- read.table(gzfile("C:/Users/Jessi/Documents/Methylation_data/Extract/09/AAALMCTHV-6434-09-0-1_S9_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"), header = FALSE)

sample10 <- read.table(gzfile("C:/Users/Jessi/Documents/Methylation_data/Extract/10/AAALMCTHV-6434-10-0-1_S10_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"), header = FALSE)

sample11 <- read.table(gzfile("C:/Users/Jessi/Documents/Methylation_data/Extract/11/AAALMCTHV-6434-11-0-1_S11_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"), header = FALSE)

sample12 <- read.table(gzfile("C:/Users/Jessi/Documents/Methylation_data/Extract/12/AAALMCTHV-6434-12-0-1_S12_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"), header = FALSE)

dss07 <- data.frame(
    chr = sample07[,1],
    pos = sample07[,2],
    N   = sample07[,5] + sample07[,6],
    X   = sample07[,5]
  )

rm(sample07)

dss08 <- data.frame(
  chr = sample08[,1],
  pos = sample08[,2],
  N   = sample08[,5] + sample08[,6],
  X   = sample08[,5]
)

rm(sample08)

dss09 <- data.frame(
  chr = sample09[,1],
  pos = sample09[,2],
  N   = sample09[,5] + sample09[,6],
  X   = sample09[,5]
)

rm(sample09)

dss10 <- data.frame(
  chr = sample10[,1],
  pos = sample10[,2],
  N   = sample10[,5] + sample10[,6],
  X   = sample10[,5]
)

rm(sample10)

dss11 <- data.frame(
  chr = sample11[,1],
  pos = sample11[,2],
  N   = sample11[,5] + sample11[,6],
  X   = sample11[,5]
)

rm(sample11)

dss12 <- data.frame(
  chr = sample12[,1],
  pos = sample12[,2],
  N   = sample12[,5] + sample12[,6],
  X   = sample12[,5]
)

rm(sample12)



dss_list <- list(dss07, dss08, dss09, dss10, dss11, dss12)

dss_list <- lapply(dss_list , function(x) {
  x[x$N >=5, ]
})


BSobj <- makeBSseqData(
  dss_list,
  sampleNames = c("C1","C2","C3","D1", "D2", "D3")
)

group <- c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")

dmlTest <- DMLtest(
  BSobj,
  group1 = which(group == "Control"),
  group2 = which(group == "Treatment"),
  smoothing = TRUE,
  smoothing.span = 500
)

dmrs <- callDMR(
  dmlTest, 
  p.threshold = 0.00001, 
  delta = 0,
  minlen= 50,
  minCG = 5
)

dmr25 <- dmrs[abs(dmrs$diff.Methy) >= 0.25, ]


write.table(dmrs, 
            file = "DSS_DMR",
            sep =  "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(dmr25,
            file = "DSS_DMR_25",
            sep =  "\t",
            quote = FALSE,
            row.names = FALSE)

APR_DSS_DMR <- read.table("C:/Users/Jessi/OneDrive/Desktop/Methylation-data/Annotated probe Reports/Annotated Probe Report for DSS_DMR.txt",
          header = TRUE,
          sep = "\t",
          quote = "",
          fill = TRUE,
          stringsAsFactors = FALSE)

colnames(APR_DSS_DMR)
colnames(Seqmonk_Gonad_F1)
colnames(dmrs)

#Creating tables with Gonad DEGs and DSS DMRs

F1_common_genes <- intersect(APR_DSS_DMR$Feature, Seqmonk_Gonad_F1$gene_id)
F0_HD_common_genes <-intersect(APR_DSS_DMR$Feature, Seqmonk_Gonad_HD$gene_id)

#Making a GRanges object for the annotated probe report of the DSS DMR (which associate the DMRs with the nearest gene)
annot_gr <- GRanges(
  seqnames = APR_DSS_DMR$Chr,
  ranges = IRanges(
    start = APR_DSS_DMR$Start,
    end   = APR_DSS_DMR$End
  ),
  Gene = APR_DSS_DMR$Feature
)

#Make a GRanges of the original DSS object, as making probe reports excludes some information we want, i.e diff.methyl
dmrs_gr <- GRanges(
  seqnames = dmrs$chr,
  ranges = IRanges(
    start = dmrs$start,
    end   = dmrs$end
  ),
  diff.methyl = dmrs$diff.Methy
)

#Applying the missing info from the original object to the probe report
Overlap <- findOverlaps(annot_gr, dmrs_gr)

DMR_full <- APR_DSS_DMR[queryHits(Overlap), ] %>%
  mutate(
    diff.Methyl = dmrs$diff.Methy[subjectHits(Overlap)]
  )

#Subsetting the full probe report + DSS object, to get only the columns we want
DMR_Full_cols <- c(
  "Chromosome",
  "Start",
  "End",
  "Feature",
  "Distance",
  "diff.Methyl"
)

DSS_DMR_subset <- DMR_full[, DMR_Full_cols]

#The joined genes and DMR table
Gonad_HD_DSS <- Seqmonk_Gonad_HD %>%
  inner_join(DSS_DMR_subset, by = c("gene_id" = "Feature"))

#Some of the genes in this table have repeats, so make a table to check how many have multiple So to get the total number of genes 
No_repeats_Gonad_HD_DSS <- Gonad_HD_DSS %>%
  group_by(gene_id) %>%
  slice_max(order_by = abs(diff.Methyl), n = 1) %>%
  ungroup()

#Write the table (with repeat genes) into excel 
write.xlsx(
 Gonad_HD_DSS,
  file = "Gonad_HD_DSS.xlsx",
  rowNames = TRUE
)

#Does the methylation make biological sense, with the current knowledge about how methylation relates to the expresssion level
Methylation_comp <- Gonad_HD_DSS %>%
  mutate(
    Meth_dir = case_when(
      diff.Methyl < 0 ~ "Hypomethylated",
      diff.Methyl > 0 ~ "Hypermethylated",
      TRUE ~ "No_change"
    )
  )

table(Methylation_comp$ExpressionGroup, Methylation_comp$Meth_dir)

### Lets try with the other samples ###


#Brain F1
Brain_F1 <- read.csv("C:/Users/Jessi/OneDrive/Desktop/Methylation-data/Gene tables/Br_res_s_F1.csv")

colnames(Brain_F1)

Brain_F1_cols <- c(
  "X",
  "log2FoldChange",
  "padj"
)

Brain_F1_subset <- Brain_F1[, Brain_F1_cols]

Brain_F1_DSS <- Brain_F1_subset %>%
  inner_join(DSS_DMR_subset, by = c("X" = "Feature"))

#Write the table into excel 
write.xlsx(
  Brain_F1_DSS,
  file = "Meth_Gene Tables/Brain_F1_DSS.xlsx",
  rowNames = TRUE
)

#Brain HD
Brain_F0_HD <- read.csv("C:/Users/Jessi/OneDrive/Desktop/Methylation-data/Gene tables/Br_res_s_HD.csv")


gene_cols <- c(
  "X",
  "log2FoldChange",
  "padj"
)

Brain_F0_HD_subset <- Brain_F0_HD[, gene_cols]

#0 genes!!! 
Brain_F0_HD_DSS <- Brain_F0_HD_subset %>%
  inner_join(DSS_DMR_subset, by = c("X" = "Feature"))

#Gonad HL
Gonad_F0_HL <- read.csv("C:/Users/Jessi/OneDrive/Desktop/Methylation-data/Gene tables/Go_res_s_HL.csv")


Gonad_F0_HL_subset <- Gonad_F0_HL[, gene_cols]

# 
Gonad_F0_HL_DSS <- Gonad_F0_HL_subset %>%
  inner_join(DSS_DMR_subset, by = c("X" = "Feature"))

#Write the table into excel 
write.xlsx(
  Gonad_F0_HL_DSS,
  file = "Meth_Gene Tables/Gonad_F0_HL_DSS.xlsx",
  rowNames = TRUE
)

#Liver F1
Liver_F1 <- read.csv("C:/Users/Jessi/OneDrive/Desktop/Methylation-data/Gene tables/Li_res_s_F1.csv")

Liver_F1_subset <- Liver_F1[, gene_cols]

# 
Liver_F1_DSS <- Liver_F1_subset %>%
  inner_join(DSS_DMR_subset, by = c("X" = "Feature"))

#Write the table into excel 
write.xlsx(
  Liver_F1_DSS,
  file = "Meth_Gene Tables/Liver_F1_DSS.xlsx",
  rowNames = TRUE
)

#Liver LD
Liver_F0_LD <- read.csv("C:/Users/Jessi/OneDrive/Desktop/Methylation-data/Gene tables/Li_res_s_LD.csv")

Liver_F0_LD_subset <- Liver_F0_LD[, gene_cols]

Liver_F0_LD_DSS <- Liver_F0_LD_subset %>%
  inner_join(DSS_DMR_subset, by = c("X" = "Feature"))

#Write the table into excel 
write.xlsx(
  Liver_F1_DSS,
  file = "Meth_Gene Tables/Liver_F0_LD_DSS.xlsx",
  rowNames = TRUE
)
