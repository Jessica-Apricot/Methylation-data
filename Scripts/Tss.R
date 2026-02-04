library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)



gtf <- import("C:/Users/Jessi/OneDrive/Desktop/seqmonk_clean.gtf")

genes <- gtf[gtf$type == "gene"]

# Determine TSS, if + than TSS is at the start of the gene, if - its at the end
tss <- ifelse(strand(genes) == "+", start(genes), end(genes))

# Create 200-bp window
tss_window <- GRanges(
  seqnames = seqnames(genes),
  ranges = IRanges(start = tss - 100, end = tss + 100),
  strand = strand(genes),
  gene_id = genes$gene_id
)

# Export as BED
export(tss_window, "TSS_windows.bed")

#Import my DEG files
Seqmonk_Gonad_HD <- read.table("C:/Users/Jessi/Documents/Microplastics_RNA-seq/Seqmonk_Gonad_HD.tsv")
Seqmonk_Gonad_F1 <- read.table("C:/Users/Jessi/Documents/Microplastics_RNA-seq/Seqmonk_Gonad_F1.tsv")

APR_DSS_DMR <- read.table("C:/Users/Jessi/OneDrive/Desktop/Methylation-data/Annotated probe Reports/Annotated Probe Report for DSS_DMR.txt")

Seqmonk_Gonad_HD$gene_id <- rownames(Seqmonk_Gonad_HD)

Seqmonk_Gonad_F1$gene_id <- rownames(Seqmonk_Gonad_F1)

Seqmonk_Gonad_HD$gene_id <- as.character(Seqmonk_Gonad_HD$gene_id)
tss_window$gene_id <- as.character(tss_window$gene_id)



#This didnt work ignore
#Bring the TSS windows to seqmonk, and filter using read counts to those with only greater than 20 reads

TSS_windows_filtered <- read.table(
  "C:/Users/Jessi/OneDrive/Desktop/Methylation-data/Annotated Probe Reports/Annotated Probe Report for TSS.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


TSS_windows_filtered$mean_meth <- rowMeans(
  TSS_windows_filtered[, c("X07","X08","X09","X10","X11","X12")],
  na.rm = TRUE   # ignore NAs
)
sum(!is.na(TSS_windows_filtered$mean_meth))

summary(TSS_windows_filtered$mean_meth[!is.na(TSS_windows_filtered$mean_meth)])


#Make the TSS filtered table shorter

TSS_filtered_clean <- TSS_windows_filtered %>% 
  select(gene_id = ID, X07, X08, X09, X10, X11, X12)

TSS_filtered_clean_pi<- TSS_filtered_clean %>%
  pivot_longer(
    cols = X07:X12,
    names_to = "sample",
    values_to = "TSS_methylation"
  )

#TSS_DEG_merged_F1 <- TSS_filtered_clean_pi %>%
  left_join(Seqmonk_Gonad_F1, by = 'gene_id')



#TSS_DEG_merged_F0 <-TSS_filtered_clean_pi %>%
  left_join(Seqmonk_Gonad_HD, by = 'gene_id')
