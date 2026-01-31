library(GenomicRanges)
library(rtracklayer)



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

Seqmonk_Gonad_HD$gene_id <- rownames(Seqmonk_Gonad_HD)

Seqmonk_Gonad_F1$gene_id <- rownames(Seqmonk_Gonad_F1)

Seqmonk_Gonad_HD$gene_id <- as.character(Seqmonk_Gonad_HD$gene_id)
tss_window$gene_id <- as.character(tss_window$gene_id)


#Bring the TSS windows to seqmonk, and filter using read counts to those with only greater than 20 reads

TSS_windows_filtered <- read.table(
  "C:/Users/Jessi/OneDrive/Desktop/Annotated Probe Report for TSS.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

