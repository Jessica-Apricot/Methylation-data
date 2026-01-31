#BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

#BiocManager::install(c("DSS"))
library(DSS)
library(bsseq)



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

write.table(dmrs, 
            file = "DSS_DMR",
            sep =  "\t",
            quote = FALSE,
            row.names = FALSE)

control <- c(82.3597, 82.8544, 82.4065)
treatment <- c(82.5454, 81.8694, 82.1792)

t.test(control, treatment)

