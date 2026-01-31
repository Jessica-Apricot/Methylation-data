#BiocManager::install("methylKit")
library("methylKit")


#Pathways to Methylation .cov files
file.list <- 
c("C:/Users/Jessi/Documents/Methylation_data/Extract/07/AAALMCTHV-6434-07-0-1_S7_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz", 
               "C:/Users/Jessi/Documents/Methylation_data/Extract/08/AAALMCTHV-6434-08-0-1_S8_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz", 
               "C:/Users/Jessi/Documents/Methylation_data/Extract/09/AAALMCTHV-6434-09-0-1_S9_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz", 
               "C:/Users/Jessi/Documents/Methylation_data/Extract/10/AAALMCTHV-6434-10-0-1_S10_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz", 
               "C:/Users/Jessi/Documents/Methylation_data/Extract/11/AAALMCTHV-6434-11-0-1_S11_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz", 
               "C:/Users/Jessi/Documents/Methylation_data/Extract/12/AAALMCTHV-6434-12-0-1_S12_R1_001_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
)


treatment <- c(0,0,0,1,1,1)

sample.id <- c("C1","C2","C3","D1","D2","D3")

file.list <- normalizePath(file.list)

# Reading each file on its own because I was getting an annoying error about receiving a single location
Meth_list <- lapply(1:length(file.list), function(i) {
  methRead(
    location = file.list[i],
    sample.id = sample.id[i],
    assembly = "GRCz11",
    treatment = treatment[i],
    context = "CpG",
    pipeline = "bismarkCoverage"
  )
})

# Convert the list into a methylRawList
Meth_list_obj <- new(
  "methylRawList",
  Meth_list,
  treatment = treatment
)


#Make into single object
Meth_files <- unite(Meth_list_obj)


# where only CpGs covered with at least 1 sample per group will be returned - retains more cpg
meth_min <- unite(Meth_list_obj, min.per.group=1L)


#Methylation statistics for each sample, make histogram
#Controls
#07
getMethylationStats(Meth_list[[1]],plot=FALSE,both.strands=FALSE)

getMethylationStats(Meth_list[[1]],plot=TRUE,both.strands=FALSE)

getCoverageStats(Meth_list[[1]],plot=TRUE,both.strands=FALSE)

#08
getMethylationStats(Meth_list[[2]],plot=FALSE,both.strands=FALSE)

getMethylationStats(Meth_list[[2]],plot=TRUE,both.strands=FALSE)

getCoverageStats(Meth_list[[2]],plot=TRUE,both.strands=FALSE)

#09
getMethylationStats(Meth_list[[3]],plot=FALSE,both.strands=FALSE)

getMethylationStats(Meth_list[[3]],plot=TRUE,both.strands=FALSE)

getCoverageStats(Meth_list[[3]],plot=TRUE,both.strands=FALSE)

#Treatment
#10
getMethylationStats(Meth_list[[4]],plot=FALSE,both.strands=FALSE)

getMethylationStats(Meth_list[[4]],plot=TRUE,both.strands=FALSE)

getCoverageStats(Meth_list[[4]],plot=TRUE,both.strands=FALSE)

#11
getMethylationStats(Meth_list[[5]],plot=FALSE,both.strands=FALSE)

getMethylationStats(Meth_list[[5]],plot=TRUE,both.strands=FALSE)

getCoverageStats(Meth_list[[5]],plot=TRUE,both.strands=FALSE)

#12
getMethylationStats(Meth_list[[6]],plot=FALSE,both.strands=FALSE)

getMethylationStats(Meth_list[[6]],plot=TRUE,both.strands=FALSE)

getCoverageStats(Meth_list[[6]],plot=TRUE,both.strands=FALSE)

#Finding correlation
getCorrelation(Meth_files,plot=TRUE)
#Make a dendrogram
clusterSamples(Meth_files, dist="correlation", method="ward", plot=TRUE)

PCASamples(Meth_files)

Diff_Meth <- calculateDiffMeth(Meth_files)

#Diff methylation in the less stringent
Diff_Meth_LS <- calculateDiffMeth(meth_min)

# get hyper methylated bases - percent difference greater than 25% - More stringent
Hyper_Diff_Meth25 <- getMethylDiff(Diff_Meth, difference=25, qvalue=0.01, type="hyper")

# get hypo methylated bases - percent difference greater than 25% - More stringent
Hypo_Diff_Meth25 <- getMethylDiff(Diff_Meth, difference=25, qvalue=0.01, type="hypo")

#Lowering the percentage difference
# get hyper methylated bases - percent difference greater than 25% - More stringent
Hyper_Diff_Meth20 <- getMethylDiff(Diff_Meth, difference=20, qvalue=0.01, type="hyper")

# get hypo methylated bases - percent difference greater than 25% - More stringent
Hypo_Diff_Meth20 <- getMethylDiff(Diff_Meth, difference=20, qvalue=0.01, type="hypo")


diffMethPerChr(Diff_Meth, plot=FALSE, qvalue.cutoff=0.01, meth.cutoff=20)



#Using the less stringent (the bases dont have to be present in all samples)
Hyper_Diff_Meth25_LS <- getMethylDiff(Diff_Meth_LS, difference=25, qvalue=0.01, type="hyper")

Hypo_Diff_Meth25_LS <- getMethylDiff(Diff_Meth_LS, difference=25, qvalue=0.01, type="hypo")


## Annotating

library(GenomicRanges)

df07 <- getData(Meth_list[[1]])

head(df07)

# read CGI file
cpg_islands <- read.table("C:/Users/Jessi/OneDrive/Desktop/cpg_islands.bed", header=FALSE)

colnames(cpg_islands) <- c("chr","start","end","strand")
cpg_gr <- GRanges(seqnames=cpg_islands$chr,
                  ranges=IRanges(start=cpg_islands$start,
                                 end=cpg_islands$end))

gr07 <- GRanges(seqnames = df07$chr,
              ranges = IRanges(start = df07$start, end = df07$end),
              coverage = df07$coverage)

overlaps <- findOverlaps(cpg_gr, gr07)

island_cov <- tapply(gr07$coverage[subjectHits(overlaps)],
                     queryHits(overlaps), mean)

mean(island_cov, na.rm = TRUE)
