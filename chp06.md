## Exercises and solutions for Chapter 6

The data for the exercises is within the `compGenomRData` package. 

Run the following to see the data files. 
```
dir(system.file("extdata",
             package="compGenomRData"))
```
You will need some of those files to complete the exercises.

### Operations on genomic intervals with the `GenomicRanges` package

1. Create a `GRanges` object using the information in the table below:[Difficulty: **Beginner**]

| chr  | start | end  |strand | score | 
| :--- |:------| :-----| :-----|:-----|
| chr1 | 10000 | 10300 |  +    | 10 |
| chr1 | 11100 | 11500 |  -    | 20 |
| chr2 | 20000 | 20030 |  +    | 15 |


**solution:**
```{r,echo=FALSE,eval=FALSE}
library(GenomicRanges)
gr <- GRanges(seqnames = c("chr1", "chr1", "chr2"),
              ranges = IRanges(start = c(10000,11100,20000),
                               end = c(10300,11500,20030)),
              strand = c("+","-","+"),
              score = c(10,20,15))
gr
 
```


2. Use the `start()`, `end()`, `strand()`,`seqnames()` and `width()` functions on the `GRanges`
object you created. Figure out what they are doing. Can you get a subset of the `GRanges` object for intervals that are only on the + strand? If you can do that, try getting intervals that are on chr1. *HINT:* `GRanges` objects can be subset using the `[ ]` operator, similar to data frames, but you may need
to use `start()`, `end()` and `strand()`,`seqnames()` within the `[]`. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
start(gr) 
end(gr)
seqnames(gr)
strand(gr)
width(gr)
gr[strand(gr)=="+",]
 
```



3. Import mouse (mm9 assembly) CpG islands and RefSeq transcripts for chr12 from the UCSC browser as `GRanges` objects using `rtracklayer` functions. HINT: Check chapter content and modify the code there as necessary. If that somehow does not work, go to the UCSC browser and download it as a BED file. The track name for Refseq genes is "RefSeq Genes" and the table name is "refGene". [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
require(rtracklayer)
session <- browserSession("UCSC",url = 'http://genome-euro.ucsc.edu/cgi-bin/')
genome(session) <- "mm9"
## choose CpG island track on chr12
query_ucsc <- ucscTableQuery(session, track="CpG Islands",table="cpgIslandExt",
                        range=GRangesForUCSCGenome("mm9", "chr12"))
## get the GRanges object for the track
track(query_ucsc)

query_ref <- ucscTableQuery(session, track="RefSeq Genes",table="refGene",
                        range=GRangesForUCSCGenome("mm9", "chr12"))
track(query_ref)
 
```

4. Following from the exercise above, get the promoters of Refseq transcripts (-1000bp and +1000 bp of the TSS) and calculate what percentage of them overlap with CpG islands. HINT: You have to get the promoter coordinates and use the `findOverlaps()` or `subsetByOverlaps()` from the `GenomicRanges` package. To get promoters, type `?promoters` on the R console and see how to use  that function to get promoters or calculate their coordinates as shown in the chapter. [Difficulty: **Beginner/Intermediate**]


**solution:**
```{r,echo=FALSE,eval=FALSE}
prom.refgene <- promoters(track(query_ref), upstream = 1000, downstream = 1000, use.names = FALSE)
prom.refgene

ov = findOverlaps(track(query_ucsc),track(query_ref))

#Calculate the percentage
prom_ov_perc <- (ov@nLnode/ov@nRnode)*100
prom_ov_perc
 
```

5. Plot the distribution of CpG island lengths for CpG islands that overlap with the 
promoters. [Difficulty: **Beginner/Intermediate**]


**solution:**
```{r,echo=FALSE,eval=FALSE}
ov_subset <- subsetByOverlaps(track(query_ucsc),track(query_ref))
ov_length <- width(ov_subset)
hist(ov_length, main="Distribution of overlaped CpG island lengths",
     col = "red",)
d <- density(width(ov_subset))
plot(d, main="Distribution of overlaped CpG island lengths")
polygon(d, col="red")
 
```

6. Get canonical peaks for SP1 (peaks that are in both replicates) on chr21. Peaks for each replicate are located in the `wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1.broadPeak.gz` and `wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep2.broadPeak.gz` files. **HINT**: You need to use `findOverlaps()` or `subsetByOverlaps()` to get the subset of peaks that occur in both replicates (canonical peaks). You can try to read "...broadPeak.gz" files using  the `genomation::readBroadPeak()` function; broadPeak is just an extended BED format. In addition, you can try to use `the coverage()` and `slice()` functions to get more precise canonical peak locations. [Difficulty: **Intermediate/Advanced**]


**solution:**
```{r,echo=FALSE,eval=FALSE}
library(genomation)
filePathPeaks=system.file("extdata",   
              "wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1.broadPeak.gz",
                      package="compGenomRData")

filePathPeaks2=system.file("extdata",   
              "wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep2.broadPeak.gz",
                      package="compGenomRData")
# read the peaks from a bed file
pk1.gr=readBroadPeak(filePathPeaks)
pk2.gr=readBroadPeak(filePathPeaks2)

# find the canonical peaks in both replicates
pk.gr = subsetByOverlaps(pk1.gr, pk2.gr)

# only chr21 peaks
pk.gr.chr21 <- pk.gr[seqnames(pk.gr)=="chr21",]
pk.gr.chr21
 
```

### Dealing with mapped high-throughput sequencing reads

1.  Count the reads overlapping with canonical SP1 peaks using the BAM file for one of the replicates. The following file in the `compGenomRData` package contains the alignments for SP1 ChIP-seq reads: `wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam`. **HINT**: Use functions from the `GenomicAlignments` package. [Difficulty: **Beginner/Intermediate**]


**solution:**
```{r,echo=FALSE,eval=FALSE}
library(Rsamtools)
library(GenomicAlignments)
bamfilePath=system.file("extdata",
            "wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam",
                      package="compGenomRData")
# get reads for regions of interest from the bam file
pk1.gr.chr21 <- pk1.gr[seqnames(pk1.gr)=="chr21",]
param <- ScanBamParam(which=pk1.gr.chr21)
counts=countBam(bamfilePath, param=param)
 
```

### Dealing with contiguous scores over the genome

1. Extract the `Views` object for the promoters on chr20 from the `H1.ESC.H3K4me1.chr20.bw` file available at `CompGenomRData` package. Plot the first "View" as a line plot. **HINT**: See the code in the relevant section in the chapter and adapt the code from there. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

2. Make a histogram of the maximum signal for the Views in the object you extracted above. You can use any of the view summary functions or use `lapply()` and write your own summary function. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Get the genomic positions of maximum signal in each view and make a `GRanges` object. **HINT**: See the `?viewRangeMaxs` help page. Try to make a `GRanges` object out of the returned object. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

### Visualizing and summarizing genomic intervals

1. Extract -500,+500 bp regions around the TSSes on chr21; there are refseq files for the hg19 human genome assembly in the `compGenomRData` package. Use SP1 ChIP-seq data in the `compGenomRData` package, access the file path via the `system.file()` function, the file name is:
`wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam`. Create an average profile of read coverage around the TSSes. Following that, visualize the read coverage with a heatmap. **HINT**: All of these are possible using the `genomation` package functions. Check `help(ScoreMatrix)` to see how you can use bam files. As an example here is how you can get the file path to refseq annotation on chr21. [Difficulty: **Intermediate/Advanced**]
```{r example,eval=FALSE}
transcriptFilechr21=system.file("extdata",
                      "refseq.hg19.chr21.bed",
                      package="compGenomRData")
```

**solution:**
```{r}
# get transcription start sites on chr21
transcriptFilechr21=system.file("extdata",
                      "refseq.hg19.chr21.bed",
                      package="compGenomRData")
feat=readTranscriptFeatures(transcriptFilechr21,
                            remove.unusual = TRUE,
                            up.flank = 500, down.flank = 500)
prom=feat$promoters # get promoters from the features


bamfilePathchr21=system.file("extdata",
                      "wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam",
                      package="compGenomRData")



sm=ScoreMatrix(bamfilePathchr21,prom,
               type="bam",strand.aware = TRUE)


# look for the average coverage
plotMeta(sm, profile.names = "SP1 ChIP-seq", xcoords = c(-500,500),
         ylab="SP1 ChIP-seq coverage",dispersion = "se",
         xlab="bases around TSS")

heatMatrix(sm,order=TRUE,xcoords = c(-500,500),
           xlab="bases around TSS")

```

2. Extract -500,+500 bp regions around the TSSes on chr20. Use H3K4me3 (`H1.ESC.H3K4me3.chr20.bw`) and H3K27ac (`H1.ESC.H3K27ac.chr20.bw`) ChIP-seq enrichment data in the `compGenomRData` package and create heatmaps and average signal profiles for regions around the TSSes.[Difficulty: **Intermediate/Advanced**]

**solution:**
```{r}
# get transcription start sites on chr20
transcriptFilechr20=system.file("extdata",
                      "refseq.hg19.chr20.bed",
                      package="compGenomRData")
feat=readTranscriptFeatures(transcriptFilechr20,
                            remove.unusual = TRUE,
                            up.flank = 500, down.flank = 500)
prom=feat$promoters # get promoters from the features

# get for H3K4me3 values around TSSes
# we use strand.aware=TRUE so - strands will
# be reversed
H3K4me3File=system.file("extdata",
                      "H1.ESC.H3K4me3.chr20.bw",
                      package="compGenomRData")


H3K27acFile=system.file("extdata",
                      "H1.ESC.H3K27ac.chr20.bw",
                      package="compGenomRData")

sml=ScoreMatrixList(c(H3K4me3=H3K4me3File,
                      H3K27ac=H3K27acFile),prom,
                      type="bigWig",strand.aware = TRUE)



plotMeta(sml,profile.names = c("H3K4me3", "H3K27ac"), xcoords = c(-500,500),
         ylab="H3K4me3 enrichment",dispersion = "se",
         xlab="bases around TSS")

set.seed(1029)
multiHeatMatrix(sml,order=TRUE,xcoords = c(-500,500),
                xlab="bases around TSS",winsorize = c(0,95),
                matrix.main = c("H3K4me3","H3K27ac"),
                column.scale=TRUE,
                clustfun=function(x) kmeans(x, centers=3)$cluster)
```

3. Download P300 ChIP-seq peaks data from the UCSC browser. The peaks are locations where P300 binds. The P300 binding marks enhancer regions in the genome. (**HINT**:  group: "regulation", track: "Txn Factor ChIP", table:"wgEncodeRegTfbsClusteredV3", you need to filter the rows for "EP300" name.) Check enrichment of H3K4me3, H3K27ac and DNase-seq (`H1.ESC.dnase.chr20.bw`) experiments on chr20 on and arounf the P300 binding-sites, use data from `compGenomRData` package. Make multi-heatmaps and metaplots. What is different from the TSS profiles? [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Cluster the rows of multi-heatmaps for the task above. Are there obvious clusters? **HINT**: Check arguments of the `multiHeatMatrix()` function. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


5. Visualize one of the -500,+500 bp regions around the TSS using `Gviz` functions. You should visualize both H3K4me3 and H3K27ac and the gene models. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

