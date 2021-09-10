## Exercises and solutions for Chapter 7

For this set of exercises, we will use the `chip_1_1.fq.bz2` and `chip_2_1.fq.bz2` files from the `QuasR` package. You can reach the folder that contains the files as follows:
```{r seqProcessEx,eval=FALSE}
folder=(system.file(package="QuasR", "extdata"))
dir(folder) # will show the contents of the folder
```
1. Plot the base quality distributions of the ChIP-seq samples `Rqc` package.
**HINT**: You need to provide a regular expression pattern for extracting the right files from the folder. `"^chip"` matches the files beginning with "chip". [Difficulty: **Beginner/Intermediate**]


**solution:**
```{r,echo=FALSE,eval=FALSE}

folder=(system.file(package="QuasR", "extdata"))
library(Rqc)

# feeds fastq.qz files in "folder" to quality check function
qcRes=rqc(path = folder, pattern = "^chip", openBrowser=FALSE)
rqcCycleQualityBoxPlot(qcRes)

```

2. Now we will trim the reads based on the quality scores. Let's trim 2-4 bases on the 3' end depending on the quality scores. You can use the `QuasR::preprocessReads()` function for this purpose. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
library(QuasR)
# obtain a list of fastq file paths
fastqFiles <- system.file(package="QuasR",
                          "extdata",
                          c("chip_1_1.fq.bz2",
                            "chip_2_1.fq.bz2")
)

# defined processed fastq file names
outfiles <- paste(tempfile(pattern=c("processed_1_",
                                     "processed_2_")),".fq",sep="")

# process fastq files
# trim 2-4 bases on the 3' end depending on the quality scores (we will remove 4 bases)

preprocessReads(fastqFiles, outfiles,
                truncateEndBases=4)
 
```

3. Align the trimmed and untrimmed reads using `QuasR` and plot alignment statistics, did the trimming improve alignments? [Difficulty: **Intermediate/Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

