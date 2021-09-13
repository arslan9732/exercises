## Exercises and solutions for Chapter 8

### Exploring the count tables

Here, import an example count table and do some exploration of the expression data. 

```{r exSetup1}
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                            package = "compGenomRData")
```

1. Normalize the counts using the TPM approach. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
geneLengths <- counts[,11]
tpm <- apply(counts[,1:10], 
             2,
             function(x) {
                     10^6 * (((x)/geneLengths)/(sum(x/geneLengths)))
             })

# Second method to count TPM
rpk <- apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000))
#normalize by the sample size using rpk values
tpm1 <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
 
```

2. Plot a heatmap of the top 500 most variable genes. Compare with the heatmap obtained using the 100 most variable genes. [Difficulty: **Beginner**]

**solution:** The clustering of genes in top 500 variable is vice verse as compared to top 100 variable genes.
```{r}
V <- apply(tpm, 1, var)
#sort the results by variance in decreasing order 
#and select the top 500 and 100 genes 
selectedGenes500 <- names(V[order(V, decreasing = T)][1:500])
selectedGenes100 <- names(V[order(V, decreasing = T)][1:100])
library(pheatmap)


pheatmap(tpm[selectedGenes500,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData)


pheatmap(tpm[selectedGenes100,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData)
```

3. Re-do the heatmaps setting the `scale` argument to `none`, and `column`. Compare the results with `scale = 'row'`. [Difficulty: **Beginner**]

**solution:**
```{r}

pheatmap(tpm[selectedGenes500,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData,
         main = "Heatmap of top 500 variable genes woth 'row' scale")
pheatmap(tpm[selectedGenes500,], scale = 'none', 
         show_rownames = FALSE, 
         annotation_col = colData,
         main = "Heatmap of top 500 variable genes with 'none' scale")

pheatmap(tpm[selectedGenes100,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData,
         main = "Heatmap of top 100 variable genes with 'row' scale")
pheatmap(tpm[selectedGenes100,], scale = 'none', 
         show_rownames = FALSE, 
         annotation_col = colData,
         main = "Heatmap of top 100 variable genes with 'none' scale")
```

4. Draw a correlation plot for the samples depicting the sample differences as 'ellipses', drawing only the upper end of the matrix, and order samples by hierarchical clustering results based on `average` linkage clustering method. [Difficulty: **Beginner**]

**solution:**
```{r}
# correlation plot

library(stats)
correlationMatrix <- cor(tpm)

library(corrplot)
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7, hclust.method = 'average',
         method = 'ellipse', type = 'upper') 

```
5. How else could the count matrix be subsetted to obtain quick and accurate clusters? Try selecting the top 100 genes that have the highest total expression in all samples and re-draw the cluster heatmaps and PCA plots. [Difficulty: **Intermediate**]

**solution:**
```{r}
Only.Counts <- counts[,1:10]
S <- apply(Only.Counts, 1, sum)
#sort the results by total increasing order 
#and select the top 100 genes 
selectedTopGenes <- names(S[order(S, decreasing = T)][1:100])

pheatmap(Only.Counts[selectedTopGenes,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData,
         main = "Heatmap of top 100 highly expressed genes")

library(stats)
library(ggplot2)
library(ggfortify)
#transpose the matrix 
M <- t(Only.Counts[selectedTopGenes,])
# transform the counts to log2 scale 
M <- log2(M + 1)
#compute PCA 
pcaResults <- prcomp(M)

#plot PCA results making use of ggplot2's autoplot function
#ggfortify is needed to let ggplot2 know about PCA data structure. 
autoplot(pcaResults, data = colData, colour = 'group')

```

6. Add an additional column to the annotation data.frame object to annotate the samples and use the updated annotation data.frame to plot the heatmaps. (Hint: Assign different batch values to CASE and CTRL samples). Make a PCA plot and color samples by the added variable (e.g. batch). [Difficulty: Intermediate]

**solution:**
```{r}
annot <- colData
annot$batch <- c("cDNA", "cDNA", "size fractionation",
                 "size fractionation","cDNA", "cDNA",
                 "size fractionation","size fractionation",
                 "cDNA", "cDNA")
pheatmap(Only.Counts[selectedTopGenes,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = annot,
         main = "Heatmap of top 100 highly expressed genes")

M <- t(Only.Counts[selectedTopGenes,])
# transform the counts to log2 scale 
M <- log2(M + 1)
#compute PCA 
pcaResults <- prcomp(M)

#plot PCA results making use of ggplot2's autoplot function
#ggfortify is needed to let ggplot2 know about PCA data structure. 
autoplot(pcaResults, data = annot, colour = 'batch')
```

7. Try making the heatmaps using all the genes in the count table, rather than sub-selecting. [Difficulty: **Advanced**]

**solution:**
```{r}
pheatmap(Only.Counts, scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData,
         main = "Heatmap of All genes")
#Error in hclust(d, method = method) : NA/NaN/Inf in foreign function call (arg 10)

# Due to above error, we have to remove all the rows which have 0 values across the samples

SelectedGenes <- Only.Counts[ rowSums(Only.Counts) > 0, ]

pheatmap(SelectedGenes, scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData,
         main = "Heatmap of All genes")

```

8. Use the [`Rtsne` package](https://cran.r-project.org/web/packages/Rtsne/Rtsne.pdf) to draw a t-SNE plot of the expression values. Color the points by sample group. Compare the results with the PCA plots. [Difficulty: **Advanced**]

**solution:**
```{r}
library(Rtsne)
library(ggplot2)
library(ggpubr)

set.seed(42) # Set a seed if you want reproducible results

#transpose the matrix 
M <- t(Only.Counts)
# transform the counts to log2 scale 
M <- log2(M + 1)
#compute PCA 

pcaResults <- prcomp(M)

#plot PCA results making use of ggplot2's autoplot function
#ggfortify is needed to let ggplot2 know about PCA data structure. 
PCA = autoplot(pcaResults, data = colData, colour = 'group')
#compute t-SNE
tsne_out <- Rtsne(M,perplexity = 3)

tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2],
col = colData$group)


tsne = ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) +
theme_classic() + ggtitle("t-SNE plot perplexity = 3") +
theme(plot.title = element_text(hjust = 0.5))

PCA
tsne

```


### Differential expression analysis

Firstly, carry out a differential expression analysis starting from raw counts.
Use the following datasets:

```
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv", 
                            package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                            package = "compGenomRData")
```

- Import the read counts and colData tables.
- Set up a DESeqDataSet object.
- Filter out genes with low counts.
- Run DESeq2 contrasting the `CASE` sample with `CONTROL` samples. 

```{r, DEG testing}
library(DESeq2)
counts <- read.table(counts_file)
counts <- subset(counts, select = c(-width))
colData <- read.table(coldata_file, header = T, sep = "\t")
designformula <- '~ group'

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = as.formula(designformula))
#and remove those that don't have at least 1 read. 
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]

dds <- DESeq(dds)
#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))
#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]
```

Now, you are ready to do the following exercises: 

1. Make a volcano plot using the differential expression analysis results. (Hint: x-axis denotes the log2FoldChange and the y-axis represents the -log10(pvalue)). [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

2. Use DESeq2::plotDispEsts to make a dispersion plot and find out the meaning of this plot. (Hint: Type ?DESeq2::plotDispEsts) [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Explore `lfcThreshold` argument of the `DESeq2::results` function. What is its default value? What does it mean to change the default value to, for instance, `1`? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. What is independent filtering? What happens if we don't use it? Google `independent filtering statquest` and watch the online video about independent filtering. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

5. Re-do the differential expression analysis using the `edgeR` package. Find out how much DESeq2 and edgeR agree on the list of differentially expressed genes. [Difficulty: **Advanced**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

6. Use the `compcodeR` package to run the differential expression analysis using at least three different tools and compare and contrast the results following the `compcodeR` vignette. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


### Functional enrichment analysis

1. Re-run gProfileR, this time using pathway annotations such as KEGG, REACTOME, and protein complex databases such as CORUM, in addition to the GO terms. Sort the resulting tables by columns `precision` and/or `recall`. How do the top GO terms change when sorted for `precision`, `recall`, or `p.value`? [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

2. Repeat the gene set enrichment analysis by trying different options for the `compare` argument of the `GAGE:gage`
function. How do the results differ? [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Make a scatter plot of GO term sizes and obtained p-values by setting the `gProfiler::gprofiler` argument `significant = FALSE`. Is there a correlation of term sizes and p-values? (Hint: Take -log10 of p-values). If so, how can this bias be mitigated? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Do a gene-set enrichment analysis using gene sets from top 10 GO terms. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

5. What are the other available R packages that can carry out gene set enrichment analysis for RNA-seq datasets? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

6.  Use the topGO package (https://bioconductor.org/packages/release/bioc/html/topGO.html) to re-do the GO term analysis. Compare and contrast the results with what has been obtained using the `gProfileR` package. Which tool is faster, `gProfileR` or topGO? Why? [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

7. Given a gene set annotated for human, how can it be utilized to work on _C. elegans_ data? (Hint: See `biomaRt::getLDS`). [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

8. Import curated pathway gene sets with Entrez identifiers from the [MSIGDB database](http://software.broadinstitute.org/gsea/msigdb/collections.jsp) and re-do the GSEA for all curated gene sets. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


### Removing unwanted variation from the expression data

For the exercises below, use the datasets at: 
```
counts_file <- system.file('extdata/rna-seq/SRP049988.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP049988.colData.tsv', 
                           package = 'compGenomRData')
```

1. Run RUVSeq using multiple values of `k` from 1 to 10 and compare and contrast the PCA plots obtained from the normalized counts of each RUVSeq run. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

2. Re-run RUVSeq using the `RUVr()` function. Compare PCA plots from `RUVs`, `RUVg` and `RUVr` using the same `k` values and find out which one performs the best. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Do the necessary diagnostic plots using the differential expression results from the EHF count table. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Use the `sva` package to discover sources of unwanted variation and re-do the differential expression analysis using variables from the output of `sva` and compare the results with `DESeq2` results using `RUVSeq` corrected normalization counts. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


