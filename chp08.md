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

DEresults.Sig <- subset(DEresults, padj < 0.1)
summary(DEresults.Sig)
```

Now, you are ready to do the following exercises: 

1. Make a volcano plot using the differential expression analysis results. (Hint: x-axis denotes the log2FoldChange and the y-axis represents the -log10(pvalue)). [Difficulty: **Beginner**]

**solution:**
```{r}
library(ggplot2)

DEGs <- data.frame(DEresults.Sig) %>% mutate()

DEGs <- DEGs[DEGs$pvalue<0.05,]
DEGs$regulation <- (DEGs$log2FoldChange) > 0
sum(is.na(DEGs))
#replace the values
DEGs$regulation <- as.character(DEGs$regulation)
DEGs$regulation[DEGs$regulation == "TRUE"] <- "Up-regulate"
DEGs$regulation[DEGs$regulation == "FALSE"] <- "Down-regulate"

volc = ggplot(data = DEGs, aes(x=log2FoldChange, y=-log10(pvalue))) +
        geom_point(aes(colour=regulation), size =1, alpha =0.3) +
        theme_classic()
volc
```

2. Use DESeq2::plotDispEsts to make a dispersion plot and find out the meaning of this plot. (Hint: Type ?DESeq2::plotDispEsts) [Difficulty: **Beginner**]

**solution:** 
```{r}
dds.disp <- estimateSizeFactors(dds)
dds.disp <- estimateDispersions(dds.disp)
plotDispEsts(dds.disp)
```

3. Explore `lfcThreshold` argument of the `DESeq2::results` function. What is its default value? What does it mean to change the default value to, for instance, `1`? [Difficulty: **Intermediate**]

**solution:** The default value is `0`. With default value total 8170 DEGs were identified and for `1` 1681 DEGs were found
```{r}
DEresults.1 <- results(dds, lfcThreshold=1)
DEresults.Sig1 <- subset(DEresults.1, padj < 0.1)

#summary of results with default value lfcThreshold=0
#summary of results with lfcThreshold=1
summary(DEresults.Sig)
summary(DEresults.Sig1)
```

4. What is independent filtering? What happens if we don't use it? Google `independent filtering statquest` and watch the online video about independent filtering. [Difficulty: **Intermediate**]

**solution:** The goal of independent filtering is to filter out those tests from the procedure that have no, or little chance of showing significant evidence, without even looking at their test statistic. Typically, this results in increased detection power at the same experiment-wide type I error. Here, we measure experiment-wide type I error in terms of the false discovery rate.
```{r,echo=FALSE,eval=FALSE}
use <- DEresults$baseMean > metadata(DEresults)$filterThreshold
h1 <- hist(DEresults$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(DEresults$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

```

5. Re-do the differential expression analysis using the `edgeR` package. Find out how much DESeq2 and edgeR agree on the list of differentially expressed genes. [Difficulty: **Advanced**] 

**solution:** `70%` common DEGs identified using both packages, `27%` extra DEGs are in DESeq2 analysis and `3%` extra genes in EdgeR analysis.
```{r}
library(edgeR)
library.sizes <- colSums(counts)
dgeFull <- DGEList(counts = counts, 
                   samples = colData, 
                   lib.size = library.sizes, 
                   group = colData$group)

dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) > 1, ],
                  samples = colData, 
                   lib.size = library.sizes, 
                   group = colData$group)

#estimate the normalization factors
dgeFull <- calcNormFactors(object = dgeFull, method = "TMM")

#estimate common and tagwise dispersion
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull

#MDS plot
plotMDS.DGEList(dgeFull,main="MDS plot",col=rep(1:2,each=5))

#
dgeTest <- exactTest(dgeFull)

resFilt <- topTags(dgeTest, n=nrow(dgeTest$table))

resFilt.table <- resFilt$table
E.DEGs<- resFilt.table[resFilt.table$FDR<0.1,]
E.DEGs <- E.DEGs[E.DEGs$PValue<0.05,]


E.DEGs$regulation <- (E.DEGs$logFC) > 0
#replace the values
E.DEGs$regulation <- as.character(E.DEGs$regulation)
E.DEGs$regulation[E.DEGs$regulation == "TRUE"] <- "Up-regulate"
E.DEGs$regulation[E.DEGs$regulation == "FALSE"] <- "Down-regulate"

E.volc = ggplot(data = E.DEGs, aes(x=logFC, y=-log10(PValue))) +
        geom_point(aes(colour=regulation), size =1, alpha =0.3) +
        theme_classic()
E.volc

library(VennDiagram)
Dseq.list <- rownames(DEGs)
E.DEGs.list <- rownames(E.DEGs)
gene_list = list(Dseq.list, E.DEGs.list)

library(ggvenn)
ggvenn(gene_list)
p1 <- ggVennDiagram(gene_list, label_alpha = 0,
                    category.names = c("DEseq genes","EdgeR genes"),)
p1 + scale_fill_distiller(palette = "RdBu")

```

6. Use the `compcodeR` package to run the differential expression analysis using at least three different tools and compare and contrast the results following the `compcodeR` vignette. [Difficulty: **Advanced**]

**solution:**
library(compcodeR)
sample.annot <- colData
colnames(sample.annot)[2] <- "condition"

info.parameters <- list(dataset = "CancerData", uID = "123456")
cpd <- compData(count.matrix = counts, 
                sample.annotations = sample.annot, 
                info.parameters = info.parameters)
check_compData(cpd)

saveRDS(cpd, file = "cpd.rds")

runDiffExp(data.file = "cpd.rds", 
           result.extent = "voom.limma", Rmdfunction = "voom.limma.createRmd", 
           output.directory = ".", norm.method = "TMM")

runDiffExp(data.file = "cpd.rds", 
           result.extent = "edgeR.exact", Rmdfunction = "edgeR.exact.createRmd", 
           output.directory = ".", norm.method = "TMM", 
           trend.method = "movingave", disp.type = "tagwise")
runDiffExp(data.file = "cpd.rds", result.extent = "ttest", 
           Rmdfunction = "ttest.createRmd", 
           output.directory = ".", norm.method = "TMM")

runDiffExp(data.file = "cpd.rds", result.extent = "DESeq2", 
           Rmdfunction = "DESeq2.createRmd",
           output.directory = ".", fit.type = "parametric", test = "Wald" )


runComparisonGUI(input.directories = ".", 
                 output.directory = ".", recursive = FALSE)
```

### Functional enrichment analysis

1. Re-run gProfileR, this time using pathway annotations such as KEGG, REACTOME, and protein complex databases such as CORUM, in addition to the GO terms. Sort the resulting tables by columns `precision` and/or `recall`. How do the top GO terms change when sorted for `precision`, `recall`, or `p.value`? [Difficulty: **Beginner**]

**solution:** If we sort the table with `recall` then the highest recall top terms is at the top. And same results, in the case of `p.value`, and `precision`.
```{r}
library(DESeq2)
library(gProfileR)
library(knitr)
genesOfInterest <- rownames(DEGs)

#calculate enriched GO terms
goResults <- gprofiler(query = genesOfInterest,  
                       organism = 'hsapiens', 
                       src_filter = c('GO','KEGG', 'REAC', 'HPA', 'CORUM'), 
                       hier_filtering = 'moderate')

#Let's define the first gene set as the list of genes from one of the
#significant GO terms found in the GO analysis. order go results by precision
goResults.prec <- goResults[order(goResults$precision),]
goResults.recall <- goResults[order(goResults$recall),]
goResults.pval <- goResults[order(goResults$p.value),]

```

2. Repeat the gene set enrichment analysis by trying different options for the `compare` argument of the `GAGE:gage`
function. How do the results differ? [Difficulty: **Beginner**]

**solution:** We can see that the random gene set shows significant up-regulation and no significant down-regulation and same in the case of  top GO terms.
```{r}

goResults <- goResults[order(goResults$p.value),]
#restrict the terms that have at most 100 genes overlapping with the query
go <- goResults[goResults$overlap.size < 100,]
# use the top term from this table to create a gene set 
geneSet1 <- unlist(strsplit(go[1,]$intersection, ','))
#Define another gene set by just randomly selecting 25 genes from the counts
#table get normalized counts from DESeq2 results
normalizedCounts <- DESeq2::counts(dds, normalized = TRUE)
geneSet2 <- sample(rownames(normalizedCounts), 25)

geneSets <- list('top_GO_term' = geneSet1,
                 'random_set' = geneSet2)

library(gage)
#use the normalized counts to carry out a GSEA. 
gseaResults <- gage(exprs = log2(normalizedCounts+1), 
                    ref = match(rownames(colData[colData$group == 'CTRL',]), 
                                colnames(normalizedCounts)), 
                    samp = match(rownames(colData[colData$group == 'CASE',]), 
                                 colnames(normalizedCounts)),
                    gsets = geneSets, compare = 'as.group')
gseaResults$greater
gseaResults$less

```

3. Make a scatter plot of GO term sizes and obtained p-values by setting the `gProfiler::gprofiler` argument `significant = FALSE`. Is there a correlation of term sizes and p-values? (Hint: Take -log10 of p-values). If so, how can this bias be mitigated? [Difficulty: **Intermediate**]

**solution:** If the term size is increased, its significantly enriched. So term size has +ve correlation with `P-values`. This biasness can be mitigated by independent filtering.
```{r}
#calculate enriched GO terms
goResults <- gprofiler(query = genesOfInterest,
                       organism = 'hsapiens', 
                       src_filter = 'GO', 
                       hier_filtering = 'moderate',
                       significant = FALSE)
goResults$GeneRatio <- goResults$term.size/ goResults$query.size

GO.plot <- ggplot(goResults, aes(x=GeneRatio, y = term.id,
                                         size = term.size, color = -log10(p.value) )) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        theme_bw()
GO.plot

```

4. Do a gene-set enrichment analysis using gene sets from top 10 GO terms. [Difficulty: **Intermediate**]

**solution:**
```{r}
#calculate enriched GO terms
goResults <- gprofiler(query = genesOfInterest,
                       organism = 'hsapiens', 
                       src_filter = 'GO', 
                       hier_filtering = 'moderate')

goResults <- goResults[order(goResults$p.value),]
#restrict the terms that have at most 100 genes overlapping with the query
go <- goResults[goResults$overlap.size < 100,]
# use the top term from this table to create a gene set 
for (i in 1:10) {
        geneSets[[i]] <-unlist(strsplit(go[i,]$intersection, ',')) %>% as.list()
}

library(gage)
#use the normalized counts to carry out a GSEA. 
gseaResults <- gage(exprs = log2(normalizedCounts+1), 
                    ref = match(rownames(colData[colData$group == 'CTRL',]), 
                                colnames(normalizedCounts)), 
                    samp = match(rownames(colData[colData$group == 'CASE',]), 
                                 colnames(normalizedCounts)),
                    gsets = as.list(geneSets))
gseaResults$greater
gseaResults$less
```

5. What are the other available R packages that can carry out gene set enrichment analysis for RNA-seq datasets? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

6.  Use the topGO package (https://bioconductor.org/packages/release/bioc/html/topGO.html) to re-do the GO term analysis. Compare and contrast the results with what has been obtained using the `gProfileR` package. Which tool is faster, `gProfileR` or topGO? Why? [Difficulty: **Advanced**]

**solution:** gProfileR is faster than topGO. And results are different, because gProfileR using a deprecated package relying on outdated data.
```{r}
library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
require(org.Hs.eg.db)


# gene list vector with pvalue
geneList <- DEGs$pvalue
names(geneList) <- rownames(DEGs)

# For Biological process
selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="BP",
                          feasibleGenes=NULL, 
                          mapping="org.Hs.eg.db", 
                          ID="symbol")
GOdata.BP <- new("topGOdata",
              ontology="BP",
              allGenes=geneList,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=geneSelectionFun,
              nodeSize=10)

#Kolmogorov-Smirnov Testing
results.ks.BP <- runTest(GOdata.BP, algorithm = "classic", statistic = "ks")

# For Molecular Function
selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="MF",
                          feasibleGenes=NULL, 
                          mapping="org.Hs.eg.db", 
                          ID="symbol")
GOdata.MF <- new("topGOdata",
              ontology="MF",
              allGenes=geneList,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=geneSelectionFun,
              nodeSize=10)


results.ks.MF <- runTest(GOdata.MF, algorithm = "classic", statistic = "ks")

# For Cellular Components
selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="CC",
                          feasibleGenes=NULL, 
                          mapping="org.Hs.eg.db", 
                          ID="symbol")
GOdata.CC <- new("topGOdata",
              ontology="CC",
              allGenes=geneList,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=geneSelectionFun,
              nodeSize=10)


results.ks.CC <- runTest(GOdata.CC, algorithm = "classic", statistic = "ks")


goEnrichment.BP <- GenTable(GOdata.BP,
                         KS = results.ks.BP,
                         orderBy = "KS",
                         topNodes = 100,
                         numChar = 99)

goEnrichment.MF <- GenTable(GOdata.MF,
                         KS = results.ks.MF,
                         orderBy = "KS",
                         topNodes = 100,
                         numChar = 99)

goEnrichment.CC <- GenTable(GOdata.CC,
                         KS = results.ks.CC,
                         orderBy = "KS",
                         topNodes = 100,
                         numChar = 99)

All.GO <- rbind(goEnrichment.MF, goEnrichment.BP, goEnrichment.CC)

All.GO$KS <- as.numeric(All.GO$KS)
All.GO <- All.GO[All.GO$KS < 0.05,] # filter terms for KS p<0.05
All.GO

par(cex = 0.3)
showSigOfNodes(GOdata.CC, score(results.ks.CC), firstSigNodes = 3, useInfo = "def")
showSigOfNodes(GOdata.BP, score(results.ks.BP), firstSigNodes = 3, useInfo = "def")
showSigOfNodes(GOdata.MF, score(results.ks.MF), firstSigNodes = 3, useInfo = "def")
# using Fisher exect test
#resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)
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

**solution:** when `K=9`  `RUVs` is better. And for `RUVg` `K=1` is better.
```{r}
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
colData <- read.table(colData_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
# remove 'width' column from counts
countData <- as.matrix(subset(counts, select = c(-width)))
# # simplify condition descriptions
colData$source_name <- ifelse(colData$group == 'CASE', 
                              'EHF_overexpression', 'Empty_Vector')

library(EDASeq)

# create a seqExpressionSet object using EDASeq package 
set <- newSeqExpressionSet(counts = countData,
                           phenoData = colData)

library(RUVSeq)

# RUVs ----
# make a table of sample groups from colData 
differences <- makeGroups(colData$group)

for(k in 1:9) {
  set_s <- RUVs(set, unique(rownames(set)), 
                k=k, differences) #all genes
  plotPCA(set_s, col=as.numeric(colData$group), 
          cex = 0.9, adj = 0.5, 
          main = paste0('with RUVs, k = ',k), 
          ylim = c(-1, 1), xlim = c(-0.6, 0.6))
}

# RUVg ----
HK_genes <- read.table(file = system.file("extdata/rna-seq/HK_genes.txt",
                                          package = 'compGenomRData'),header = FALSE)
# let's take an intersection of the house-keeping genes with the genes available
# in the count table
house_keeping_genes <- intersect(rownames(set), HK_genes$V1)

# now, we use these genes as the empirical set of genes as input to RUVg.
# we try different values of k and see how the PCA plots look 

for(k in 1:9) {
  set_g <- RUVg(x = set, cIdx = house_keeping_genes, k = k)
  plotPCA(set_g, col=as.numeric(colData$group), cex = 0.9, adj = 0.5, 
          main = paste0('with RUVg, k = ',k), 
          ylim = c(-1, 1), xlim = c(-1, 1), )
}

```

2. Re-run RUVSeq using the `RUVr()` function. Compare PCA plots from `RUVs`, `RUVg` and `RUVr` using the same `k` values and find out which one performs the best. [Difficulty: **Intermediate**]

**solution:** `RUVs` is perform better as compare to `RUVr` and `RUVg`.
```{r}
library(edgeR)

# Residuals from negative binomial GLM regression of UQ-normalized
# counts on covariates of interest, with edgeR
design <- model.matrix(~ colData$group)
dgeFull <- DGEList(counts = countData, 
                   samples = colData, 
                   group = colData$group)

dgeFull <- calcNormFactors(dgeFull, method="upperquartile")

dgeFull <- estimateGLMCommonDisp(dgeFull, design)
dgeFull <- estimateGLMTagwiseDisp(dgeFull, design)

fit <- glmFit(dgeFull, design)
res <- residuals(fit, type="deviance")
# RUVr normalization (after UQ)
setUQ <- betweenLaneNormalization(set, which="upper")
controls <- rownames(set)
set_r <- RUVr(set, controls, k=1, res)

par(mfrow = c(3, 3))
for(k in 1:9) {
  set_r <- RUVr(set, controls, k=k, res)
  plotPCA(set_r, col=as.numeric(colData$group), cex = 0.9, adj = 0.5, 
          main = paste0('with RUVr, k = ',k), 
          ylim = c(-1, 1), xlim = c(-1, 1), )
}

```

3. Do the necessary diagnostic plots using the differential expression results from the EHF count table. [Difficulty: **Intermediate**]

**solution:**
```{r}
library(DESeq2)
#set up DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData, 
                              design = ~ group)
# filter for low count genes
dds <- dds[rowSums(DESeq2::counts(dds)) > 10]

# insert the covariates W1 and W2 computed using RUVs into DESeqDataSet object
colData(dds) <- cbind(colData(dds), 
                      pData(set_s)[rownames(colData(dds)), 
                                   grep('W_[0-9]', 
                                        colnames(pData(set_s)))])
# update the design formula for the DESeq analysis (save the variable of
# interest to the last!)
design(dds) <- ~ W_1 + W_2 + group 
# repeat the analysis 
dds <- DESeq(dds)
# extract deseq results 
res <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
res <- subset(res, padj < 0.1)
summary(res)

# diagnostic plots
DESeq2::plotMA(object = res, ylim = c(-5, 5))

library(ggplot2)
ggplot(data = as.data.frame(res), aes(x = pvalue)) + 
  geom_histogram(bins = 100)

countsNormalized <- DESeq2::counts(dds, normalized = TRUE)
selectedGenes <- names(sort(apply(countsNormalized, 1, var), 
                            decreasing = TRUE)[1:500])

plotPCA(countsNormalized[selectedGenes,], 
        col = as.numeric(colData$group), adj = 0.5, 
        xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.6))

library(EDASeq)
par(mfrow = c(1, 2))
plotRLE(countData, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), 
        main = 'Raw Counts')
plotRLE(DESeq2::counts(dds, normalized = TRUE), 
        outline=FALSE, ylim=c(-4, 4), 
        col = as.numeric(colData$group), 
        main = 'Normalized Counts')
 
```

4. Use the `sva` package to discover sources of unwanted variation and re-do the differential expression analysis using variables from the output of `sva` and compare the results with `DESeq2` results using `RUVSeq` corrected normalization counts. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


