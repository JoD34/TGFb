# Import libraries ----
library(DESeq2)
library(EnhancedVolcano)
library(flow)
library(ggplot2)
library(ggrepel)
library(plotly)
library(pheatmap)
library(RColorBrewer)
library(ReportingTools)
library(scales)
library(stringr)
library(tidyverse)
library(vsn)

# functions ----
getNames <- function(list){
  name <- names(list) %>% 
    str_extract_all("(yoda\\.(1|2)u(m|M)|tgfb|wt)")%>% 
    unlist() %>% 
    sapply(function(elem){
      elem <- if(elem %in% "tgfb") c("TGF\u03B2") else 
        if((elem %in% "wt")) "Wild Type" else
          paste0("YODA","-",sub("^.*\\.(.*)$", "\\1",elem))
    }) %>% 
    rev() 
}

# Pdf data collection
pdf("graph_tnfb_yoda.pdf", width = 11, height = 8.5, onefile = T)

# Make the design table ----
designs <- lapply(column.names, function(x){
  myConditions <- str_extract(x,"(Yoda|TGFb)(_(WT|[0-9]+uM))?") %>% 
    gsub("_",".",.) %>% 
    ifelse(. == "TGFb",paste0("TGFb.treated"), .) %>% 
    tolower(.)
  
  dataSet <- str_extract(x[1],"Yoda|TGF")
  myLevel <- if(dataSet == "Yoda") c("Yoda.wt", "Yoda.1uM", "Yoda.2uM") 
  else c("TGFb.wt", "TGFb.treated")

  myDesign <- data.frame(
    condition = myConditions,
    sample = str_extract(x, "(1|2)$")) %>% 
    mutate(
      condition = factor(condition, levels = tolower(myLevel)),
      sample = factor(sample)) %>% 
    'rownames<-'(x)
})

# DESeq2 data generations ----
dds <- lapply(seq_along(txi.kallistos), function(x){
  x <- DESeqDataSetFromTximport(txi.kallistos[[x]], 
                                  colData = designs[[x]], 
                                  design = ~ condition + sample) %>% 
    DESeq(.)
}) %>% 
  setNames(names(column.names))

# Implement contrast between conditions + correct for NAs ----
graph_list <- list(
  con_tgfb_wt = results(dds$TGFb, 
                      contrast = c("condition", "tgfb.treated","tgfb.wt")) %>% 
    na.omit(),
  con_yoda.2uM_yoda.1uM = results(dds$YODA,
                         contrast = c("condition", "yoda.2um", "yoda.1um")) %>%
    na.omit(),
  con_yoda.1uM_wt = results(dds$YODA, 
                          contrast = c("condition", "yoda.1um", "yoda.wt")) %>% 
    na.omit(),
  con_yoda.2uM_wt = results(dds$YODA, 
                          contrast = c("condition", "yoda.2um", "yoda.wt")) %>% 
    na.omit()
)

# Organize data for functional analysis ----
contrast_list <- list()
GSEA_data_list <- list()
GO_list <- list()

for(i in seq_along(graph_list)){
  
  # Extract list name
  listName <- names(graph_list[i]) %>% 
    str_extract("(yoda|tgfb)((.|_)((1|2)uM|wt))?(_yoda)?((.|_)((1|2)uM|wt))?")
  subList <- graph_list[[i]]
  
  # Filtering for p-adj < 0.01 & log2(FC) > 1.5
  filtered <- subList[subList$padj < 0.01 & 
                        abs(subList$log2FoldChange) > log2(1.5), ]
  # Ordered by p-adj
  ordered <- filtered[order(filtered$padj), ]
  
  # Get Go data
  GO_list <- append(GO_list, list(subList))
  names(GO_list)[i] <- listName
  
  # Get data for GSEA analysis
  GSEA_data_list <- append(GSEA_data_list, 
                           list(ordered["log2FoldChange"]))
  names(GSEA_data_list)[i] <- listName
  
  # Seperate the upregulated from the downregulated
  decisionSet <- ordered$log2FoldChange > 0 
  myList <- list(
    pos = ordered[decisionSet > 0, ],
    neg = ordered[!decisionSet > 0, ]
  )

  # Get data for GO analysis
  contrast_list <- append(contrast_list, list(myList))
  names(contrast_list)[i] <- listName
}

# Volcano plots ----
for(i in seq_along(graph_list)){
  
  names <- getNames(graph_list[i])

  value <- graph_list[[i]]
  labels <- data.frame(value[order(log10(value$padj)), ][1:25,])
  
  myPlot <- ggplot(as_tibble(value)) +
    aes(x = log2FoldChange, 
        y = -log10(padj)) +
    geom_point(aes(color = log2FoldChange > 0),
               show.legend = F, 
               shape = 1, 
               size = 2) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               alpha = 0.5) +
    labs(title = paste0("Volcano plot: ", names[[2]], 
                        " in relation to ", names[[1]]),
         x = expression(log[2]("Fold Change")),
         y = expression(-log[10]("adjusted pValue")),
         caption = paste("Produced on", Sys.time())) +
    theme_bw()
  
  y <- diff(layer_scales(myPlot)$y$range$range)/2
  x <- min(abs(layer_scales(myPlot)$x$range$range)) *(1/22)

  
  myPlot <- myPlot +
    annotate("text",
             label = "Downregulated", 
             alpha = 0.5, 
             angle = 90, 
             x = -x, 
             y = y,
             size = 3) +
    annotate("text",
             label = "Upregulated", 
             alpha = 0.5, 
             angle = -90, 
             x = x, 
             y = y,
             size = 3)
    print(myPlot)
}
# Apply log fold change shrinkage -> lower noises from log2(FC)----
ma_data <- list(
  con_tgfb.wt.noNA.filt = lfcShrink(dds$TGFb,
                                coef=c("condition_tgfb.treated_vs_tgfb.wt"), 
                                type = "apeglm") %>% 
    na.omit(),
  con_yoda.1uM_wt.noNA.filt = lfcShrink(dds$YODA,
                                 coef=c("condition_yoda.1um_vs_yoda.wt"), 
                                 type = "apeglm") %>% 
    na.omit(),
  con_yoda.2uM_wt.noNA.filt = lfcShrink(dds$YODA,
                                 coef=c("condition_yoda.2um_vs_yoda.wt"), 
                                 type = "apeglm") %>% 
    na.omit()
)

for (i in seq_along(ma_data)){
  names <- getNames(ma_data[i])
  value <- ma_data[[i]]

  myMA <- plotMA(value, colSig= "royalblue4", alpha = 0.01, 
         main=paste0("MAplot: ", names[[2]], " in relation to ", names[[1]]),
         ylab="log (FC)") 
  print(myMA)
}

# Effects of transformations on the variance ----
lapply(seq_along(dds), function(i){
  smallDDS <- dds[[i]]
  dataName <- names(dds)[i] %>% 
    ifelse(.=="TGFb", "TGF\u03B2",.)

# vst -> variance stabilizing transformation
vsd <- vst(smallDDS, blind=F)
meanSD <- meanSdPlot(assay(vsd))
print(meanSD)

# PCA analysis
pcaData <- plotPCA(vsd, 
                   intgroup=c("condition", "sample"),
                   returnData=T,
                   ntop=50000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
namesPCA <- attr(pcaData,"names")

# PCA plot
pca <- ggplot(pcaData) +
  aes(PC1, PC2, color=condition, shape=sample) +
  geom_point(size=3) +
  labs(title = paste("Principal Component Analysis (PCA) - ", dataName," data"),
       x = paste0(namesPCA[1], ": ", percentVar[1], " %"),
       y = paste0(namesPCA[2], ": ", percentVar[2], " %"),
       caption = paste("Produced on", Sys.time())) +
  coord_fixed() +
  theme_bw()

print(pca)


# Heatmap of the count matrix
select <- order(rowMeans(counts(smallDDS,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(smallDDS))

heat <- pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=F, annotation_col=df)
print(heat)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Generate heatmap of the distance matrix 
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample, sep=":")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(name="BuPu", n=9)))(255)
dist_mat <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Sample-to-sample Distance visualization")

print(dist_mat)
})

# Print datas in PDF
grDevices::dev.off()