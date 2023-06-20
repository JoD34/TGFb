library(Biobase)
library(clusterProfiler)
library(DT)
library(enrichplot)
library(limma)
library(gplots)
library(gprofiler2)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(tidyverse)
library(xlsx)

# GO analysis ----
contrast_list <- contrast_list[2:4]

pdf_width <- 11
pdf_height <- 8.5
pdf("go_graph.pdf", width=pdf_width, height=pdf_height)

lapply(contrast_list, function(list){
  lapply(list, function(sublist){
    gost.res <- gost(rownames(sublist), 
                     organism = "hsapiens", 
                     correction_method = "fdr")
    
    gostyplot <- gostplot(gost.res, interactive = F, capped = F)
    
    publish_gostplot(
      gostyplot,
      filename= NULL,
      width = pdf_width,
      height = pdf_height)
    
  })
})

# Extract data for GO enrichment

for( i in seq_along(GO_list)){
  name = names(GO_list)[4]
  write.xlsx(x = GO_list[[4]], 
             file = "TGFb_geneList.xlsx", 
             sheetName = paste0(name, "_fullData"),
             append = T)
}
for (i in seq_along(txi.kallistos)){
  name = names(txi.kallistos)[i]
  write.xlsx(x = txi.kallistos[[2]]$counts,
             file = "background_TGFb.xlsx",
             sheetName = names(txi.kallistos)[2],
             append=T) 
}

# GSEA analysis ----

# Get gsea for Homo sapiens
hs_gsea <- msigdbr(species = "Homo sapiens")

gsea_list <- hs_gsea %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat)

# GSEA data input
GSEA_data_list <- lapply(GSEA_data_list, function(x) {
  data <- x$log2FoldChange %>% 
    'names<-'(rownames(x)) %>% 
    sort(decreasing = T)
})

# Make list categories
myGSEA.res <- list()
cat_list <- c("C2","C5")

# Make GSEA data
for (i in 1:length(cat_list)) {
  mySig <- msigdbr(species = "Homo sapiens", category = cat_list[i]) %>%
    dplyr::select(gs_name, gene_symbol) 
  
  myGSEA.res <- lapply(GSEA_data_list, function(gsea){
    gsea <- GSEA(gsea, TERM2GENE = mySig,verbose = F)
  })
}
  
myGSEA.df <- lapply(myGSEA.res, function(gsea){as_tibble(gsea@result)})

# Send GSEA to excel table
for(i in 1:length(myGSEA.df)){
  write.xlsx(myGSEA.df[i], 
             file= "GSEAdata.xlsx",
             sheetName = names(myGSEA.df)[i],
             col.names = T,
             row.names = T,
             append = T)
  }

# Visualize GSEA output - interactive table 
datatable(myGSEA.df[[1]],
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in TNF\u03B1 inflammation',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2) 
 
# Plot GSEA output + graph 
for(i in 1:length(myGSEA.res)){
  PDFname <- paste0("GSEAplot_", names(myGSEA.res)[i],".pdf")
  pdf(PDFname)
  elem <- myGSEA.res[[i]]
  
  for (j in 1:nrow(myGSEA.df[[i]])){
    myplotg <- gseaplot2(elem, 
                         geneSetID = j, 
                         pvalue_table = F, 
                         title = elem$Description[j])
    print(myplotg)
    }
  dev.off()
}

# DataFrame for "bubble plot"
for(i in 1:length(myGSEA.df)){
  name <- str_split_1(str_extract(names(myGSEA.df)[1],"[0-9]+h_[0-9]+h"),"_") %>% 
    paste0("TNF\u03B1-", .)
  
  myGSEA.df[[i]] <- myGSEA.df[[i]] %>% 
    mutate(phenotype = case_when(NES > 0 ~ name,
                                 NES < 0 ~ "Wild Type"))
  }

# Bubble plot summarizing 'y' signatures across 'x' phenotypes - C2
pdf("GSEA_bubble.pdf", width = 11, height=8.5)

for (i in 1:length(myGSEA.df)){
  
  bubbly <- ggplot(myGSEA.df[[i]][1:20,]) +
    aes(x = phenotype, y = ID) +
    geom_point(aes(size = setSize, 
                   color = NES, 
                   alpha = -log10(p.adjust))) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw()
  
  print(bubbly)
}
dev.off()
