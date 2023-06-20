# Import libraries ----
library(AnnotationHub)
library(rhdf5)
library(stringr)
library(tidyverse)
library(tximport)

# Get .h5 paths ----

files.cwd <- list.files(pattern = "\\.h5$")
paths <- file.path(files.cwd,"abundance.h5")

# Verify the existence of the paths
if (!all(file.exists(paths))) stop("All given paths aren't working!")

# Get column names for DESeq2
col.names <- sapply(X = files.cwd,
                    FUN = str_extract,
                    pattern = "(TGFb|Yoda)+_(WT_|[0-9]uM_)?(1|2)") %>% 
  as.character() %>% 
  '['(1:length(files.cwd))

# Annotation ----

# Retrieve db for annotation
hub <- AnnotationHub()
req <- query(hub, "EnsDb.Hsapiens.v103")
ensdb <- hub[[names(req)]]

# get annotations 
suppressWarnings({
  tx <- transcripts(ensdb, columns = c("tx_id", "gene_name")) %>% 
    as_tibble() %>% 
    dplyr::rename(target_id = tx_id) %>% 
    dplyr::select(target_id, gene_name)
})

# Import kallisto .h5 ----
decisionSep <- grepl("TGFb", paths)

myPaths <- list(
  paths.tgfb <- paths[decisionSep],
  paths.yoda <- paths[!decisionSep]
)

column.names <- list(
  col.names.tgfb <- col.names[decisionSep],
  col.names.yoda <- col.names[!decisionSep]
) %>% 
  setNames(c("TGFb","YODA"))

txi.kallistos <- lapply(myPaths, function(x){
  new.kallisto <- tximport(x, 
                           type = "kallisto", 
                           tx2gene = tx, 
                           ignoreTxVersion = T)
})

for (i in seq_along(txi.kallistos)){
  colnames(txi.kallistos[[i]]$abundance) <- column.names[[i]]
  colnames(txi.kallistos[[i]]$counts) <- column.names[[i]]
  colnames(txi.kallistos[[i]]$length) <- column.names[[i]]
}
