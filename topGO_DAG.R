##########################################
#     prj: Thesis
#     Assignment: GO DAG
#     Author: Shwan Wang
#     Date: Mar. 19, 2021
##########################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require('topGO')) BiocManager::install('topGO',update = FALSE)
if (!require('getopt')) install.packages('getopt')
options(stringsAsFactors = F)
library(getopt)
library(topGO)
# args --------------------------------------------------------------------

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'background', 'b', 1, 'character', 'inputfile: background file, 1st column is GeneID, 2nd is collapsed GOID, Separated by commas',
  'querylist', 'q', 1, 'character', 'inputfile: Interested gene list, eg: DEGs.',
  'ont', 'o', '1', 'integer', 'Ontology: "BP","MF","CC".',
  'firstSigNodes', 'f', '1', 'integer', 'firstSigNodes: number of significant nodes (retangle nodes in the graph).'
),byrow = T, ncol = 5)
args = getopt(command)

## help information
if (!is.null(args$help)) {
  cat(paste(getopt(command, usage = T), "\n"))
  #  q(status=1)
}

if (is.null(args$ont)){
  args$ont = "BP"
}

if (is.null(args$firstSigNodes)){
  args$firstSigNodes <- 10
}
bg = args$background
qlist = args$querylist
ont = args$ont
fSN = args$firstSigNodes

# Data cleaning -----------------------------------------------------------

gene2go = readMappings(bg)
genenames <- names(gene2go)
gene_data <- read.table(qlist)
gene_id <- gene_data[,1]
genelist <- factor(as.integer(genenames %in% gene_id))
names(genelist) <- genenames

GOdata <- new("topGOdata", allGenes = genelist, ontology = ont, 
              annot = annFUN.gene2GO, gene2GO = gene2go)

printGraph(GOdata, resultFisher, firstSigNodes = fSN, fn.prefix = "sampleFile", useInfo = "all", pdfSW = TRUE)