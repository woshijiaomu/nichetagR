# NicheTagR: Reconstruction of Spatial and Lineage Relationships from CellTag-labeled Single-cell Data

**Date:** May 2025  

---

## Introduction

Thank you for your interest in the **nichetagR** approach. This guide describes how to apply the package to your intercellular RNA-barcode delivery experiment (*Multiplexed Oligonucleotide Spatial Ancestry Interaction Coding*, **MOSAIC**). MOSAIC utilizes conditionally expressed, proximally transmissible `sLP-mCherry-MCP-MS2-CellTag` barcodes to trace the spatiotemporal lineages and communicating networks within intact tissues at single-cell resolution.  

The goal of the **nichetagR** package is to use CellTag labels from single-cell transcriptomics data to reconstruct the lineage and spatial proximity relationships of mCherry-positive cells. To learn more about this type of experiment and details about the **nichetagR** approach, please consult 👉 *Bin et al., 2025* ([insert paper URL here]).

This vignette covers:

1. Package installation  
2. Setting up your input data  
3. Running the `nichetagR` computational pipeline  
4. Exporting data and session information  

---

## Package Installation

You can install the package and dependencies using:

```r
# Install dependencies
install.packages(c("Seurat", "Polychrome", "stringdist", "stringr",
                   "ggplot2", "cowplot", "igraph", "openxlsx",
                   "reshape2", "viridis", "fields"))

if (!require(devtools)) {
  install.packages("devtools")
}

# Install nichetagR from GitHub
if (!require(nichetagR)) {
  devtools::install_github("woshijiaomu/nichetagR")
}
```

## setting up your input data

Now that the package is installed, we load the package and attach its example data. We also load Seurat that will be used to manipulate Seurat Objects in R.

```r
library(nichetagR)
data(example)
```

The example data is a single-cell transcriptome (10XGenomics) which was performed in the mCherry+ cells sorted from early (2-week) metastatic lung tissues 4 days after Dox exposure. In addition to gene expression data, cell classification information and CellTag expression profiles within the object are utilized to identify spatial proximity relationships. In this guide we will use the current directory, but you can use the directory of your choice. That is all we need in terms of inputs!  Below is the pre-processsing workflow for obtaining this data:

1. To enable the detection of exogenous genes and tags, a custom mouse reference genome was generated for Cell Ranger (version 8.0.1). The coding sequences of mCherry, BFP, and all relevant CellTag sequences were appended to the mouse genome FASTA file, and their corresponding gene annotations were added to the GTF file under a new gene category. Using these modified files, a new Cell Ranger reference was built using the “cellranger mkref” command. Raw data were then processed using the cellranger count pipeline, aligning reads to the custom reference genome and generating filtered gene-barcode matrices for downstream analysis.
2. The matrix data was analyzed using the Seurat package (version 4.0.0), and cell clusters were annotated using the SingleR package (version 2.10.0).

## Running the nichetagR computational pipeline

Now that we have our gene expresion data, celltag expression profiles and cell type annotation, we have everything we need to run the pipeline.

All the pipeline is mainly packaged into a single function, `Dnichenetwork()`.

Here is a description of the basic arguments it requires:

1. scObject: Seurat object, input seurat object in which tag expression matrix is merged in gene expression matrix
2. groupby: character, a column name in meta data of the Seurat object, used for set definition
3. share_method: min,max,mean, method to compute connection strength between different sets
4. highconfidence: TorF, if or not to distinguish high confident tag-sender
5. send_cutoff1: int, count cutoff for primary celltag asignment
6. send_cutoff2: numeric, percentage cutoff for primary celltag asignment
7. direction: True or False, whether or not to distinguish tag-sender or tag-receiver

We have everything we need in our input object to fill these arguments:
```r
#compute connectome
nnt=Dnichenetwork(scObject,groupby="cell_clusters")
summary(nnt)
```

This function returns a list containing connectome information for all niches or FMUs, which mainly includes:

1. tag_matrix: data.frame, CellTag expression matrix of all cells
2. cell_type: factor, cell type or state used for set definition
3. cell_tags: numeric, total expression counts of each CellTag
4. niche_size: integer, covered cell number of each CellTag
5. niche_celltypes: integer, covered cell types or cell states by each CellTag
6. niche_interaction: data.frame, overlapped cell numbers of FMU-FMU pairs
7. cell_barcode: character, code of each cell, celltagcode+celltypecode
8. set list: each element represent a CellTag expression matrix of each set
9. set_size: integer, covered cell number of each set
10. niche list: mean CellTag expression of covered cells of all covered set for each CellTag
11. set_interaction: data.frame, overlapped tag expression of set-set pairs
12. setID: integer, order of sets when make graph with graph_from_data_frame in igraph

Visualization for connectome data encompass the following functions:

1. print_nichenetwork: set network - each vertex represents a set, and edges between vertices indicate shared CellTags between two sets
2. print_nichetag: Bar plot - displays the cell number and cell types covered by each FMU
3. print_clustertag: Density plot - shows the distribution of CellTag expression levels within each cell cluster
4. tag_cancer_noncancer: Bar plot - comparatively displays the expression levels of each CellTag in cancer cells and stromal cells
5. tag_cci: Bar plot - shows the distribution of sender and receiver cell counts corresponding to different CellTags
6. settype: Bar plot - displays the cell types and quantities of sets within each FMU
7. tag_cellsettype: Bar plot - shows the types and quantities of sets or cells covered by each CellTag
8. print_nichenet: set network - each vertex represents an FMU, and edges between vertices indicate the presence of shared CellTags between two FMUs

We run these functions for visualization:
```r
#graph set-set network
print_nichenetwork(nnt,file="Set_Network_Shared_CellTags.pdf",vertex.label.cex=0.1,vsize=3,esize=1,margin=c(0,0.3,0,0))
#graph quanlity control figures
print_nichetag(nnt, file = "Niche_CellNumber_CellType.pdf")
print_clustertag(nnt, file = "Cluster_CellTag_Expression_Density.pdf")
tag_cancer_noncancer(nnt, file = "CellTag_Expression_Cancer_vs_Noncancer.pdf")
tag_cci(nnt, file = "CellTag_Sender_Receiver_Distribution.pdf")
settype(nnt, file="Niche_Set_CellType_Composition.pdf")
tag_cellsettype(nnt, file="CellTag_Coverage_Set_Cell.pdf")
#graph niche-niche network
print_nichenet(nnt,file="Niche_Network_Shared_CellTags.pdf",vsize = 10)
```

## Export Data

In addition to the visualization parts above, we can also extract from the connectome object: 
1. the list of sets contained in each FMU
2. the list of cells contained in each set
3. the average expression matrix for each set

The information can be written to local CSV files, which can be useful to further analyze the data in R:
```r
#sets and the cells it contains
set2cell=nnt$set
set2cell.vec=unlist(lapply(names(set2cell),function(x){
  y=set2cell[[x]]
  res=rep(x,nrow(y))
  names(res)=rownames(y)
  res
}))
set2cell.df=data.frame(set_code=set2cell.vec,cell=names(set2cell.vec))
set2cell.df$set_ID=nnt$setID[set2cell.df$set_code]
write.csv(set2cell.df,file = "Set_Cell_Mapping.csv")

#FMUs and the sets it contains
niche2set=nnt$niche
niche2set.vec=unlist(lapply(names(niche2set),function(x){
  res=paste(x,names(niche2set[[x]]),sep=":")
  res
}))
niche2set.df=as.data.frame(stringr::str_split_fixed(niche2set.vec,":",n=2))
colnames(niche2set.df) = c("niche","set_code")
niche2set.df$set_ID=nnt$setID[niche2set.df$set_code]
write.csv(niche2set.df,file = "Niche_Set_Membership.csv")

#set expression matrix
setmatrix=set_expr(scObject,nnt,seurat_layer="counts")
write.csv(setmatrix,"Set_Expression_Matrix.csv")
```

## Conclusion

Thank you for following this guide, I hope you made it through the end without too much headache and that it was informative. If you encounter bugs, feel free to raise an issue on nichetagR's [github](https://github.com/woshijiaomu/nichetagR/issues).

