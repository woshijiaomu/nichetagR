## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installation, eval=FALSE-------------------------------------------------
# #Install dependencies
# install.packages(c("Seurat","Polychrome","stringdist","stringr","ggplot2","cowplot","igraph","openxlsx","reshape2","viridis","fields"))
# if(!require(devtools)){
# 	install.packages("devtools")
# }
# #install nichetagR
# if(!require(nichetagR)){
#   devtools::install_github("woshijiaomu/nichetagR")
# }

## ----load_package, eval = TRUE------------------------------------------------
library(nichetagR)
data(example)

## ----inspect_input_scObject---------------------------------------------------
#calculate connectome
nnt=Dnichenetwork(scObject,groupby="cell_clusters")
summary(nnt)

## ----visualization------------------------------------------------------------
#draw set-set network
print_nichenetwork(nnt,file="nichenetwork.pdf",vertex.label.cex=0.1,vsize=3,esize=1,margin=c(0,0.3,0,0))
#draw quanlity control figures
print_nichetag(nnt, file = "nichetag.pdf")
print_clustertag(nnt, file = "cluster.pdf")
tag_cancer_noncancer(nnt, file = "celltype.pdf")
tag_cci(nnt, file = "cci.pdf")
settype(nnt, file="settype.pdf")
tag_cellsettype(nnt, file="tag2celltype_settype.pdf")
#draw niche-niche network
print_nichenet(nnt,file="niche2nichenetwork.pdf",vsize = 10)

## ----export_data--------------------------------------------------------------
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
write.csv(set2cell.df,file = "set2cell.csv")

#FMUs and the sets it contains
niche2set=nnt$niche
niche2set.vec=unlist(lapply(names(niche2set),function(x){
  res=paste(x,names(niche2set[[x]]),sep=":")
  res
}))
niche2set.df=as.data.frame(stringr::str_split_fixed(niche2set.vec,":",n=2))
colnames(niche2set.df) = c("niche","set_code")
niche2set.df$set_ID=nnt$setID[niche2set.df$set_code]
write.csv(niche2set.df,file = "niche2set.csv")

#set expression matrix
setmatrix=set_expr(scObject,nnt,seurat_layer="counts")
write.csv(setmatrix,"set_matrix.csv")

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

