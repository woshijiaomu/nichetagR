#' Calculate niche-niche connections with connectome data from nichetag method
#'
#' This fuction use tag expression matrix to infer the connection among niches labeled with different nichetags
#'
#' @param scObject    Seurat, input seurat object in which tag expression matrix is merged in gene expression matrix
#' @param groupby    character, a colname in meta data of seurat object, used for set defination
#' @param share_method    min,max,mean, method to calculate connection strength between different sets
#' @param direction    TorF, if or not to distinguish tag-sender or tag-receiver
#'
#' @return    a list contains connectome information of all niches:
#' \describe{
#'   \item{tag_matrix}{data.frame, tag expression matrix of all cells}
#'   \item{cell_type}{factor, cell type or state used for set defination}
#'   \item{cell_tags}{numeric, total expression of each tag}
#'   \item{niche_size}{integer, covered cell number of each tag}
#'   \item{niche_celltypes}{integer, covered cell types or cell states of each tag}
#'   \item{niche_interaction}{data.frame, overlapped cell numbers of tag-tag pairs}
#'   \item{code2celltype}{character,corresponding celltype of each celltypes code positons}
#'   \item{code2tag}{character, corresponding tag of each tag code positons}
#'   \item{cell_barcode}{character, code of each cell, tagcode+celltypecode}
#'   \item{set}{list, each element represent a tag expression matrix of the set}
#'   \item{set_size}{integer, covered cell number of each set}
#'   \item{niche}{list, mean tag expression of covered cells of all covered set of each tag}
#'   \item{set_interaction}{data.frame, overlapped tag expression of set-set pairs}
#'   \item{setID}{integer, order of sets when make graph with graph_from_data_frame in igraph }
#' }
#' @importFrom SeuratObject FetchData
#' @importFrom stats complete.cases sd ecdf
#' @importFrom igraph graph_from_data_frame
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork_v1(scObject,groupby="cell_clusters")
#' summary(nnt)
#' nnt=Dnichenetwork_v1(scObject,groupby="cell_clusters",share_method="mean")
#' dnnt=Dnichenetwork_v1(scObject,groupby="cell_clusters",direction = TRUE)
Dnichenetwork_v1<-function(scObject,groupby="cell_clusters",share_method="min",direction=FALSE){
  res=list()
  tags <- grep("^[ATCG]{8}$", rownames(scObject), value = TRUE)
  tag_expression <- FetchData(scObject, vars = tags,layer= "counts")
  tag_expression=quality_control(tag_expression)
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
  cell_clusters=scObject@meta.data[,groupby]
  names(cell_clusters)=colnames(scObject)
  cell_clusters=droplevels(cell_clusters)
  res[["tag_matrix"]]=tag_expression
  res[["cell_type"]]=cell_clusters
  niche=apply(tag_expression,2,function(x){x[x>0]})
  niche.size=sapply(niche,function(x){length(x)})
  res[["niche_tags"]]=apply(tag_expression,2,sum)
  res[["niche_size"]]=niche.size
  niche.celltypes=sapply(niche,function(x){length(unique(cell_clusters[names(x)]))})
  res[["niche_celltypes"]]=niche.celltypes
  niche_pairs=combn(names(niche),2)
  niche_share=data.frame()
  for(n in 1:ncol(niche_pairs)){
    n1=niche_pairs[1,n]
    n2=niche_pairs[2,n]
    novlp=intersect(names(niche[[n1]]),names(niche[[n2]]))
    dfovlp=cbind(niche[[n1]][novlp],niche[[n2]][novlp])
    niche_share=rbind(niche_share,df2ovlp=data.frame(niche1=n1,niche2=n2,value=sum(apply(dfovlp,1,min))))
  }
  rownames(niche_share)=NULL
  niche_share=niche_share[niche_share$value>0,]
  res[["niche_interaction"]]=niche_share
  cell.typeBarcode<-function(cell_clusters){
    cell.cluster=droplevels(cell_clusters)
    cell.cluster.code=as.data.frame(matrix(rep(0,length(cell.cluster)*length(unique(cell.cluster))),
                                           ncol=length(unique(cell.cluster))))
    colnames(cell.cluster.code)=sort(unique(cell.cluster))
    rownames(cell.cluster.code)=names(cell_clusters)
    for(cc in colnames(cell.cluster.code)){
      index=which(cell.cluster==cc)
      cell.cluster.code[index,cc]=1
    }
    cell.typeBarcode=apply(cell.cluster.code,1,function(x){paste0(x,collapse="")})
    res=list()
    res[["typeBarcode"]]=cell.typeBarcode
    res[["celltype"]]=colnames(cell.cluster.code)
    res
  }
  cell.typeBarcode.list=cell.typeBarcode(cell_clusters)
  cell.typeBarcode=cell.typeBarcode.list[["typeBarcode"]]
  res[["code2celltype"]]=cell.typeBarcode.list[["celltype"]]
  cell.setBarcode=apply(tag_expression,1,function(x){paste0(as.integer(x>0),collapse="")})
  res[["code2tag"]]=colnames(tag_expression)
  cell.2barcode=paste(cell.setBarcode,cell.typeBarcode,sep="_")
  res[["cell_barcode"]]=cell.2barcode
  set=list()
  for(setBarcode in unique(cell.2barcode)){
    set[[setBarcode]]=tag_expression[cell.2barcode==setBarcode,]
  }
  set.size=sapply(set,function(x){nrow(x)})
  res[["set"]]=set
  res[["set_size"]]=set.size
  set.mean=data.frame(t(sapply(set,function(x){apply(x,2,mean)})))
  set.mean.barcode.cell=apply(set.mean,2,function(x){sort(x[x>0],decreasing = T)})
  res[["niche"]]=set.mean.barcode.cell
  if(direction==T){
    overlap=lapply(set.mean.barcode.cell,findoverlap)
    sr2value=data.frame()
    for(ol in names(overlap)){
      if(!all(is.na(overlap[[ol]]))){
        sr2value=rbind(sr2value,data.frame(sr=names(overlap[[ol]]),value=overlap[[ol]],tag=ol))
      }
    }
    rownames(sr2value)=NULL
    shareData=data.frame(setPair=sr2value$sr,share=sr2value$value)
  }

  if(direction==F){
    print(direction)
    print(share_method)
    if(share_method=="min"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,min)}else{return(NA)}
      }))
    }else if(share_method=="mean"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,mean)}else{return(NA)}
      }))
    }else if(share_method=="max"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,max)}else{return(NA)}
      }))
    }else{
      stop("share_method wrong")
    }
    shareset=unlist(lapply(set.mean.barcode.cell,function(a){if(length(a)>1){apply(combn(names(a),2),2,function(x){paste(sort(x),collapse="-")})}else{return(NA)}}))
    shareNumber=shareNumber[!is.na(shareNumber)]
    shareset=shareset[!is.na(shareset)]
    shareData=data.frame(setPair=shareset,share=shareNumber)
  }

  shareData2=aggregate(share ~ setPair, data = shareData, sum)
  res[["set_interaction"]]=shareData2
  sharesetAB=as.data.frame(str_split_fixed(shareData2$setPair,"-",n=2))
  g <- graph_from_data_frame(sharesetAB, directed = direction)
  setID=1:length(V(g)$name)
  names(setID)=V(g)$name
  res[["setID"]]=setID
  return(res)
}

#' Calculate niche-niche connections with connectome data from nichetag method
#'
#' This fuction use tag expression matrix to infer the connection among niches labeled with different nichetags
#'
#' @param scObject    Seurat, input seurat object in which tag expression matrix is merged in gene expression matrix
#' @param groupby    character, a colname in meta data of seurat object, used for set defination
#' @param share_method    min,max,mean, mothed to calculate connection strength between different sets
#' @param direction    TorF, if or not to distinguish tag-sender or tag-receiver
#'
#' @return    a list contains connectome information of all niches:
#' \describe{
#'   \item{tag_matrix}{data.frame, tag expression matrix of all cells}
#'   \item{cell_type}{factor, cell type or state used for set defination}
#'   \item{cell_tags}{numeric, total expression of each tag}
#'   \item{niche_size}{integer, covered cell number of each tag}
#'   \item{niche_celltypes}{integer, covered cell types or cell states of each tag}
#'   \item{niche_interaction}{data.frame, overlapped cell numbers of tag-tag pairs}
#'   \item{code2celltype}{character,corresponding celltype of each celltypes code positons}
#'   \item{code2tag}{character, corresponding tag of each tag code positons}
#'   \item{cell_barcode}{character, code of each cell, tagcode+celltypecode}
#'   \item{set}{list, each element represent a tag expression matrix of the set}
#'   \item{set_size}{integer, covered cell number of each set}
#'   \item{niche}{list, mean tag expression of covered cells of all covered set of each tag}
#'   \item{set_interaction}{data.frame, overlapped tag expression of set-set pairs}
#'   \item{setID}{integer, order of sets when make graph with graph_from_data_frame in igraph }
#' }
#' @importFrom SeuratObject FetchData
#' @importFrom stats complete.cases sd ecdf
#' @importFrom igraph graph_from_data_frame
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' summary(nnt)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters",share_method="mean")
#' dnnt=Dnichenetwork(scObject,groupby="cell_clusters",direction = TRUE)
Dnichenetwork_v2<-function(scObject,groupby="cell_clusters",share_method="min",direction=FALSE){
  res=list()
  tags <- grep("^[ATCG]{8}$", rownames(scObject), value = TRUE)
  tag_expression <- FetchData(scObject, vars = tags,layer= "counts")
  tag_expression=quality_control(tag_expression)
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
  cell_clusters=scObject@meta.data[,groupby]
  names(cell_clusters)=colnames(scObject)
  cell_clusters=droplevels(cell_clusters)
  res[["tag_matrix"]]=tag_expression
  res[["cell_type"]]=cell_clusters
  niche=apply(tag_expression,2,function(x){x[x>0]})
  niche.size=sapply(niche,function(x){length(x)})
  res[["niche_tags"]]=apply(tag_expression,2,sum)
  res[["niche_size"]]=niche.size
  niche.celltypes=sapply(niche,function(x){length(unique(cell_clusters[names(x)]))})
  res[["niche_celltypes"]]=niche.celltypes
  niche_pairs=combn(names(niche),2)
  niche_share=data.frame()
  for(n in 1:ncol(niche_pairs)){
    n1=niche_pairs[1,n]
    n2=niche_pairs[2,n]
    novlp=intersect(names(niche[[n1]]),names(niche[[n2]]))
    dfovlp=cbind(niche[[n1]][novlp],niche[[n2]][novlp])
    niche_share=rbind(niche_share,df2ovlp=data.frame(niche1=n1,niche2=n2,value=sum(apply(dfovlp,1,min))))
  }
  rownames(niche_share)=NULL
  niche_share=niche_share[niche_share$value>0,]
  res[["niche_interaction"]]=niche_share
  cell.typeBarcode<-function(cell_clusters){
    cell.cluster=droplevels(cell_clusters)
    cell.cluster.code=as.data.frame(matrix(rep(0,length(cell.cluster)*length(unique(cell.cluster))),
                                           ncol=length(unique(cell.cluster))))
    colnames(cell.cluster.code)=sort(unique(cell.cluster))
    rownames(cell.cluster.code)=names(cell_clusters)
    for(cc in colnames(cell.cluster.code)){
      index=which(cell.cluster==cc)
      cell.cluster.code[index,cc]=1
    }
    cell.typeBarcode=apply(cell.cluster.code,1,function(x){paste0(x,collapse="")})
    res=list()
    res[["typeBarcode"]]=cell.typeBarcode
    res[["celltype"]]=colnames(cell.cluster.code)
    res
  }
  cell.typeBarcode.list=cell.typeBarcode(cell_clusters)
  cell.typeBarcode=cell.typeBarcode.list[["typeBarcode"]]
  res[["code2celltype"]]=cell.typeBarcode.list[["celltype"]]
  cell.setBarcode=apply(tag_expression,1,function(x){paste0(as.integer(x>0),collapse="")})
  res[["code2tag"]]=colnames(tag_expression)
  cell.2barcode=paste(cell.setBarcode,cell.typeBarcode,sep="_")
  res[["cell_barcode"]]=cell.2barcode
  set=list()
  for(setBarcode in unique(cell.2barcode)){
    set[[setBarcode]]=tag_expression[cell.2barcode==setBarcode,]
  }
  set.size=sort(sapply(set,function(x){nrow(x)}),decreasing = T)#correction1
  res[["set"]]=set
  res[["set_size"]]=set.size
  tag_set=names(set.size)[is.set(names(set.size))]#correction2
  tag_set_ids=1:length(tag_set)#correction2
  names(tag_set_ids)=tag_set#correction2
  res[["setID"]]=tag_set_ids#correction2
  set.mean=data.frame(t(sapply(set,function(x){apply(x,2,mean)})))
  set.mean.barcode.cell=apply(set.mean,2,function(x){sort(x[x>0],decreasing = T)})
  res[["niche"]]=set.mean.barcode.cell
  if(direction==T){
    overlap=lapply(set.mean.barcode.cell,findoverlap)
    sr2value=data.frame()
    for(ol in names(overlap)){
      if(!all(is.na(overlap[[ol]]))){
        sr2value=rbind(sr2value,data.frame(sr=names(overlap[[ol]]),value=overlap[[ol]],tag=ol))
      }
    }
    rownames(sr2value)=NULL
    shareData=data.frame(setPair=sr2value$sr,share=sr2value$value)
  }

  if(direction==F){
    print(direction)
    print(share_method)
    if(share_method=="min"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,min)}else{return(NA)}
      }))
    }else if(share_method=="mean"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,mean)}else{return(NA)}
      }))
    }else if(share_method=="max"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,max)}else{return(NA)}
      }))
    }else{
      stop("share_method wrong")
    }
    shareset=unlist(lapply(set.mean.barcode.cell,function(a){if(length(a)>1){apply(combn(names(a),2),2,function(x){paste(sort(x),collapse="-")})}else{return(NA)}}))
    shareNumber=shareNumber[!is.na(shareNumber)]
    shareset=shareset[!is.na(shareset)]
    shareData=data.frame(setPair=shareset,share=shareNumber)
  }

  shareData2=aggregate(share ~ setPair, data = shareData, sum)
  res[["set_interaction"]]=shareData2
  # sharesetAB=as.data.frame(str_split_fixed(shareData2$setPair,"-",n=2))
  # g <- graph_from_data_frame(sharesetAB, directed = direction)
  # setID=1:length(V(g)$name)
  # names(setID)=V(g)$name
  # res[["setID"]]=setID
  return(res)
}

#' Calculate niche-niche connections with connectome data from nichetag method
#'
#' This fuction use tag expression matrix to infer the connection among niches labeled with different nichetags
#'
#' @param scObject    Seurat, input seurat object in which tag expression matrix is merged in gene expression matrix
#' @param groupby    character, a colname in meta data of seurat object, used for set defination
#' @param share_method    string, min,max,mean, the method to calculate connection strength between different sets
#' @param highconfidence    TorF, if or not to distinguish high confident tag-sender
#' @param send_cutoff1    int, count cutoff for primary celltag asignment
#' @param send_cutoff2    numeric, percentage cutoff for primary celltag asignment
#' @param direction    TorF, if or not to distinguish tag-sender or tag-receiver
#'
#' @return    a list contains connectome information of all niches:
#' \describe{
#'   \item{tag_matrix}{data.frame, tag expression matrix of all cells}
#'   \item{cell_type}{factor, cell type or state used for set defination}
#'   \item{cell_tags}{numeric, total expression of each tag}
#'   \item{niche_size}{integer, covered cell number of each tag}
#'   \item{niche_celltypes}{integer, covered cell types or cell states of each tag}
#'   \item{niche_interaction}{data.frame, overlapped cell numbers of tag-tag pairs}
#'   \item{code2celltype}{character,corresponding celltype of each celltypes code positons}
#'   \item{code2tag}{character, corresponding tag of each tag code positons}
#'   \item{cell_barcode}{character, code of each cell, tagcode+celltypecode}
#'   \item{set}{list, each element represent a tag expression matrix of the set}
#'   \item{set_size}{integer, covered cell number of each set}
#'   \item{niche}{list, mean tag expression of covered cells of all covered set of each tag}
#'   \item{set_interaction}{data.frame, overlapped tag expression of set-set pairs}
#'   \item{setID}{integer, order of sets when make graph with graph_from_data_frame in igraph }
#' }
#' @importFrom SeuratObject FetchData
#' @importFrom stats complete.cases sd ecdf
#' @importFrom igraph graph_from_data_frame
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' summary(nnt)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters",share_method="mean")
#' dnnt=Dnichenetwork(scObject,groupby="cell_clusters",direction = TRUE)
Dnichenetwork<-function(scObject,groupby="cell_clusters",share_method="min",
                        highconfidence=TRUE,send_cutoff1=10,send_cutoff2=0.8,
                        direction=FALSE){
  res=list()
  tags <- grep("^[ATCG]{8}$", rownames(scObject), value = TRUE)
  tag_expression <- FetchData(scObject, vars = tags,layer= "counts")
  tag_expression=quality_control(tag_expression)
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
  cell_clusters=scObject@meta.data[,groupby]
  names(cell_clusters)=colnames(scObject)
  cell_clusters=droplevels(cell_clusters)
  res[["tag_matrix"]]=tag_expression
  res[["cell_type"]]=cell_clusters
  niche=apply(tag_expression,2,function(x){x[x>0]})
  niche.size=sapply(niche,function(x){length(x)})
  res[["niche_tags"]]=apply(tag_expression,2,sum)
  res[["niche_size"]]=niche.size
  niche.celltypes=sapply(niche,function(x){length(unique(cell_clusters[names(x)]))})
  res[["niche_celltypes"]]=niche.celltypes
  niche_pairs=combn(names(niche),2)
  niche_share=data.frame()
  for(n in 1:ncol(niche_pairs)){
    n1=niche_pairs[1,n]
    n2=niche_pairs[2,n]
    novlp=intersect(names(niche[[n1]]),names(niche[[n2]]))
    dfovlp=cbind(niche[[n1]][novlp],niche[[n2]][novlp])
    niche_share=rbind(niche_share,df2ovlp=data.frame(niche1=n1,niche2=n2,value=sum(apply(dfovlp,1,min))))
  }
  rownames(niche_share)=NULL
  niche_share=niche_share[niche_share$value>0,]
  res[["niche_interaction"]]=niche_share
  cell.typeBarcode<-function(cell_clusters){
    cell.cluster=droplevels(cell_clusters)
    cell.cluster.code=as.data.frame(matrix(rep(0,length(cell.cluster)*length(unique(cell.cluster))),
                                           ncol=length(unique(cell.cluster))))
    colnames(cell.cluster.code)=sort(unique(cell.cluster))
    rownames(cell.cluster.code)=names(cell_clusters)
    for(cc in colnames(cell.cluster.code)){
      index=which(cell.cluster==cc)
      cell.cluster.code[index,cc]=1
    }
    cell.typeBarcode=apply(cell.cluster.code,1,function(x){paste0(x,collapse="")})
    res=list()
    res[["typeBarcode"]]=cell.typeBarcode
    res[["celltype"]]=colnames(cell.cluster.code)
    res
  }
  cell.typeBarcode.list=cell.typeBarcode(cell_clusters)
  cell.typeBarcode=cell.typeBarcode.list[["typeBarcode"]]
  res[["code2celltype"]]=cell.typeBarcode.list[["celltype"]]
  if(highconfidence){
    cell.setBarcode=apply(tag_expression,1,function(x){
      cb=paste0(as.integer(x>0),collapse="")
      if((max(x)>=send_cutoff1)&(max(x)/sum(x)>=send_cutoff2)){
        position=order(x,decreasing = T)[1]
        substr(cb,position,position)="2"
      }
      cb
    })
  }else{
    cell.setBarcode=apply(tag_expression,1,function(x){paste0(as.integer(x>0),collapse="")})
  }
  res[["code2tag"]]=colnames(tag_expression)
  cell.2barcode=paste(cell.setBarcode,cell.typeBarcode,sep="_")

  res[["cell_barcode"]]=cell.2barcode
  set=list()
  for(setBarcode in unique(cell.2barcode)){
    set[[setBarcode]]=tag_expression[cell.2barcode==setBarcode,]
  }
  set.size=sort(sapply(set,function(x){nrow(x)}),decreasing = T)#correction1
  res[["set"]]=set
  res[["set_size"]]=set.size
  tag_set=names(set.size)[is.set(names(set.size))]#correction2
  tag_set_ids=1:length(tag_set)#correction2
  names(tag_set_ids)=tag_set#correction2
  res[["setID"]]=tag_set_ids#correction2
  set.mean=data.frame(t(sapply(set,function(x){apply(x,2,mean)})))
  set.mean.barcode.cell=apply(set.mean,2,function(x){sort(x[x>0],decreasing = T)})
  res[["niche"]]=set.mean.barcode.cell
  if(direction==T){
    overlap=lapply(set.mean.barcode.cell,findoverlap)
    sr2value=data.frame()
    for(ol in names(overlap)){
      if(!all(is.na(overlap[[ol]]))){
        sr2value=rbind(sr2value,data.frame(sr=names(overlap[[ol]]),value=overlap[[ol]],tag=ol))
      }
    }
    rownames(sr2value)=NULL
    shareData=data.frame(setPair=sr2value$sr,share=sr2value$value)
  }

  if(direction==F){
    print(direction)
    print(share_method)
    if(share_method=="min"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,min)}else{return(NA)}
      }))
    }else if(share_method=="mean"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,mean)}else{return(NA)}
      }))
    }else if(share_method=="max"){
      shareNumber=unlist(lapply(set.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,max)}else{return(NA)}
      }))
    }else{
      stop("share_method wrong")
    }
    shareset=unlist(lapply(set.mean.barcode.cell,function(a){if(length(a)>1){apply(combn(names(a),2),2,function(x){paste(sort(x),collapse="-")})}else{return(NA)}}))
    shareNumber=shareNumber[!is.na(shareNumber)]
    shareset=shareset[!is.na(shareset)]
    shareData=data.frame(setPair=shareset,share=shareNumber)
  }

  shareData2=aggregate(share ~ setPair, data = shareData, sum)
  res[["set_interaction"]]=shareData2
  # sharesetAB=as.data.frame(str_split_fixed(shareData2$setPair,"-",n=2))
  # g <- graph_from_data_frame(sharesetAB, directed = direction)
  # setID=1:length(V(g)$name)
  # names(setID)=V(g)$name
  # res[["setID"]]=setID
  return(res)
}


#' find overlap of sets
#'
#' @param settag     sets
#' @param share_method    min,max,mean, mothed to calculate connection strength between different sets
#'
#' @return    overlap
#' @noRd
findoverlap<-function(settag,share_method="min"){
  settag=sort(settag,decreasing = T)
  sender=names(settag)[settag>=10]
  receivor=setdiff(names(settag),sender)
  if(length(sender)>0 & length(receivor)>0){
    pairs=combn(names(settag),2)
    setaside=which(pairs[1,] %in% sender)
    pairs2=pairs[,setaside,drop=F]
    print(dim(pairs))
    print(dim(pairs2))
    pairstag=rbind(settag[pairs2[1,]],settag[pairs2[2,]])
    if(share_method=="min"){shareNumber=apply(pairstag,2,min)}
    if(share_method=="mean"){shareNumber=apply(pairstag,2,mean)}
    if(share_method=="max"){shareNumber=apply(pairstag,2,max)}
    names(shareNumber)=apply(pairs2,2,function(x){paste(x,collapse="-")})
    return(shareNumber)
  }else{
    return(NA)
  }
}

#' Extract set connections from nnt for Gephi analysis
#'
#' @param nnt a list contains connectome information of all niches
#' @param file output file name
#'
#' @returns a file for Gephi to get more circular layout
#' @importFrom stringr str_split_fixed
#' @import openxlsx
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' Gephi_prepare(nnt)
Gephi_prepare<-function(nnt,file="Gephi_input.xlsx"){
  id.df=str_split_fixed(nnt$set_interaction$setPair,"-",n=2)
  id.df2=apply(id.df,2,function(x){nnt$setID[x]})
  res=cbind(id.df2,nnt$set_interaction$share)
  colnames(res)=c("source","target","weight")
  wb <- createWorkbook()
  addWorksheet(wb, "edges")
  writeData(wb, "edges", res,rowNames = F)
  addWorksheet(wb, "nodes")
  writeData(wb, "nodes", data.frame(Id=nnt$setID),rowNames = F)
  saveWorkbook(wb, file, overwrite = TRUE)
}


#' Get set expression matrix, for each gene use mean value of all cells of the same set
#'
#' @param scObject Seurat, input seurat object in which tag expression matrix is merged in gene expression matrix
#' @param nnt a list contains connectome information of all niches
#' @param seurat_layer layer in scObject used to calculate expression level of genes
#'
#' @returns set expression matrix of all nontag genes
#' @importFrom SeuratObject FetchData
#' @export
#'
#' @examples
#' data(example)
#' nnt <- Dnichenetwork(scObject, groupby = "cell_clusters")
#' setmatrix=set_expr(scObject,nnt,seurat_layer="counts")
set_expr<-function(scObject,nnt,seurat_layer="counts"){
  nontags <- grep("^[ATCG]{8}$", rownames(scObject), value = TRUE,invert = TRUE)
  expr=FetchData(scObject, vars = nontags,layer= seurat_layer)
  res=sapply(nnt[["set"]],function(x){
    cells=rownames(x)
    apply(expr[cells,],2,mean)
  })
  res=as.data.frame(res)
  res=res[,colnames(res) %in% names(nnt$setID)]
  colnames(res)=paste0("set",nnt[["setID"]][colnames(res)])
  res
}


#' Get logic vector, whether input sets express tags
#'
#' @param code vector, set codes
#'
#' @returns TorF logic vector showing whether sets express tags
#' @importFrom stringr str_split_fixed
#' @noRd
is.set<-function(code){
  cell.setBarcode=str_split_fixed(code,"_",n=2)[,1]
  code_matrix=sapply(strsplit(cell.setBarcode, split = ""),as.integer)
  apply(code_matrix,2,sum)>0
}


