#' print niche network
#'
#' @param nnt    a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file    file names to be saved
#' @param vsize    vertex size
#' @param esize    edge size
#' @param vertex.label.cex    vertex label size
#' @param axes   TorF, if or not draw axes
#' @param width    file size
#' @param height    file size
#' @param weighted    layout caculation parameter
#' @param direction    layout caculation parameter
#'
#' @return    pdf file
#' @import igraph
#' @importFrom scales hue_pal
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' print_nichenet(nnt,file="nichenet-w2.pdf",vsize = 10)
print_nichenet<-function(nnt,file="nichenet.pdf",vsize=1,esize=1,
                         vertex.label.cex =1, axes=FALSE,
                         width=7,height=7,weighted=FALSE,direction=FALSE){

  niche.size=nnt[["niche_size"]]
  niche_share=nnt[["niche_interaction"]]
  g <- graph_from_data_frame(niche_share[,1:2], directed = direction)
  # 设置顶点大小
  V(g)$size <- log10(1+niche.size[V(g)$name])  # 通过节点名称匹配大小
  E(g)$size <- log10(1+niche_share[,3])

  nicheID=1:length(V(g)$name)
  names(nicheID)=V(g)$name

  #library(matlab)
  #color_pallete=jet.colors(length(V(g)$name))
  #library(scales)
  color_pallete=hue_pal()(length(V(g)$name))

  pdf(file,width = width,height = height)
  #E(g)$curved <- seq(-0.5, 0.5, length.out=ecount(g))
  #E(g)$curved <- rep( 1, ecount(g))
  # 绘制网络图
  set.seed(666)
  if(weighted){
    lo=layout_with_fr(g,  grid = "nogrid", niter = 100000, weights = E(g)$size)
    #print(lo)
    #print(norm_coords(lo))
  }else{
    lo=layout_with_fr(g, grid = "nogrid", niter = 100000)
  }
  plot(g,
       axes = axes,
       layout = lo,
       vertex.frame.width = 0.5,
       vertex.label = V(g)$name,
       vertex.label.cex = vertex.label.cex,
       edge.arrow.size = 0.2,
       edge.width = esize*E(g)$size,  # 使用 size 作为边的宽度
       vertex.size = vsize*V(g)$size, # 使用 vector 作为顶点大小
       edge.color = "grey",
       vertex.color = color_pallete)
  dev.off()
}

#' print set network
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf file name to be saved
#' @param file_layout file name to save layout of sets
#' @param ilayout file used to set external layout
#' @param highlight whether to hightlight sets for each tag in pdf file
#' @param seed seed for random set colour setting
#' @param mark_groups igrap parameter
#' @param margin igraph parameter
#' @param vertex.label.cex vertex label size
#' @param vsize vertex size
#' @param esize edge size
#' @param ecolor edge color
#' @param width file size
#' @param height file size
#' @param weighted layout parameter
#' @param direction layout parameter
#' @param niter layout parameter
#' @param axes TorF, if or not draw axes
#'
#' @return df.niche2setID
#' @import igraph
#' @importFrom Polychrome createPalette
#' @importFrom scales hue_pal
#' @importFrom stringr str_split_fixed str_split
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' print_nichenetwork(nnt,file="nichenetwork-w3-26.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=adjustcolor("grey80", alpha.f = 0.5),mark_groups = TRUE)
#' print_nichenetwork(nnt,file="nichenetwork-w3-27.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=NA,mark_groups = TRUE)
#' print_nichenetwork(nnt,file="nichenetwork-w3-28.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=NA)
#' print_nichenetwork(nnt,file="nichenetwork-w3-29.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1)
print_nichenetwork_v1<-function(nnt,file="nichenetwork.pdf",
                             file_layout="layout.csv",ilayout=NULL,
                             highlight=FALSE,seed=888,
                             mark_groups = FALSE,margin = c(0, 0, 0, 0),
                             vertex.label.cex=1,vsize=1,esize=1,ecolor="grey80",
                             width=7,height=7,
                             weighted=FALSE,direction=FALSE,niter = 100000,
                             axes=FALSE){
  shareData2=nnt[["set_interaction"]]
  set.size=nnt[["set_size"]]
  celltypes=nnt[["code2celltype"]]
  niche=nnt[["niche"]]
  #将互作set分开
  sharesetAB=as.data.frame(str_split_fixed(shareData2$setPair,"-",n=2))
  #shareData3=data.frame(c1=sharesetAB[,1],c2=sharesetAB[,2],edge)
  #取出所有包含交叉的set,无互作的克隆不在igraph图中显示
  #sharesetID=unique(c(sharesetAB[,1],sharesetAB[,2]))#correction2
  #包含交叉的set的size用其中的细胞数表示
  #sharesetSize=set.size[sharesetID]#correction2
  #将sharesetID，sharesetSize等作为输出结果输出到list
  #使用igraph画set互作network，input是sharesetAB，shareData2$share，sharesetSize
  # 创建图对象
  g <- graph_from_data_frame(d=sharesetAB, vertices=names(nnt[["setID"]]),directed = direction)#correction1
  if(identical(V(g)$name,names(nnt[["setID"]]))){
     print("Identical")
     setID=nnt[["setID"]]
   }else{
     setID=nnt[["setID"]][V(g)$name]
   }
  #setID=1:length(V(g)$name)
  #names(setID)=V(g)$name
  #niche2setID=lapply(niche,function(x){setID[intersect(names(x),names(setID))]})
  #newniche=lapply(niche,function(x){x[intersect(names(x),names(setID))]})
  # 设置顶点大小
  V(g)$size <- vsize*log10(1+set.size[V(g)$name])  # 通过节点名称匹配大小correction3
  #E(g)$size <- log2(1+shareData2$share)
  E(g)$size <- shareData2$share

  pdf(file,width = width,height = height)
  last_code<- str_split_fixed(V(g)$name, "_",n=2)[,2]
  set.seed(seed)
  color_pallete <- createPalette(nchar(last_code[1]),
                                 seedcolors = c("#E69F00" ,"#56B4E9" ,"#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7"))
  char_vectors <- strsplit(last_code, split = "")
  #print(char_vectors)
  set_df=sapply(char_vectors,as.integer)
  vcolors=apply(set_df,2,function(x){color_pallete[as.logical(x)]})
  #set_cluster=apply(set_df,2,function(x){nnt[["code2celltype"]][as.logical(x)]})
  # 绘制网络图
  if(is.null(ilayout)){
    set.seed(seed)
    if(weighted){
      lo=layout_with_fr(g,  grid = "nogrid", niter = niter, weights = E(g)$size)
      #print(lo)
      #print(norm_coords(lo))
    }else{
      lo=layout_with_fr(g, grid = "nogrid", niter = niter)
    }
  }else{
    lo=as.matrix(ilayout)
  }
  write.csv(lo,file = file_layout,row.names=F)
  if(mark_groups){
    #library(Polychrome)
    first_code <-str_split_fixed(V(g)$name, "_",n=2)[,1]
    char_vectors <- strsplit(first_code, split = "")
    set_df=as.data.frame(sapply(char_vectors,as.integer))
    colnames(set_df)=1:ncol(set_df)
    set_1tag=set_df[,which(apply(set_df,2,sum)==1)]
    set_tag=apply(set_1tag,2,function(x){nnt[["code2tag"]][as.logical(x)]})
    set1tag_list=list()
    for(tag in nnt[["code2tag"]]){
      set1tag_list[[tag]]=as.integer(names(set_tag[set_tag==tag]))
    }
    mark.groups=set1tag_list[sapply(set1tag_list,length)>0]
    #library(scales)
    mark.col=hue_pal()(length(set1tag_list))
  }else{
    mark.groups = list()
    mark.col = rainbow(length(mark.groups), alpha = 0.3)
  }

  plot(g,
       axes = axes,
       layout = lo,
       mark.groups = mark.groups,
       mark.col=mark.col,
       #mark.border=mark.colors,
       mark.shape =1/2,
       vertex.frame.width = 0.5,
       vertex.label = setID,
       vertex.label.cex = vertex.label.cex,
       edge.arrow.size = 0.2,
       edge.width = esize*log10(1+E(g)$size),  # 使用 size 作为边的宽度
       vertex.size = V(g)$size, # 使用 vector 作为顶点大小
       edge.color = ecolor,
       vertex.color = vcolors,
       margin = margin
  )

  legend("topleft", legend = celltypes,col = color_pallete,pch = 21, pt.bg = color_pallete, pt.cex = 1,cex = 0.5,bty = "n")

  if(mark_groups){legend("topright", legend = names(set1tag_list),col = mark.col,
                         pch = 21, pt.bg = mark.col, pt.cex = 1,cex = 0.5,bty = "n")}
  if(highlight==T){
    nichetags=names(nnt$niche)[sapply(nnt$niche,length)>0]
    for(tag in nichetags){
      vcolors2=vcolors
      vcolors2[!(V(g)$name %in% names(nnt[["niche"]][[tag]]))]="grey80"
      subtitle=paste(setID[(V(g)$name %in% names(nnt[["niche"]][[tag]]))],collapse=",")
      #print(tag)
      #print(sum(V(g)$name %in% names(nnt[["niche"]][[tag]])))
      plot(g,
	   main=tag,sub=subtitle,
           axes = axes,
           layout = lo,
           mark.groups = mark.groups,
           mark.col=mark.col,
           #mark.border=mark.colors,
           mark.shape =1/2,
           vertex.frame.width = 0.5,
           vertex.label = setID,
           vertex.label.cex = vertex.label.cex,
           edge.arrow.size = 0.2,
           edge.width = esize*log10(1+E(g)$size),  # 使用 size 作为边的宽度
           vertex.size = V(g)$size, # 使用 vector 作为顶点大小
           edge.color = ecolor,
           vertex.color = vcolors2,
           margin = margin
      )
      legend("topleft", legend = celltypes,col = color_pallete,pch = 21, pt.bg = color_pallete, pt.cex = 1,cex = 0.5,bty = "n")
    }
  }
  dev.off()

  #df.niche2setID=data.frame()
  #for(name in names(niche2setID)){
   # if(length(newniche[[name]])>0){
    #  df=data.frame(tag=name,tagnum=newniche[[name]],
     #               setID=niche2setID[[name]],
      #              code=names(niche2setID[[name]]))
      #df.niche2setID=rbind(df.niche2setID,df)
    #}
  #}
  #rownames(df.niche2setID)=NULL
  #return(lo)
}

#' print set network
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf file name to be saved
#' @param file_layout file name to save layout of sets
#' @param ilayout file used to set external layout
#' @param vertex2label TorF, if or not lable vertex
#' @param edge2label TorF, if or not lable edge weight
#' @param highlight.color string, epxr,sr,num, an attribute of the set to hightlight for each tag in pdf file
#' @param seed seed for random set colour setting
#' @param mark_groups igrap parameter
#' @param margin igraph parameter
#' @param vertex.label.cex vertex label size
#' @param edge.label.cex edge label size
#' @param vsize vertex size
#' @param esize edge size
#' @param ecolor edge color
#' @param width file size
#' @param height file size
#' @param weighted layout parameter
#' @param direction layout parameter
#' @param niter layout parameter
#' @param axes TorF, if or not draw axes
#'
#' @return pdf
#' @import igraph
#' @importFrom Polychrome createPalette
#' @importFrom scales hue_pal rescale
#' @importFrom viridis viridis
#' @importFrom stringr str_split_fixed str_split str_count
#' @importFrom fields image.plot
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' print_nichenetwork(nnt,file="nichenetwork-w3-26.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=adjustcolor("grey80", alpha.f = 0.5),mark_groups = TRUE)
#' print_nichenetwork(nnt,file="nichenetwork-w3-27.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=NA,mark_groups = TRUE)
#' print_nichenetwork(nnt,file="nichenetwork-w3-28.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=NA)
#' print_nichenetwork(nnt,file="nichenetwork-w3-29.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1)
print_nichenetwork<-function(nnt,file="nichenetwork.pdf",
                             file_layout="layout.csv",ilayout=NULL,vertex2label=TRUE,edge2label=FALSE,
                             highlight.color=NULL,seed=123,
                             mark_groups = FALSE,margin = c(0, 0, 0, 0),
                             vertex.label.cex=1,edge.label.cex=1,vsize=1,esize=1,ecolor="grey80",
                             width=7,height=7,
                             weighted=FALSE,direction=FALSE,niter = 100000,
                             axes=FALSE){
  shareData2=nnt[["set_interaction"]]
  set.size=nnt[["set_size"]]
  celltypes=nnt[["code2celltype"]]
  niche=nnt[["niche"]]
  #将互作set分开
  sharesetAB=as.data.frame(str_split_fixed(shareData2$setPair,"-",n=2))
  ishighconfident_c1=str_detect(sharesetAB[,1],"2")
  iscancer_c1=str_detect(codesplit(sharesetAB[,1],nnt$code2tag,nnt$code2celltype)$types,"Cancer")
  ishighconfident_c2=str_detect(sharesetAB[,2],"2")
  iscancer_c2=str_detect(codesplit(sharesetAB[,2],nnt$code2tag,nnt$code2celltype)$types,"Cancer")
  keep_row= !(iscancer_c1 & (!ishighconfident_c1)) & !(iscancer_c2 & (!ishighconfident_c2))
  #shareData3=data.frame(c1=sharesetAB[,1],c2=sharesetAB[,2],edge)
  #取出所有包含交叉的set,无互作的克隆不在igraph图中显示
  #sharesetID=unique(c(sharesetAB[,1],sharesetAB[,2]))#correction2
  #包含交叉的set的size用其中的细胞数表示
  #sharesetSize=set.size[sharesetID]#correction2
  #将sharesetID，sharesetSize等作为输出结果输出到list
  #使用igraph画set互作network，input是sharesetAB，shareData2$share，sharesetSize
  # 创建图对象
  g <- graph_from_data_frame(d=sharesetAB[keep_row,],directed = direction)#correction1
  if(identical(V(g)$name,names(nnt[["setID"]]))){
    print("Identical")
    setID=nnt[["setID"]]
  }else{
    setID=nnt[["setID"]][V(g)$name]
  }
  #setID=1:length(V(g)$name)
  #names(setID)=V(g)$name
  #niche2setID=lapply(niche,function(x){setID[intersect(names(x),names(setID))]})
  #newniche=lapply(niche,function(x){x[intersect(names(x),names(setID))]})
  # 设置顶点大小
  V(g)$size <- vsize*log10(1+set.size[V(g)$name])  # 通过节点名称匹配大小correction3
  #E(g)$size <- log2(1+shareData2$share)
  E(g)$size <- shareData2[keep_row,]$share

  pdf(file,width = width,height = height)
  last_code<- str_split_fixed(V(g)$name, "_",n=2)[,2]
  set.seed(seed)
  color_pallete <- createPalette(nchar(last_code[1]),
                                 seedcolors = c("#E69F00" ,"#56B4E9" ,"#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7"))
  char_vectors <- strsplit(last_code, split = "")
  #print(char_vectors)
  set_df=sapply(char_vectors,as.integer)
  vcolors=apply(set_df,2,function(x){color_pallete[as.logical(x)]})
  #set_cluster=apply(set_df,2,function(x){nnt[["code2celltype"]][as.logical(x)]})
  # 绘制网络图
  if(is.null(ilayout)){
    set.seed(seed)
    if(weighted){
      lo=layout_with_fr(g,  grid = "nogrid", niter = niter, weights = E(g)$size)
      #print(lo)
      #print(norm_coords(lo))
    }else{
      lo=layout_with_fr(g, grid = "nogrid", niter = niter)
    }
  }else{
    lo=as.matrix(ilayout[paste0("id_",setID),])
  }
  rownames(lo)=V(g)$name
  V(g)$x=lo[,1]
  V(g)$y=lo[,2]
  #write.csv(lo,file = file_layout,row.names=F)
  if(mark_groups){
    #library(Polychrome)
    first_code <-str_split_fixed(V(g)$name, "_",n=2)[,1]
    char_vectors <- strsplit(first_code, split = "")
    set_df=as.data.frame(sapply(char_vectors,as.integer))
    colnames(set_df)=1:ncol(set_df)
    set_1tag=set_df[,which(apply(set_df,2,sum)==1)]
    set_tag=apply(set_1tag,2,function(x){nnt[["code2tag"]][as.logical(x)]})
    set1tag_list=list()
    for(tag in nnt[["code2tag"]]){
      set1tag_list[[tag]]=as.integer(names(set_tag[set_tag==tag]))
    }
    mark.groups=set1tag_list[sapply(set1tag_list,length)>0]
    #library(scales)
    mark.col=hue_pal()(length(set1tag_list))
  }else{
    mark.groups = list()
    mark.col = rainbow(length(mark.groups), alpha = 0.3)
  }
  if(vertex2label==T){vertex.label=setID}else{vertex.label=NA}
  if(edge2label==T){edge.label=round(E(g)$size,2)}else{edge.label=NA}

  edge_df_modified=data.frame()
  nichetags=names(nnt$niche)[sapply(nnt$niche,length)>0]
  for(tag in nichetags){
    vcolors2=vcolors
    vcolors2[!(V(g)$name %in% names(nnt[["niche"]][[tag]]))]="grey80"
    if(!is.null(highlight.color)){
      if(highlight.color=="expr"){
        set_tag_mean=nnt[["niche"]][[tag]][V(g)$name]
        set_tag_mean=set_tag_mean[!is.na(set_tag_mean)]
        if(length(set_tag_mean)<=1){next}#注意
        set_tag_color=viridis(100)[rescale(log(set_tag_mean), to = c(1,100)) %>% round()]
        vcolors2[(V(g)$name %in% names(nnt[["niche"]][[tag]]))]= set_tag_color
      }
      if(highlight.color=="sr"){
        idx=which(names(nnt$niche)==tag)
        cname=intersect(V(g)$name,names(nnt[["niche"]][[tag]]))
        char_x <- substr(cname, idx, idx)
        set_tag_color=rep("#440154FF",length(cname))
        set_tag_color[char_x=="2"]="#FDE725FF"
        vcolors2[(V(g)$name %in% names(nnt[["niche"]][[tag]]))]= set_tag_color
      }
      if(highlight.color=="num"){
        cname=intersect(V(g)$name,names(nnt[["niche"]][[tag]]))
        tagcode=str_split_fixed(cname,"_",n=2)[,1]
        n_zeros <- str_count(tagcode, "0")
        n_types=length(nnt$code2tag)
        n_tags =n_types - n_zeros
        n_tags[n_tags>7]=7
        set_tag_color=viridis(7)[n_tags]
        vcolors2[(V(g)$name %in% names(nnt[["niche"]][[tag]]))]= set_tag_color
      }
    }
    print(tag)
    set_tag=names(nnt[["niche"]][[tag]])
    keep_nodes=intersect(V(g)$name,set_tag)
    keep_edges <- E(g)[
      ends(g, E(g))[,1] %in% keep_nodes &
        ends(g, E(g))[,2] %in% keep_nodes
      ]
    g2 <- subgraph.edges(
      g,
      keep_edges,
      delete.vertices = FALSE
    )
    idx=which(names(nnt$niche)==tag)
    char_x <- substr(keep_nodes, idx, idx)
    sender=keep_nodes[char_x=="2"]
    keep_edges <- incident_edges(g2, sender)
    g3 <- subgraph_from_edges(
      g2,
      unlist(keep_edges),
      delete.vertices = FALSE
    )
    edges_to_remove <- E(g3)[
      ends(g3, E(g3))[,1] %in% sender &
        ends(g3, E(g3))[,2] %in% sender
      ]
    g4 <- delete_edges(g3, edges_to_remove)
    E(g4)$size <- apply(ends(g4, E(g4)),1,function(x){min(nnt[["niche"]][[tag]][x])})
    edge_df <- igraph::as_data_frame(g4, what = "edges")
    if(nrow(edge_df)>0){
      edge_df$celltag=tag
      edge_df_modified=rbind(edge_df_modified,edge_df)
    }
    #vcolors2=vcolors
    #vcolors2[!(V(g)$name %in% names(nnt[["niche"]][[tag]]))]="grey80"
    subtitle=paste(setID[(V(g)$name %in% names(nnt[["niche"]][[tag]]))],collapse=",")
    #print(tag)
    #print(sum(V(g)$name %in% names(nnt[["niche"]][[tag]])))
    if(edge2label==T){edge.label=round(E(g4)$size,2)}else{edge.label=NA}
    plot(g4,
         main=tag,sub=subtitle,
         axes = axes,
         mark.groups = mark.groups,
         mark.col=mark.col,
         #mark.border=mark.colors,
         mark.shape =1/2,
         vertex.frame.width = 0.5,
         vertex.label = vertex.label,
         vertex.label.cex = vertex.label.cex,
         edge.arrow.size = 0.1,
         edge.width = esize*log10(1+E(g4)$size),  # 不使用 size 作为边的宽度
         edge.label = edge.label,
         edge.label.cex = edge.label.cex,
         vertex.size = V(g)$size, # 使用 vector 作为顶点大小
         edge.color = ecolor,
         vertex.color = vcolors2,
         margin = margin
    )
    if(is.null(highlight.color)){
      legend("topleft", legend = celltypes,col = color_pallete,pch = 21, pt.bg = color_pallete, pt.cex = 1,cex = 0.5,bty = "n")
    }else{
      if(highlight.color=="sr"){
        legend("topleft", legend = c("Sender","Receiver"),
               col = c("#FDE725FF","#440154FF"),pch = 21, pt.bg = c("#FDE725FF","#440154FF"), pt.cex = 1,cex = 0.5,bty = "n")
      }
      if(highlight.color=="num"){
        legend("topleft", legend = rev(1:7),
               col = rev(viridis(7)),pch = 21, pt.bg = rev(viridis(7)), pt.cex = 1,cex = 0.5,bty = "n")
      }
      if(highlight.color=="expr"){
        if(max(range(log(set_tag_mean), na.rm = TRUE))>0){
          #print(range(log(set_tag_mean), na.rm = TRUE))
          image.plot(
            legend.only = TRUE,
            zlim = range(log(set_tag_mean), na.rm = TRUE),
            col = viridis(100),
            # 左上角
            smallplot = c(0.18, 0.20, 0.68, 0.88),
            legend.mar = 2,
            axis.args = list(
              cex.axis = 0.6
            ),
            legend.args = list(
              text = "Log(mean CellTag UMI)",
              side = 4,
              line = -2,
              cex = 0.7
            )
          )
        }
      }
    }
  }
  shareData=data.frame(setPair=paste(edge_df_modified[,1],edge_df_modified[,2],sep="-"),
                       share=edge_df_modified[,3])
  shareData2=aggregate(share ~ setPair, data = shareData, sum)
  sharesetAB=as.data.frame(str_split_fixed(shareData2$setPair,"-",n=2))
  g <- graph_from_data_frame(d=sharesetAB,directed = direction)#correction1
  lostnode=setdiff(rownames(lo),V(g)$name)
  if(length(lostnode)>0){g=g+vertices(lostnode)}
  if(identical(V(g)$name,names(nnt[["setID"]]))){
    print("Identical")
    setID=nnt[["setID"]]
  }else{
    setID=nnt[["setID"]][V(g)$name]
  }
  #setID=1:length(V(g)$name)
  #names(setID)=V(g)$name
  #niche2setID=lapply(niche,function(x){setID[intersect(names(x),names(setID))]})
  #newniche=lapply(niche,function(x){x[intersect(names(x),names(setID))]})
  # 设置顶点大小
  V(g)$size <- vsize*log10(1+set.size[V(g)$name])  # 通过节点名称匹配大小correction3
  #E(g)$size <- log2(1+shareData2$share)
  E(g)$size <- shareData2$share
  last_code<- str_split_fixed(V(g)$name, "_",n=2)[,2]
  set.seed(seed)
  color_pallete <- createPalette(nchar(last_code[1]),
                                 seedcolors = c("#E69F00" ,"#56B4E9" ,"#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7"))
  char_vectors <- strsplit(last_code, split = "")
  #print(char_vectors)
  set_df=sapply(char_vectors,as.integer)
  vcolors=apply(set_df,2,function(x){color_pallete[as.logical(x)]})
  #set_cluster=apply(set_df,2,function(x){nnt[["code2celltype"]][as.logical(x)]})
  # 绘制网络图
  lo=lo[V(g)$name,]
  V(g)$x=lo[,1]
  V(g)$y=lo[,2]
  write.csv(lo,file = file_layout,row.names=F)
  if(mark_groups){
    #library(Polychrome)
    first_code <-str_split_fixed(V(g)$name, "_",n=2)[,1]
    char_vectors <- strsplit(first_code, split = "")
    set_df=as.data.frame(sapply(char_vectors,as.integer))
    colnames(set_df)=1:ncol(set_df)
    set_1tag=set_df[,which(apply(set_df,2,sum)==1)]
    set_tag=apply(set_1tag,2,function(x){nnt[["code2tag"]][as.logical(x)]})
    set1tag_list=list()
    for(tag in nnt[["code2tag"]]){
      set1tag_list[[tag]]=as.integer(names(set_tag[set_tag==tag]))
    }
    mark.groups=set1tag_list[sapply(set1tag_list,length)>0]
    #library(scales)
    mark.col=hue_pal()(length(set1tag_list))
  }else{
    mark.groups = list()
    mark.col = rainbow(length(mark.groups), alpha = 0.3)
  }
  if(vertex2label==T){vertex.label=setID}else{vertex.label=NA}
  if(edge2label==T){edge.label=round(E(g)$size,2)}else{edge.label=NA}
  plot(g,
       axes = axes,
       mark.groups = mark.groups,
       mark.col=mark.col,
       #mark.border=mark.colors,
       mark.shape =1/2,
       vertex.frame.width = 0.5,
       vertex.label = vertex.label,
       vertex.label.cex = vertex.label.cex,
       edge.arrow.size = 0.1,
       edge.width = esize*log10(1+E(g)$size),  # 使用 size 作为边的宽度
       edge.label = edge.label,
       edge.label.cex = edge.label.cex,
       vertex.size = V(g)$size, # 使用 vector 作为顶点大小
       edge.color = ecolor,
       vertex.color = vcolors,
       margin = margin
  )

  legend("topleft", legend = celltypes,col = color_pallete,pch = 21, pt.bg = color_pallete, pt.cex = 1,cex = 0.5,bty = "n")

  if(mark_groups){legend("topright", legend = names(set1tag_list),col = mark.col,
                         pch = 21, pt.bg = mark.col, pt.cex = 1,cex = 0.5,bty = "n")}

  dev.off()
}

#' print tag vs cell number and types
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#'
#' @return pdf
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' print_nichetag(nnt)
print_nichetag<-function(nnt,file="nichetag.pdf"){
  #library(cowplot)
  niche2tagnum=nnt[["niche_tags"]]
  niche2cellnum=nnt[["niche_size"]]
  niche2typenum=nnt[["niche_celltypes"]]
  data0=data.frame(celltag=names(niche2tagnum),tagnumber=log10(niche2tagnum))
  data0$celltag=factor(data0$celltag,levels = data0$celltag)
  data1=data.frame(celltag=names(niche2cellnum),cellnumber=niche2cellnum)
  data1$celltag=factor(data1$celltag,levels = data1$celltag)
  data2=data.frame(celltag=names(niche2typenum),celltypes=niche2typenum)
  data2$celltag=factor(data2$celltag,levels = data2$celltag)
  p0 <- ggplot(data0, aes(x = celltag, y =tagnumber)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal(base_size = 7) + ylab("Log10(CellTag UMI)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p1 <- ggplot(data1, aes(x = celltag, y =cellnumber)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p2 <- ggplot(data2, aes(x = celltag, y =celltypes)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pg=plot_grid(p0,p1,p2,ncol=1,nrow = 3,hjust = "hv")
  ggsave(plot = pg,filename=file)
}


#' print cluster vs tag number and types
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#'
#' @return pdf
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' print_clustertag(nnt)
print_clustertag_v1<-function(nnt,file="clustertag.pdf"){
  tag_expression=nnt[["tag_matrix"]]
  cell_clusters=droplevels(nnt[["cell_type"]])
  plot.list=list()
  for (cluster_id in unique(cell_clusters)){
    print(cluster_id)
    cluster_cells <- names(cell_clusters)[cell_clusters == cluster_id]
    tag_expression_cluster <- tag_expression[cluster_cells, ]
    tag_counts_per_cell <- rowSums(tag_expression_cluster)
    unique_tag_counts  <- rowSums(tag_expression_cluster > 0)
    tag_counts_df <- data.frame(
      tag_count = tag_counts_per_cell,
      unique_tag_count = unique_tag_counts)
    tag_count_plot <- ggplot(tag_counts_df, aes(x = tag_count)) +
      geom_density(fill = "blue", alpha = 0.5) +
      labs(title = cluster_id,
           x = "Number of Tags per Cell",
           y = "Density") +
      theme_minimal()
    plot.list[[paste0(cluster_id,"_1")]]=tag_count_plot
    unique_tag_count_plot <- ggplot(tag_counts_df, aes(x = unique_tag_count)) +
      geom_density(fill = "red", alpha = 0.5) +
      labs(title = cluster_id,
           x = "Kinds of Tags per Cell",
           y = "Density") +
      theme_minimal()
    plot.list[[paste0(cluster_id,"_2")]]=unique_tag_count_plot
  }
  #library(cowplot)
  pg=plot_grid(plotlist = plot.list,ncol=2,align = "hv")
  ggsave(plot = pg,filename = file,height = 1.5*length(unique(cell_clusters)),limitsize = F)
}

#' print cluster vs tag number and types
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#'
#' @return pdf
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' print_clustertag(nnt)
print_clustertag<-function (nnt, file = "clustertag.pdf"){
  library(cowplot)
  tag_expression = nnt[["tag_matrix"]]
  cell_clusters = droplevels(nnt[["cell_type"]])
  plot.list = list()
  for (cluster_id in unique(cell_clusters)) {
    print(cluster_id)
    cluster_cells <- names(cell_clusters)[cell_clusters == cluster_id]
    tag_expression_cluster <- tag_expression[cluster_cells, ]
    tag_counts_per_cell <- rowSums(tag_expression_cluster)
    maxtag_counts_per_cell <- apply(tag_expression_cluster,1,max)
    othertag_counts_per_cell <- apply(tag_expression_cluster,1,function(x){sum(x)-max(x)})

    unique_tag_counts <- rowSums(tag_expression_cluster > 0)

    tag_counts_df <- data.frame(tag_count = tag_counts_per_cell,
                                unique_tag_count = unique_tag_counts,
                                max_tag_count = maxtag_counts_per_cell,
                                other_tag_count = othertag_counts_per_cell
                                #max_N2_FC = log1p(maxtag_counts_per_cell)-log1p(N2tag_counts_per_cell)
    )

    tag_count_plot <- ggplot(tag_counts_df, aes(x = tag_count)) +
      geom_density(fill = "blue", alpha = 0.5) +
      labs(title = cluster_id, x = "Number of Tags per Cell", y = "Density") + theme_minimal()
    plot.list[[paste0(cluster_id, "_1")]] = tag_count_plot
    unique_tag_count_plot <- ggplot(tag_counts_df, aes(x = unique_tag_count)) +
      geom_density(fill = "red", alpha = 0.5) +
      labs(title = cluster_id, x = "Kinds of Tags per Cell", y = "Density") + theme_minimal()
    plot.list[[paste0(cluster_id, "_2")]] = unique_tag_count_plot
    max_tag_count_plot <- ggplot(tag_counts_df, aes(x = max_tag_count)) +
      geom_density(fill = "green", alpha = 0.5) +
      labs(title = cluster_id, x = "Number of Primary Tags per Cell", y = "Density") + theme_minimal()
    plot.list[[paste0(cluster_id, "_3")]] = max_tag_count_plot
    other_tag_count_plot <- ggplot(tag_counts_df, aes(x = other_tag_count)) +
      geom_density(fill = "yellow", alpha = 0.5) +
      #geom_vline(xintercept = q95_density, color = "blue", linetype = "dashed", size = 0.5) +
      #geom_vline(xintercept = q99_density, color = "red", linetype = "dashed", size = 0.5) +
      labs(title = cluster_id, x = "Number of Non-primary Tags per Cell", y = "Density") +
      scale_x_continuous(limits = c(0,min(max(tag_counts_df$other_tag_count),20)))+
      theme_minimal()
    plot.list[[paste0(cluster_id, "_4")]] = other_tag_count_plot

  }
  pg = plot_grid(plotlist = plot.list, ncol = 4, align = "hv")
  ggsave(plot = pg, filename = file, height = 1.5 * length(unique(cell_clusters)), width=14,
         limitsize = F)
}


#' print cancer and non-maligant cell numbers for each tag
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#'
#' @return pdf
#' @importFrom stringr str_detect
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' tag_cancer_noncancer(nnt)
tag_cancer_noncancer<-function(nnt,file="tag2celltypes.pdf"){
  tag_expression=nnt[["tag_matrix"]]
  cell_clusters=droplevels(nnt[["cell_type"]])
  niche=apply(tag_expression,2,function(x){x[x>0]})
  #cancer_noncancer=ifelse(str_detect(tolower(cell_clusters),"cancer"),"cancer","noncancer")
  cancer=names(cell_clusters)[str_detect(tolower(cell_clusters),"cancer")]
  noncancer=names(cell_clusters)[!str_detect(tolower(cell_clusters),"cancer")]
  corn=function(x){
    a=sum(names(x) %in% cancer)
    b=sum(names(x) %in% noncancer)
    c=c(a,b)
    names(c)=c("cancer","noncancer")
    c
  }
  niche_celltype=as.data.frame(t(sapply(niche,corn)))
  niche_celltype$tag=rownames(niche_celltype)
  #library(reshape2)
  #library(ggplot2)
  data=melt(niche_celltype)
  data$tag=factor(data$tag,levels = niche_celltype$tag)
  gp=ggplot(data, aes(x = tag, y = value, fill =variable )) +
    geom_bar(stat = "identity") +theme_minimal(base_size = 7)+
    labs(x = "CellTag Sequence", y = "Cell Count", fill = "Category") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top")+
    scale_fill_manual(values = c("#0072B2","#009E73"))
  ggsave(plot=gp,filename=file,height = 3,width = 7)
}

#' print sender and receiver cell numbers for each tag, cutoff is 3
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#'
#' @return pdf
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' tag_cci(nnt)
tag_cci_v1<-function(nnt,file="tag2ccitypes.pdf"){
  tag_expression=nnt[["tag_matrix"]]
  cell_clusters=droplevels(nnt[["cell_type"]])
  niche=apply(tag_expression,2,function(x){x[x>0]})
  #cancer_noncancer=ifelse(str_detect(tolower(cell_clusters),"cancer"),"cancer","noncancer")
  sorr=function(x){
    a=sum(x>3)
    b=sum(x<=3)
    c=c(a,b)
    names(c)=c("sender","receiver")
    c
  }
  niche_celltype=as.data.frame(t(sapply(niche,sorr)))
  niche_celltype$tag=rownames(niche_celltype)
  #library(reshape2)
  #library(ggplot2)
  data=melt(niche_celltype)
  data$tag=factor(data$tag,levels = niche_celltype$tag)
  gp=ggplot(data, aes(x = tag, y = value, fill =variable )) +
    geom_bar(stat = "identity") +theme_minimal(base_size = 7)+
    labs(x = "CellTag Sequence", y = "Cell Count", fill = "Category") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top")+
    scale_fill_manual(values = c("#009E73","#E69F00"))
  ggsave(plot=gp,filename=file,height = 3,width = 7)
}

#' print sender and receiver cell numbers for each tag, cutoff for sender is 10 and 0.8
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param send_cutoff1    int, count cutoff for primary celltag asignment
#' @param send_cutoff2    numeric, percentage cutoff for primary celltag asignment
#' @param file pdf
#'
#' @return pdf
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' tag_cci(nnt)
tag_cci<-function(nnt, send_cutoff1=10, send_cutoff2=0.8, file = "tag2ccitypes.pdf")
{
  tag_expression = nnt[["tag_matrix"]]
  tag_expr=tag_expression[rowSums(tag_expression)>0,]
  cell2maxtag=apply(tag_expr,1,function(x){colnames(tag_expr)[which(x==max(x))[1]]})
  horl=apply(tag_expr, 1, function(x){(max(x)>=send_cutoff1)&(max(x)/sum(x)>=send_cutoff2)})
  horl2=horl | (!str_detect(nnt$cell_type[rownames(tag_expr)],"Cancer"))
  tag_expr2=tag_expr[horl2,]
  tag_expr3=tag_expr[!horl2,]
  niche = apply(tag_expr2, 2, function(x) {
    x[x > 0]
  })

  niche_celltype = as.data.frame(t(sapply(1:length(niche),function(x){
    niche_tag=niche[[x]]
    tag=names(niche)[x]
    non_cancer_cell=names(niche_tag)[!str_detect(nnt$cell_type[names(niche_tag)],"Cancer")]
    cancer_cell=names(niche_tag)[str_detect(nnt$cell_type[names(niche_tag)],"Cancer")]
    sender=cancer_cell[cell2maxtag[cancer_cell]==tag]
    receiver=cancer_cell[cell2maxtag[cancer_cell]!=tag]
    c=c(length(sender),length(receiver)+length(non_cancer_cell))
    names(c) = c("Celltag-sender cancer cells", "Celltag-receiver cells")
    c
  })))
  niche_celltype$tag = names(niche)
  level = niche_celltype$tag
  data = melt(niche_celltype)
  data$tag = factor(data$tag, levels=level)

  nonHtag2cells=apply(tag_expr3,2,function(x){sum(x>0)})
  lowdata=data.frame(tag=names(nonHtag2cells) ,variable="Celltag-indeterminate cancer cells", value=nonHtag2cells)
  #Celltag-sender cancer cells
  #Celltag-receiver cells
  #Celltag-indeterminate cancer cells

  data2=rbind(data,lowdata)

  gp1 = ggplot(data, aes(x = tag, y = value, fill = variable)) +
    geom_bar(stat = "identity") + theme_minimal(base_size = 7) +
    labs(x = "CellTag Sequence", y = "Cell Count", fill = "Category") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    scale_fill_manual(values = c("#009E73","#E69F00"))

  gp2 = ggplot(data2, aes(x = tag, y = value, fill = variable)) +
    geom_bar(stat = "identity") + theme_minimal(base_size = 7) +
    labs(x = "CellTag Sequence", y = "Cell Count", fill = "Category") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    scale_fill_manual(values = c("#009E73","#E69F00", "#CCCCCC"))

  pg=plot_grid(gp1,gp2,ncol=1,align = "hv")
  ggsave(plot = pg, filename = file, height = 7, width = 7)
}


#' print set type distribution for each niche
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#'
#' @return pdf
#' @import cowplot
#' @import ggplot2
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' settype(nnt)
settype_v1<-function(nnt,file="settype.pdf"){
  niche2set=nnt$niche
  niche2set.vec=unlist(lapply(names(niche2set),function(x){
    res=paste(x,names(niche2set[[x]]),sep=":")
    res
  }))
  niche2set.df=as.data.frame(stringr::str_split_fixed(niche2set.vec,":",n=2))
  colnames(niche2set.df) = c("niche","set_code")
  niche2set.df$set_ID=nnt$setID[niche2set.df$set_code]
  niche2set.df$celltype=codesplit(niche2set.df$set_code,nnt[["code2tag"]],nnt[["code2celltype"]])$types
  set2celltype=as.data.frame(sort(table(niche2set.df$celltype),decreasing=T))
  plot.list=list()
  for(tag in names(niche2set)){
    df=niche2set.df[niche2set.df$niche==tag,]
    data=as.data.frame(table(df$celltype))
    data2=data.frame(Var1=setdiff(set2celltype$Var1,data$Var1),Freq=0)
    data=rbind(data,data2)
    data$Var1=factor(data$Var1,levels=set2celltype$Var1)
    plot.list[[tag]]<- ggplot(data, aes(x = Var1, y =Freq)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      theme_minimal(base_size = 6) + ylab("Number")+ xlab(tag)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  pg=plot_grid(plotlist = plot.list,ncol=1,align = "hv")
  ggsave(plot=pg,filename=file,
         height =length(unique(niche2set.df$niche)),
         width = length(levels(set2celltype$Var1))/5,
         limitsize = FALSE)
}


#' print set type distribution for each niche
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#' @param onlyHighconfident TorF, if or not use remove sets of celltag indeterminate cancer cells
#'
#' @return pdf
#' @import cowplot
#' @import ggplot2
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' settype(nnt)
settype<-function(nnt, file = "settype.pdf",onlyHighconfident=TRUE){
  niche2set = nnt$niche
  niche2set.vec = unlist(lapply(names(niche2set), function(x) {
    res = paste(x, names(niche2set[[x]]), sep = ":")
    res
  }))
  niche2set.df = as.data.frame(stringr::str_split_fixed(niche2set.vec, ":", n = 2))
  colnames(niche2set.df) = c("niche", "set_code")
  niche2set.df$set_ID = nnt$setID[niche2set.df$set_code]
  niche2set.df$celltype = codesplit(niche2set.df$set_code,
                                      nnt[["code2tag"]], nnt[["code2celltype"]])$types

  if(onlyHighconfident==TRUE){
    niche2set.df=niche2set.df[str_detect(niche2set.df$set_code,"2")|(!str_detect(niche2set.df$celltype,"Cancer")),]
  }

  set2celltype = as.data.frame(sort(table(niche2set.df$celltype),
                                      decreasing = T))
  plot.list = list()
  plot.list[["AllFMU"]]=ggplot(set2celltype, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", fill = "skyblue") + theme_minimal(base_size = 6) +
    ylab("Set Count") + xlab("AllFCU") + theme(axis.text.x = element_text(angle = 45,
                                                                           hjust = 1))
  for (tag in names(niche2set)) {
    df = niche2set.df[niche2set.df$niche == tag, ]
    data = as.data.frame(table(df$celltype))
    data2 = data.frame(Var1 = setdiff(set2celltype$Var1,
                                      data$Var1), Freq = 0)
    data = rbind(data, data2)
    data$Var1 = factor(data$Var1, levels = set2celltype$Var1)
    plot.list[[tag]] <- ggplot(data, aes(x = Var1, y = Freq)) +
      geom_bar(stat = "identity", fill = "skyblue") + theme_minimal(base_size = 6) +
      ylab("Set Count") + xlab(tag) + theme(axis.text.x = element_text(angle = 45,
                                                                        hjust = 1))
  }
  pg = plot_grid(plotlist = plot.list, ncol = 1, align = "hv")
  ggsave(plot = pg, filename = file, height = (length(unique(niche2set.df$niche))+1),
         width = length(levels(set2celltype$Var1))/5, limitsize = FALSE)
}


#' print cell type and set type distribution for each niche
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#' @param seed seed for color palette creation
#'
#' @return pdf
#' @import cowplot
#' @import ggplot2
#' @import reshape2
#' @export
#'
#' @examples
#' data(example)
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' tag_cellsettype(nnt)
tag_cellsettype<-function(nnt,file="tag2celltype_settype.pdf",seed=888){
  tag_expression=nnt[["tag_matrix"]]
  #tag_expression=tag_expression[apply(tag_expression,1,sum)>0,]
  cell_clusters=droplevels(nnt[["cell_type"]])
  set.seed(seed)
  color_pallete <- createPalette(length(levels(cell_clusters)),
                                 seedcolors = c("#E69F00" ,"#56B4E9" ,"#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7"))
  names(color_pallete)=NULL
  niche=apply(tag_expression,2,function(x){x[x>0]})
  niche_celltype=as.data.frame(t(sapply(niche,function(x){table(cell_clusters[names(x)])})))
  niche_celltype$tag=rownames(niche_celltype)
  #library(reshape2)
  #library(ggplot2)
  #library(cowplot)
  data=melt(niche_celltype)
  data$tag=factor(data$tag,levels = niche_celltype$tag)
  gp1=ggplot(data, aes(x = tag, y = value, fill =variable )) +
    geom_bar(stat = "identity") +theme_minimal(base_size = 7)+
    labs(x = "CellTag Sequence", y = "Cell Count", fill = "Category") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top")+
    scale_fill_manual(values = color_pallete)+
    theme(
      legend.text = element_text(size = 6),        # 图例文字大小
      legend.title = element_text(size = 7),       # 图例标题大小（可选）
      legend.key.size = unit(0.2, "cm")            # 图例图形大小
    )

  sets=names(nnt$setID)
  set_types=codesplit(sets,nnt[["code2tag"]],nnt[["code2celltype"]])$types
  set_types=factor(set_types,levels =nnt[["code2celltype"]])
  names(set_types)=sets

  set.seed(seed)
  color_pallete <- createPalette(length(levels(set_types)),
                                 seedcolors = c("#E69F00" ,"#56B4E9" ,"#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7"))
  names(color_pallete)=NULL
  niche=nnt$niche
  niche_settype=as.data.frame(t(sapply(niche,function(x){table(set_types[names(x)])})))
  niche_settype$tag=rownames(niche_settype)
  data=melt(niche_settype)
  data$tag=factor(data$tag,levels = niche_settype$tag)
  gp2=ggplot(data, aes(x = tag, y = value, fill =variable )) +
    geom_bar(stat = "identity") +theme_minimal(base_size = 7)+
    labs(x = "CellTag Sequence", y = "Set Count", fill = "Category") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top")+
    scale_fill_manual(values = color_pallete)+
    theme(
      legend.text = element_text(size = 6),        # 图例文字大小
      legend.title = element_text(size = 7),       # 图例标题大小（可选）
      legend.key.size = unit(0.2, "cm")            # 图例图形大小
    )
  pg=plot_grid(gp1,gp2,ncol=1,align = "hv")
  ggsave(plot=pg,filename=file,height = 7,width = 7)
}
