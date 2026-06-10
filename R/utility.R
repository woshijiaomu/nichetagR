#' Title
#'
#' @param tag_expression    dataframe, input tag expression matrix in which cells are rownames and tag expression are colnames
#'
#' @return    dataframe, removed uncorrected tag
#' @importFrom stringdist stringdist
#' @noRd
quality_control_v1<-function(tag_expression){
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
  #library(stringdist)
  cbn=combn(colnames(tag_expression),2)
  similar=cbn[,apply(cbn,2,function(x){stringdist(x[1], x[2], method = "lv")})<=2]
  tag_sum=apply(tag_expression,2,sum)
  print(tag_sum[similar[1,]])
  print(tag_sum[similar[2,]])
  tag_expression=tag_expression[,!(colnames(tag_expression) %in% similar[2,])]
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
  print(apply(tag_expression,2,sum))
  tag_expression
}

#' Remove closely related CellTags and CellTags with no resource cells
#'
#' @param tag_expression    dataframe, input tag expression matrix in which cells are rownames and tag expression are colnames
#' @param highconfidence    TorF, if or not to distinguish high confident tag-sender
#' @param send_cutoff1    int, count cutoff for primary celltag asignment
#' @param send_cutoff2    numeric, percentage cutoff for primary celltag asignment
#' @return    dataframe, removed uncorrected tag
#' @importFrom stringdist stringdist
#' @noRd
quality_control<-function(tag_expression,highconfidence=TRUE,send_cutoff1=10,send_cutoff2=0.8){
  print(dim(tag_expression))
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
  #library(stringdist)
  cbn=combn(colnames(tag_expression),2)
  similar=cbn[,apply(cbn,2,function(x){stringdist(x[1], x[2], method = "lv")})<=2]
  #tag_sum=apply(tag_expression,2,sum)
  #print(tag_sum[similar[1,]])
  #print(tag_sum[similar[2,]])
  tag_expression=tag_expression[,!(colnames(tag_expression) %in% similar[2,])]
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
  hpos=c()
  if(highconfidence==TRUE){
    for(tagposition in 1:ncol(tag_expression)){
      s=apply(tag_expression,1,function(x){
        x[tagposition]>=send_cutoff1 & x[tagposition]/sum(x)>=send_cutoff2
      })
      if(sum(s)>0){hpos=c(hpos,tagposition)}
    }
    #print(hpos)
    tag_expression=tag_expression[,hpos]
  }
  tag_expression
}


#' Translate code to tag and cell type
#'
#' @param code    code of cells
#' @param code2tag    each position of code corresponds to a tag
#' @param code2type    each position of code corresponds to a celltype
#'
#' @return     text of tag and celltype
#' @importFrom stringr str_split_fixed
#' @noRd
codesplit_v1<-function(code,code2tag,code2type){
  res=list()
  #library(stringr)
  codes=str_split_fixed(code,"_",n=2)
  codes_tag=codes[,1]
  codes_type=codes[,2]
  char_vectors <- strsplit(codes_tag, split = "")
  set_df=sapply(char_vectors,as.integer)
  tags=apply(set_df,2,function(x){code2tag[as.logical(x)]})
  res[["tags"]]=tags
  char_vectors <- strsplit(codes_type, split = "")
  set_df=sapply(char_vectors,as.integer)
  types=apply(set_df,2,function(x){code2type[as.logical(x)]})
  res[["types"]]=types
  res
}

#' Translate code to tag and cell type
#'
#' @param code    code of cells
#' @param code2tag    each position of code corresponds to a tag
#' @param code2type    each position of code corresponds to a celltype
#'
#' @return    tag, lineage tag and celltype
#' @importFrom stringr str_split_fixed
#' @noRd
codesplit<-function(code, code2tag, code2type){
  res = list()
  codes = str_split_fixed(code, "_", n = 2)
  codes_tag = codes[, 1]
  codes_type = codes[, 2]
  char_vectors <- strsplit(codes_tag, split = "")
  set_df = sapply(char_vectors, as.integer)
  tags = apply(set_df, 2, function(x) {
    code2tag[as.logical(x)]
  })
  res[["tags"]] = tags
  lineage=apply(set_df, 2, function(x) {
    code2tag[x==2]
  })
  res[["lineage"]] = lineage
  char_vectors <- strsplit(codes_type, split = "")
  set_df = sapply(char_vectors, as.integer)
  types = apply(set_df, 2, function(x) {
    code2type[as.logical(x)]
  })
  res[["types"]] = types
  res
}
