.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
#writer:slm
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(DOSE))
suppressMessages(library(clusterProfiler))

#function:according to the expression of ENS and ENS-SYS,merge the symbol and expression
#input:ENS-expression ; ENS-sys
#output:SYS-expression
sys_expression <- function(ens_sys,ens_expression){
  colnames(ens_sys) <- c("ENS","SYS")
  colnames(ens_expression)[1] <- "ENS"
  rownames(ens_expression) <- ens_expression[,1]
  ens_expression <- ens_expression[,-1]
  ens_expression <- ens_expression[which(rowSums(ens_expression) > 0),]
  ens_expression$ENS <- rownames(ens_expression)
  all_merge <- merge(ens_sys,ens_expression,by="ENS")
  all_merge <- all_merge[,-1]
  all_mer <- as.data.frame(lapply(all_merge, as.numeric))
  all_mer[,1] <- all_merge[,1]
  n=ncol(all_mer)
  all_la <- aggregate(all_mer[2:n],by=all_mer[1],FUN=function(X)sum(X))
  return(all_la)
}

#genelist : just one column,no header,ENSEMBL,ENTREZID or symbol
#type and species are all just one word
#function : genlist(one column)  --->  new_genelist(ENSEMBL,ENTREZID,SYMBOL)
genelist_to_allgene <- function(genelist,type,species){
  if (species=="mouse") {
    org=org.Mm.eg.db
  }else if (species=="human") {
    org=org.Hs.eg.db
  }
  if(type=="ENSEMBL"){
    colnames(genelist) <- "ENSEMBL"
    genelist$ENSEMBL <- as.character(genelist$ENSEMBL)
    genelist$symbol <- mapIds(x=org,keys = genelist$ENSEMBL,keytype = "ENSEMBL",column ="SYMBOL" )
    genelist <- na.omit(genelist)
    genelist$ENTREZID <- mapIds(x=org,keys = genelist$ENSEMBL,keytype = "ENSEMBL",column ="ENTREZID" )
    genelist <- na.omit(genelist)
  }else if (type=="SYMBOL") {
    colnames(genelist) <- "symbol"
    genelist$symbol <- as.character(genelist$symbol)
    genelist$ENSEMBL <- mapIds(x=org,keys = genelist$symbol,keytype = "SYMBOL",column ="ENSEMBL" )
    genelist <- na.omit(genelist)
    genelist$ENTREZID <- mapIds(x=org,keys = genelist$symbol,keytype = "SYMBOL",column ="ENTREZID" )
    genelist <- na.omit(genelist)
  }else if (type=="ENTREZID") {
    colnames(genelist) <- "ENTREZID"
    genelist$ENTREZID <- as.character(genelist$ENTREZID)
    genelist$ENSEMBL <- mapIds(x=org,keys = genelist$ENTREZID,keytype = "ENTREZID",column ="ENSEMBL" )
    genelist <- na.omit(genelist)
    genelist$symbol <- mapIds(x=org,keys = genelist$ENTREZID,keytype = "ENTREZID",column ="SYMBOL" )
    genelist <- na.omit(genelist)
  }
  genelist$symbol <- make.unique(genelist$symbol)
  new_genelist <- data.frame(genelist$ENSEMBL,genelist$ENTREZID,genelist$symbol)
  return(new_genelist)
}

#function : extract the result of function-enrichment from genelist
#input : the list (gene(ENTREZID)) , the character (treatment[order])
#output : the list (gofunction top10 and ad.pvalue)
top_gofun_pva <- function(genelist,metalist,type,species,top){
  if (species=="mouse") {
    org=org.Mm.eg.db
  }else if (species=="human") {
    org=org.Hs.eg.db
  }
  enrich <-list()
  func <- c()
  for(i in 1:(length(metalist))){
    test <- as.data.frame(genelist[i])
    colnames(test) <- "gene"
    test$gene <- as.character(test$gene)
    go <- enrichGO(gene = test$gene,ont=type,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05,OrgDb = org.Hs.eg.db,readable = TRUE)
    test_go <- data.frame(go@result$Description)
    topn <- as.data.frame(test_go[1:top,1])
    colnames(topn) <- "enrichment"
    func <- c(func,as.character(topn$enrichment))
    enrich[i] <- go
  }
  func <- unique(func)
  func <- as.data.frame(func)
  for (i in 1:length(metalist)) {
    tmp_p <- data.frame(enrich[[i]]@result$Description,enrich[[i]]@result$p.adjust)
    colnames(tmp_p) <- c("func","adjustpvalue")
    if (i==1) {
      matri <- merge(func,tmp_p,by="func")
    }else if (i!=1) {
      matri <- merge(matri,tmp_p,by="func")
    }
    colnames(matri)[i+1] <- metalist[i]
  }
  matri <- reshape2::melt(matri,id.var="func")
  colnames(matri) <- c("enrichment","treatment","adjust.value")
  return(matri)
}

