.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
#writer:slm
suppressMessages(library(DESeq2))
suppressMessages(library(DOSE))
suppressMessages(library(igraph))
suppressMessages(library(reshape2))
suppressMessages(library(enrichplot))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(getopt))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggraph))
suppressMessages(library(ggridges))
suppressMessages(library(ggpubr))
suppressMessages(library(limma))
suppressMessages(library(clusterProfiler))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(forcats))
suppressMessages(library(ggupset))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(plotrix))
suppressMessages(library(topGO))
suppressMessages(library(pathview))
suppressMessages(library(ggstream))
suppressMessages(library(RColorBrewer))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))

#input:data(type:matrix | col:description,variable[treatment],value) | meta(type:list)
#variable : include three character,treat1[+],treat2[-],cover[+ and -]
#meta :  treat1 treat2 cover 
dumbbell_chart_cover <- function(data,meta){
  colnames(data) <- c("description","variable","value")
  tomax <- max(abs(data$value))
  max=ceiling(tomax/10)*10
  treat1 <- subset(data,variable==meta[1])
  treat2 <- subset(data,variable==meta[2])
  cover <- subset(data,variable==meta[3])
  cover_n <- nrow(cover)
  cover_n_ban <- cover_n/2
  cover1 <- cover[1:cover_n_ban,]
  cover2 <- cover[(cover_n_ban+1):cover_n,]
  ggplot(data)+
    geom_point(aes(x=value,y = description,color=variable),
               size=5,show.legend = TRUE)+
    scale_x_continuous(limits = c(-max,max),breaks = seq(-max,max,10))+
    geom_segment(data=treat1,aes(x = value,y = description,
                                yend=treat2$description,
                                xend=treat2$value),
                 color="#aeb6bf",size=4.5,alpha=.5)+
    geom_segment(data=cover1,aes(x=value,y=description,
                               xend=cover2$value,
                               yend=cover2$description),
                 color="green",size=4.5,alpha=.5)+
    labs(x="Gene number",y="function")
}

#input : genelist (one column : ENTREZID ; no header)
#output : the list include dotplot(all,cc,bp,mf),barplot(all,cc,bp,mf)
GO_dotbar <- function(genelist,species){
  colnames(genelist) <- "gene"
  genelist$gene <- as.character(genelist$gene)
  if (species=="mouse") {
    org=org.Mm.eg.db
  }else if (species=="human") {
    org=org.Hs.eg.db
  }
  ALL <- enrichGO(gene = genelist$gene,OrgDb = org,ont="ALL",pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  CC <- enrichGO(gene = genelist$gene,OrgDb = org,ont="CC",pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  BP <- enrichGO(gene = genelist$gene,OrgDb = org,ont="BP",pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  MF <- enrichGO(gene = genelist$gene,OrgDb = org,ont="MF",pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  aldot <- dotplot(ALL,showCategory=20,title="Enrichment_GO_ALL_dot")
  ccdot <- dotplot(CC,showCategory=20,title="Enrichment_GO_CC_dot")
  bpdot <- dotplot(BP,showCategory=20,title="Enrichment_GO_BP_dot")
  mfdot <- dotplot(MF,showCategory=20,title="Enrichment_GO_MF_dot")
  albar <- barplot(ALL,showCategory=20,title="Enrichment_GO_ALL_bar")
  ccbar <- barplot(CC,showCategory=20,title="Enrichment_GO_CC_bar")
  bpbar <- barplot(BP,showCategory=20,title="Enrichment_GO_BP_bar")
  mfbar <- barplot(MF,showCategory=20,title="Enrichment_GO_MF_bar")
  return(list(aldot,ccdot,bpdot,mfdot,albar,ccbar,bpbar,mfbar))
}

#input : genelist (one column : ENTREZID ; no header)
#output : dot , bar and pathway
KEGG_path_one <- function(genelist,species){
  colnames(genelist) <- "gene"
  genelist$gene <- as.character(genelist$gene)
  if(species=="mouse"){
    org="mmu"
  }else if (species=="human") {
    org="hsa"
  }
  kk <- enrichKEGG(gene = genelist$gene,organism = org,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  dot <- dotplot(kk,title="Enrichment_KEGG_pathway_dot",showCategory=10)
  bar <- barplot(kk,title="Enrichment_KEGG_pathway_bar",showCategory=10)
  kk <- enrichKEGG(gene=genelist$gene,organism=org,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  path <- kk@result$ID[1:5]
  result <- c()
  for(pa in path){
    print(pa)
    pp <- pathview(gene.data=genelist$gene,pathway.id = pa,species = org)
    result <- c(result,pp)
  }
  return(list(dot,bar,result))
}

#input : genelist(include three columns : treatment,gene,fpkm)
#ouput : river plot([1]=river,[2]=%)
river_plot <- function(genelist){
  colnames(genelist) <- c("time","gene","value")
  genelist$gene <- factor(genelist$gene,levels = unique(genelist$gene))
  p <- ggplot(genelist,aes(time,value,
                           fill=gene,
                           color=gene,
                           label=gene))
  pri <- p+geom_stream(extra_span = 1.1,bw=0.75,type = c("ridge"))
  pbai <- p+geom_stream(extra_span = 1.1,bw=0.75,type = c("proportional"))
  return(list(pri,pbai))
}

#function : to show the heatmap
#input : matrix (include three columns : Y, X, value),as the reshap2
#output : heatmap 
heatmap <- function(list){
  colnames(list) <- c("varible","condition","value")
  ggplot(data = list,aes(x = condition,y = varible,fill=value))+
    geom_tile()+
    theme_minimal()+
    scale_y_discrete(position = "right")+
    scale_fill_continuous(low="white",high="red")+
    xlab(NULL) + ylab(NULL)+
    theme_minimal()
}



