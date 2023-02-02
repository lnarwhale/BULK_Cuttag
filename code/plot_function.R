.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
suppressMessages(library(getopt))
suppressMessages(library(DESeq2))
suppressMessages(library(DOSE))
suppressMessages(library(igraph))
suppressMessages(library(reshape2))
suppressMessages(library(ggraph))
suppressMessages(library(enrichplot))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(forcats))
suppressMessages(library(ggupset))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Ss.eg.db))
suppressMessages(library(org.Rn.eg.db))
suppressMessages(library(plotrix))

  
bases_content_lineplot <- function(stat_melt,len){
  ggplot(stat_melt, aes(x=position, y=percent, group=type)) +
  geom_line(aes(color=type),size=0.8) +
  theme_bw()+theme(panel.grid = element_blank()) +
  xlab("Position along reads") + ylab("Percent of bases (%)") +
  geom_vline(xintercept=len,lty=2,col="black",lwd=1)
}

read_dis_lineplot <- function(quali_data){
  ggplot(quali_data, aes(x=position, y=log10_coverage,fill=type)) +
  geom_line(size=0.8,aes(color=type)) +
  theme_bw()+theme(panel.grid = element_blank()) +
  xlab("Gene body percentile") + ylab("log10_coverage") 
}

error_barplot <- function(error_meanqual,len){
  ggplot(error_meanqual,aes(x=column,y=error))+
  geom_bar(stat='identity',width = 0.5,color='grey')+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Position along reads") + ylab("Error")+
  geom_vline(xintercept=len,lty=2,col="black",lwd=1)
}

pieplot <- function(filter_df){
  ggplot(data = filter_df, mapping = aes(x = 'Content', y = count, fill = perc)) + 
  geom_bar(stat = 'identity', position = 'stack',width=1)+
  coord_polar(theta = 'y')+
  labs(x = '', y = '', title = '')+
  theme_bw()+theme(axis.text = element_blank(),panel.grid = element_blank())+
  theme(axis.ticks = element_blank())
}

rseqc_barplot <- function(rseqc_stat){
  ggplot(rseqc_stat,aes(x=Group,y=Tags.Kb))+
  geom_bar(stat='identity',width = 0.5,color='grey')+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Group") + ylab("Tags/Kb")
}

fpkm_dis_boxplot <- function(fpkm_melt){
  ggplot(fpkm_melt, aes(x=sample, y=fpkm,fill=sample)) + 
  geom_boxplot()+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Sample") + ylab("FPKM(log10)")+
  scale_x_discrete(breaks=NULL)+
  scale_fill_discrete(name="Sample")
}
fpkm_dis_boxplot_rename <- function(fpkm_melt_tmp){
  ggplot(fpkm_melt_tmp, aes(x = variable, y=value,fill=variable)) + 
  geom_boxplot()+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Sample") + ylab("FPKM(log10)")+
  scale_fill_discrete(name="Sample")+
  theme(axis.text.x=element_text(angle=90, hjust=1))
}

cluster_plot <- function(dds){
    rld <- rlog(dds,blind=FALSE)
    sampleDist <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDist)
    rownames(sampleDistMatrix) <- rld$condition
    colnames(sampleDistMatrix) <- rld$condition
    cluster_plot <- pheatmap(sampleDistMatrix,display_numbers = FALSE,number_color = "black",clustering_distance_rows=sampleDist,clustering_distance_cols=sampleDist,color = colorRampPalette(c("blue","white"))(256))
    #fontsize=40,treeheight_row=200,treeheight_col=200,
}

scatter_plot <- function(fpkm_matrix){
  pairs(gene_fpkm_matrix,upper.panel = panel.cor,lower.panel = panel.lm,diag.panel = panel.hist)
}

PCA_plot <- function(dds){
    rld <- rlog(dds,blind=FALSE)
    PCA_data <- plotPCA(rld,returnData=TRUE)
    percentVar <- round(100*attr(PCA_data,"percentVar"))
    ggplot(PCA_data,aes(PC1,PC2,color=condition))+
    geom_point(size=3)+
    xlab(paste0("PC1:",percentVar[1],"% variance"))+ylab(paste0("PC2:",percentVar[2],"% variance"))+
    theme_bw()+theme(panel.grid = element_blank())+geom_text_repel(aes(label=PCA_data$name),size = 3)
}

MA_plot <- function(res){
    ggmaplot(res,
	fdr=0.05,fc=2,size=0.3,
	palette=c("red","blue","grey"),
	genenames=as.vector(rownames(res)),
	xlab="mean of normalized counts",
	ylab="log2FoldChange",
	legend="top",top=25,
	font.label=c("plain",7,"black"),
	font.legend=c("plain",7,"black"),
	ggtheme=ggplot2::theme_classic())
}

MAtrans_plot <- function(res){
    ggmaplot(res,
	fdr=0.05,fc=2,size=0.3,
	palette=c("red","blue","grey"),
	genenames=as.vector(res$Ensembl_Transcript_symbol),
	#genenames=as.vector(rownames(res)),
	xlab="mean of normalized counts",
	ylab="log2FoldChange",
	legend="top",top=25,
	font.label=c("plain",7,"black"),
	font.legend=c("plain",7,"black"),
	ggtheme=ggplot2::theme_classic())
}

Volcano_plot <- function(res){
    res$threshold <- factor(ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 1, ifelse(res$log2FoldChange>= 1 ,'Up','Down'),'NS'),levels=c('Up','Down','NS'))
    res_data = res[res$padj<0.05&abs(res$log2FoldChange)>1,]
	print("03")
	print(head(res_data))
	ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
	  geom_point(size=1)+
	  scale_color_manual(values=c("Up"="#DC143C","Down"="#00008B","NS"="#808080"))+
	  geom_text_repel(
	    data = res_data[1:25,],
	    aes(label = Gene),
        size = 3,
        segment.color = "black", show.legend = FALSE )+
	  theme_bw()+theme(panel.grid = element_blank())+
	  theme(
	    legend.title = element_blank()
	  )+
	  ylab('-log10 (p-adj)')+
	  xlab('log2 (FoldChange)')+
	  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5)+
	  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
	#print("04")
}

diff_Heatmap_plot <- function(mat,metadata){
    pheatmap(mat,annotation_col=metadata,fontsize=10,treeheight_row=100,treeheight_col=100)
}
diff_Heatmap_plot_uncluster <- function(mat,metadata){
    pheatmap(mat,annotation_col=metadata,fontsize=10,cluster_rows=FALSE)
}

GO_dot_plot <- function(go_data,go_number){
  
  if(go_number < 20){
    df <- fortify(go_data,showCategory=go_number,by=x)
	labels=(sapply(levels(df$Description)[as.numeric(df$Description)],shorten_names))
    names(labels) = rev(1:nrow(df))
    p <- ggplot(df,aes(x = Count,y = Description))+
         geom_point(aes(color = p.adjust,size=GeneRatio))+
         theme_bw()+
         scale_x_discrete(labels=labels)+
		 theme(panel.grid = element_blank())+
		 theme_dose(12)+
         scale_color_gradient(low = "red", high = "blue")+
		 scale_fill_continuous(low="red", high="blue", name = "p.adjust", guide=guide_colorbar(reverse=TRUE))+
         xlab("Count")
    #dotplot(go_data,showCategory=go_number)
  }else{
    df <- fortify(go_data,showCategory=20,by=x)
	labels=(sapply(levels(df$Description)[as.numeric(df$Description)],shorten_names))
    names(labels) = rev(1:nrow(df))
    p <- ggplot(df,aes(x = Count,y = Description))+
         geom_point(aes(color = p.adjust,size=GeneRatio))+
         theme_bw()+
         scale_x_discrete(labels=labels)+
		 theme(panel.grid = element_blank())+
		 theme_dose(12)+
         scale_color_gradient(low = "red", high = "blue")+
		 scale_fill_continuous(low="red", high="blue", name = "p.adjust", guide=guide_colorbar(reverse=TRUE))+
         xlab("Count")
    #dotplot(go_data,showCategory=go_number)
  }
  
}

GO_bar_plot <- function(go_data,go_number){
  if(go_number < 20){
    go_data <- fortify(go_data,showCategory=go_number,by=x)
    go_data$number <- factor(rev(1:nrow(go_data)))
	shorten_names <- function(x, n_word=4, n_char=200){
      if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 200))
      {
        if (nchar(x) > 200) x <- substr(x, 1, 200)
        x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],collapse=" "), "...", sep="")
        return(x)
      } 
      else
      {
        return(x)
      }
    }
	labels=(sapply(levels(go_data$Description)[as.numeric(go_data$Description)],shorten_names))
    df <- fortify(go_data,showCategory=go_number,by=x)
	labels=(sapply(levels(df$Description)[as.numeric(df$Description)],shorten_names))
    names(labels) = rev(1:nrow(go_data))
	print(labels)
	p <- ggplot(go_data, aes_string(x = go_data$number,y="Count",fill="p.adjust"))+
	     geom_bar(stat = "identity")+coord_flip()+
		 theme_bw()+
		 scale_x_discrete(labels=labels)+
	     theme(panel.grid = element_blank())+
		 theme_dose(12)+
		 scale_fill_continuous(low="red", high="blue", name = "p.adjust", guide=guide_colorbar(reverse=TRUE))+
         xlab("Description") + ylab("GeneCount")
	#p <- ggplot(data=df, aes_string(x = "Description",y="Count",fill="p.adjust")) +
         #geom_bar(stat="identity", width=0.8) + coord_flip() +
         #scale_fill_manual(values = CPCOLS) + theme_bw() +
         #scale_x_discrete(labels=labels) +
         #xlab("GO term") +
         #theme(axis.text=element_text(face = "bold", color="gray50")) +
         #labs(title = "The Most Enriched GO Terms")
  }else{
    go_data <- fortify(go_data,showCategory=20,by=x)
    go_data$number <- factor(rev(1:nrow(go_data)))
	shorten_names <- function(x, n_word=4, n_char=200){
      if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 200))
      {
        if (nchar(x) > 200) x <- substr(x, 1, 200)
        x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],collapse=" "), "...", sep="")
        return(x)
      } 
      else
      {
        return(x)
      }
    }
	labels=(sapply(levels(go_data$Description)[as.numeric(go_data$Description)],shorten_names))
    df <- fortify(go_data,showCategory=go_number,by=x)
	labels=(sapply(levels(df$Description)[as.numeric(df$Description)],shorten_names))
    names(labels) = rev(1:nrow(go_data))
	print(labels)
	p <- ggplot(go_data, aes_string(x = go_data$number,y="Count",fill="p.adjust"))+
	     geom_bar(stat = "identity")+coord_flip()+
		 theme_bw()+
		 scale_x_discrete(labels=labels)+
	     theme(panel.grid = element_blank())+
		 theme_dose(12)+
		 scale_fill_continuous(low="red", high="blue", name = "p.adjust", guide=guide_colorbar(reverse=TRUE))+
         xlab("Description") + ylab("GeneCount")
	#p <- ggplot(data=df, aes_string(x = "Description",y="Count",fill="p.adjust")) +
         #geom_bar(stat="identity", width=0.8) + coord_flip() +
         #scale_fill_manual(values = CPCOLS) + theme_bw() +
         #scale_x_discrete(labels=labels) +
         #xlab("GO term") +
         #theme(axis.text=element_text(face = "bold", color="gray50")) +
         #labs(title = "The Most Enriched GO Terms")
  }
}



GO_combine_bar_plot <- function(go_enrich_df){
    go_enrich_df <- fortify(go_enrich_df)
    go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
	labels=go_enrich_df$Description
    names(labels) = rev(1:nrow(go_enrich_df))
    CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
    ggplot(data=go_enrich_df, aes(x=number, y=-log(qvalue),fill=type)) +
    geom_bar(stat="identity", width=0.8) + coord_flip() + 
	geom_text(aes(label=GeneNumber),hjust=-0.1)+
    scale_color_gradient(low = "blue", high = "red")+
    scale_fill_manual(values = CPCOLS) + 
    theme_bw() + 
    scale_x_discrete(labels=labels) +
    xlab("GO term") + 
    ylab("-log10(qvalue)") + 
    theme_bw()+theme(panel.grid = element_blank())
}

GO_combine_dotplot <- function(go_enrich_df){
    go_enrich_df <- fortify(go_enrich_df)
    go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
	labels=go_enrich_df$Description
    names(labels) = rev(1:nrow(go_enrich_df))
	go_enrich_df$GeneRatio <- parse_ratio(go_enrich_df$GeneRatio)
	CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
    #y=factor(goinput$GO_Name,levels = goinput$GO_Name
    ggplot(data=go_enrich_df, aes(x=number, y=GeneRatio)) +
    geom_point(aes(size=GeneNumber,color=-log(qvalue),shape=type,))+
    coord_flip()+
    scale_color_gradient(low = "blue", high = "red")+
    #scale_fill_manual(values = CPCOLS) + 
    scale_x_discrete(labels=labels) +
    xlab("GO term") +
	theme_bw()
}

GO_net_plot <- function(go_data){
    plotGOgraph(go_data)
}

GO_enMap_plot <- function(go_data){
    emapplot(go_data, pie_scale=1.5,layout="kk")
	
}

GO_cent_plot <- function(go_data,genelist){
    cnetplot(go_data,categorySize="pvalue", foldChange=genelist,node_label="all")
}
GO_cent_plot_circle <- function(go_data,genelist){
    cnetplot(go_data,foldChange=genelist,circular = TRUE, colorEdge = TRUE,,categorySize="pvalue",node_label="all")
}
GO_heatmap_plot <- function(go_data,genelist){
    heatplot(go_data,foldChange=genelist)
}
SNP_dis_plot <- function(SNP_data){
  ggplot(SNP_data, aes(x = sample,y = value,fill = type))+
  geom_bar(stat ="identity")+
  labs(x = "Sample",y = "Counts")+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=90, hjust=1))
}
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                       collapse=" "), "...", sep="")

    return(x)
  }
  else
  {
    return(x)
  }
}

panel.cor<-function(x,y){
  usr <- par('usr')
  on.exit(par(usr))
  par(usr=c(0,1,0,1))
  r<-cor(x,y)
  text<-format(r,digits = 2)
  text(0.5,0.5,text,cex=abs(r)*5)
}

panel.hist<-function(x,y){
  usr <- par('usr')
  on.exit(par(usr))
  par(usr = c(usr[1:2],0, 1.5) )
  h<-hist(x,plot=F)
  breaks<-h$breaks
  b_length<-length(breaks)
  y<-h$counts/max(h$counts)
  rect(breaks[-b_length],0,breaks[-1],y,col='blue')
}

panel.lm<-function(x,y){
  points(x,y,col='black',cex=1)
  abline(lm(y~x),col='red')
}


#modified from clusterProfiler

GSEA_summary <- function(objectall,object,gsedata,geneSetID){
  gsea_data_frame <- as.data.frame(object)
  enrichment_name <- as.character(gsea_data_frame$ID[geneSetID])
  #print(enrichment_name)
  msigdf <- objectall
  msigdf_geneSet <- msigdf[which(msigdf$geneset==enrichment_name),]
  colnames(msigdf_geneSet) <- c("geneset","symbol")
  gsea_core_gene <- unlist(strsplit(gsea_data_frame$core_enrichment[geneSetID],split="/"))
  gsea_core_gene <- data.frame(gene=gsea_core_gene)
  gsea_core_gene$core_enrich_gene <- "YES"
  #print(head(msigdf_geneSet))
  
  msigdf_geneSet$ENSEMBL <- mapIds(spi_database,keys=as.character(msigdf_geneSet$symbol),column="ENSEMBL",keytype="SYMBOL",multiVals="first")
  msigdf_geneSet$ENTREZID <- mapIds(spi_database,keys=as.character(msigdf_geneSet$symbol),column="ENTREZID",keytype="SYMBOL",multiVals="first")
  msigdf_geneSet$GENENAME <- mapIds(spi_database,keys=as.character(msigdf_geneSet$symbol),column="GENENAME",keytype="SYMBOL",multiVals="first")
  #colnames(gsedata) <- c("ENTREZID","Log2FoldChange")
  #gsea_core_gene <- merge(gsea_core_gene,gsedata,by="ENTREZID",all.x=TRUE,all.y=FALSE)
  gsdata <- gsScore(object, geneSetID)

  #print(head(gsdata))
  msigdf_geneSet <- merge(msigdf_geneSet,gsdata,by.x="symbol",by.y="gene",all.x=TRUE,all.y=FALSE)
  #print(head(gsea_core_gene))
  names(msigdf_geneSet)[names(msigdf_geneSet) == 'x'] <- 'Rank'
  names(msigdf_geneSet)[names(msigdf_geneSet) == 'geneList'] <- 'Log2FoldChange'
  msigdf_geneSet$ymin <- NULL
  msigdf_geneSet$ymax <- NULL
  msigdf_geneSet$position <- NULL
  msigdf_geneSet$geneset <- NULL
  msigdf_geneSet <- na.omit(msigdf_geneSet[order(msigdf_geneSet$Rank),])
  msigdf_geneSet <- merge(msigdf_geneSet,gsea_core_gene,by.x="symbol",by.y="gene",all.x=TRUE,all.y=FALSE)
  msigdf_geneSet$core_enrich_gene[which(is.na(msigdf_geneSet$core_enrich_gene))] <- "NO"
  #msigdf_geneSet$core_enrich_gene <- which(gsea_core_gene$gene%in%msigdf_geneSet$symbol)
  return(msigdf_geneSet)
  #print(head(msigdf_geneSet))
}

GSEA_plot <- function (x, geneSetID, title = "", color = "green", 
    base_size = 11, rel_heights = c(1.5, 0.5, 1), subplots = 1:3) 
{
    ES_geom <- "line"
    geneList <- position <- NULL
    gsdata <- gsInfo(x, geneSetID)
	#print(head(gsdata))
	anno <- x[geneSetID, c('enrichmentScore','NES', 'pvalue', 'p.adjust')]
	#print(head(anno))
	colnames(anno) <- c('ES','NES','pvalue','p.adjust')
	lab <- paste0(names(anno), '=', round(anno, 3), collapse='\n')
	title = x$Description[geneSetID]
    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
        theme(panel.border = element_rect(linetype = "solid", color="black",fill=NA),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
        scale_x_continuous(expand = c(0, 0))
	
        es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
            size = 1)
	#print(p)
    p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
        legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))+ 
		geom_hline(yintercept = c(0),linetype="dotted")
	#print(p.res)
    p.res <- p.res + ylab("Enrichment Score(ES)") + 
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.line.x = element_blank(), plot.margin = margin(t = 0.2, 
                r = 0.2, b = 0, l = 0.2, unit = "cm"))
	
    i <- 0
    for (term in unique(gsdata$Description)) {
        idx <- which(gsdata$ymin != 0 & gsdata$Description == 
            term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
    }
	

    p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(size=0.1,aes_(ymin = ~ymin, 
        ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
        theme_classic(base_size) + theme(panel.border = element_rect(linetype = "solid", color="black",fill=NA),legend.position = "none", 
        plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
        axis.ticks = element_blank(), axis.text = element_blank(), 
        axis.line.x = element_blank()) + scale_x_continuous(expand = c(0, 
        0)) + scale_y_continuous(expand = c(0, 0))
	
        v <- seq(1, sum(gsdata$position), length.out = 9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) 
            inv <- inv + 1
        col = c(rev(brewer.pal(5, "Blues")), brewer.pal(5, 
            "Reds"))
        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * 0.2
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
            xmax = xmax, col = col[unique(inv)])
        p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
            ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
            inherit.aes = FALSE)+theme(panel.background = element_rect(fill="white"))
	
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
	cross_value <- max(which(df2$y>0))+1
	cross_text <- paste("Zero cross at",cross_value,sep=" ")
    p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
        y = ~y, yend = 0), color = "grey")+geom_vline(xintercept = cross_value,linetype="dotted")+annotate('text', as.numeric(cross_value+10), 1.3,label = cross_text, hjust=0, vjust=0)+
		annotate('text', as.numeric(0+max(which(df2$y<0))*0.02), as.numeric(max(df2$y)-1),label = "Positively correlated", hjust=0, vjust=0,color="red")+annotate('text',as.numeric(max(which(df2$y<0))-max(which(df2$y<0))*0.25), as.numeric(max(df2$y)-1),label = "Negatively correlated", hjust=0, vjust=0,color="blue")+geom_hline(yintercept = c(-1,1),linetype="dotted")
    p.pos <- p.pos + ylab("Ranked list metric(LFC)") + xlab("Rank in Ordered Dataset") + 
        theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
            l = 0.2, unit = "cm"))
    if (!is.null(title) && !is.na(title) && title != "") 
        p.res <- p.res + ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values = color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        }
        else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }
	

	
    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
        axis.ticks.x = element_line(), axis.text.x = element_text())
    if (length(subplots) == 1) 
        return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
            r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
    if (length(rel_heights) > length(subplots)) 
        rel_heights <- rel_heights[subplots]
    plot_grid(plotlist = plotlist, ncol = 1, align = "v", 
        rel_heights = rel_heights)+annotate('text', 0.83, 0.83,label = lab, hjust=0, vjust=0)
	
	
}


GSEA_BAR_PLOT <- function(object) {
    gsea_data_frame <- as.data.frame(object)
	#print(head(gsea_data_frame))
	gsea_bar_plot_data <- data.frame(Description=gsea_data_frame$ID,NES=gsea_data_frame$NES,padj=gsea_data_frame$p.adjust)
	#print(head(gsea_bar_plot_data))
	gsea_bar_plot_data$Description <- factor(gsea_bar_plot_data$Description, levels = gsea_bar_plot_data[order(-gsea_bar_plot_data$NES), 1])
	#print(head(gsea_bar_plot_data))
	ggplot(gsea_bar_plot_data, aes(Description, NES, fill = padj)) + 
	coord_flip() +
	theme_dose(12) +
	geom_bar(stat = 'identity') +
	xlab("") +ylab("Normalized Enrichment Score")+
	theme_bw() +theme(panel.grid =element_blank())+theme(panel.border = element_rect(size = 0.6))+
	scale_fill_continuous(low="red", high="blue", name = "p.adj", guide=guide_colorbar(reverse=TRUE))
}

GSEA_DOT_PLOT <- function(object){
    gsea_data_frame <- as.data.frame(object)
	gsea_dot_plot_data <- data.frame(Description=gsea_data_frame$ID,NES=gsea_data_frame$NES,padj=gsea_data_frame$p.adjust,core_enrichment=gsea_data_frame$core_enrichment)
	for(i in 1:dim(gsea_dot_plot_data)[1]){
	    gsea_dot_plot_data$core_number[i] <- length(unlist(strsplit(as.character(gsea_dot_plot_data$core_enrichment[i]),split="/")))
	}
	gsea_dot_plot_data$core_enrichment <- NULL
	#print(head(gsea_dot_plot_data))
	gsea_dot_plot_data$Description <- factor(gsea_dot_plot_data$Description, levels = gsea_dot_plot_data[order(-gsea_dot_plot_data$NES), 1])
	#labels=gsea_dot_plot_data$Description
    ggplot(gsea_dot_plot_data,aes(x = NES,y = Description))+
    geom_point(aes(x = NES,y = Description,color = padj,size=core_number))+
    theme_bw()+theme(panel.grid = element_blank())+
    theme_dose(12)+
    scale_color_gradient(low = "red", high = "blue")+
    scale_fill_continuous(low="red", high="blue", name = "p.adj", guide=guide_colorbar(reverse=TRUE))+
	theme(panel.border = element_rect(size = 0.6))+
    xlab("Normalized Enrichment Score")
    #dotplot(go_data,showCategory=go_number)
}

GSEA_HEATMAP_PLOT <- function(objectall,object,gsedata,metadata,gene_rld,geneSetID){
    gsea_data_frame <- as.data.frame(object)
	#print(head(gsea_data_frame))
    enrichment_name <- as.character(gsea_data_frame$ID[geneSetID])
    #print(enrichment_name)
    msigdf <- objectall
	#print(head(msigdf))
    msigdf_geneSet <- msigdf[which(msigdf$geneset==enrichment_name),]
    colnames(msigdf_geneSet) <- c("geneset","symbol")
    gsea_core_gene <- unlist(strsplit(gsea_data_frame$core_enrichment[geneSetID],split="/"))
    gsea_core_gene <- data.frame(gene=gsea_core_gene)
    gsea_core_gene$core_enrich_gene <- "YES"
    #print(head(msigdf_geneSet))
  
    msigdf_geneSet$ENSEMBL <- mapIds(spi_database,keys=as.character(msigdf_geneSet$symbol),column="ENSEMBL",keytype="SYMBOL",multiVals="first")
    msigdf_geneSet$ENTREZID <- mapIds(spi_database,keys=as.character(msigdf_geneSet$symbol),column="ENTREZID",keytype="SYMBOL",multiVals="first")
    msigdf_geneSet$GENENAME <- mapIds(spi_database,keys=as.character(msigdf_geneSet$symbol),column="GENENAME",keytype="SYMBOL",multiVals="first")
    #colnames(gsedata) <- c("ENTREZID","Log2FoldChange")
    #gsea_core_gene <- merge(gsea_core_gene,gsedata,by="ENTREZID",all.x=TRUE,all.y=FALSE)
    gsdata <- gsScore(object, geneSetID)

    #print(head(gsdata))
    msigdf_geneSet <- merge(msigdf_geneSet,gsdata,by.x="symbol",by.y="gene",all.x=TRUE,all.y=FALSE)
    #print(head(gsea_core_gene))
    names(msigdf_geneSet)[names(msigdf_geneSet) == 'x'] <- 'Rank'
    names(msigdf_geneSet)[names(msigdf_geneSet) == 'geneList'] <- 'Log2FoldChange'
    msigdf_geneSet$ymin <- NULL
    msigdf_geneSet$ymax <- NULL
    msigdf_geneSet$position <- NULL
    msigdf_geneSet$geneset <- NULL
    msigdf_geneSet <- na.omit(msigdf_geneSet[order(msigdf_geneSet$Rank),])
    msigdf_geneSet <- merge(msigdf_geneSet,gsea_core_gene,by.x="symbol",by.y="gene",all.x=TRUE,all.y=FALSE)
    msigdf_geneSet$core_enrich_gene[which(is.na(msigdf_geneSet$core_enrich_gene))] <- "NO"
    
	
	core_gene_number <- sum(msigdf_geneSet$core_enrich_gene=="YES")
	#print(msigdf_geneSet)
	
    gene_rld <- as.data.frame(gene_rld)
    gene_rld$symbol <- mapIds(spi_database,keys=row.names(gene_rld),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    
    gene_rld<-na.omit(gene_rld)
    gene_rld$symbol <- make.unique(as.character(gene_rld$symbol))
	
	gene_rld_tmp <- merge(msigdf_geneSet,gene_rld,by="symbol",all.x=TRUE,all.y=FALSE)
	if(gsea_data_frame$NES[geneSetID]>0){
	  gene_rld_tmp <- na.omit(gene_rld_tmp[order(gene_rld_tmp$Rank),])
	  gene_rld_tmp <- na.omit(gene_rld_tmp[order(gene_rld_tmp$core_enrich_gene,decreasing=T),])
	}else{
	  gene_rld_tmp <- na.omit(gene_rld_tmp[order(gene_rld_tmp$Rank,decreasing=T),])
	  gene_rld_tmp <- na.omit(gene_rld_tmp[order(gene_rld_tmp$core_enrich_gene,decreasing=T),])
	}
	gene_rld_tmp$symbol_expend <- paste(gene_rld_tmp$Rank,"-",gene_rld_tmp$symbol,"-",gene_rld_tmp$core_enrich_gene,sep="")
	gene_rld_tmp <- gene_rld_tmp[,-(2:9)]
    rownames(gene_rld_tmp)<-make.unique(as.character(gene_rld_tmp$symbol_expend))
    gene_rld_tmp$symbol_expend <- NULL
	gene_rld_tmp$symbol <- NULL
	#print(gene_rld_tmp)
    
	
    #gene_rld_sub<-subset(gene_rld,gene_rld$symbol%in%msigdf_geneSet$symbol)
    #print(gene_rld_sub)
	#gene_rld_sub$symbol <- NULL
    gene_rld_sub_mat <- gene_rld_tmp - rowMeans(gene_rld_tmp)
	pheatmap(gene_rld_sub_mat,annotation_col=metadata,fontsize=10,cluster_rows=FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cellwidth=15,cellheight=15)
    #print(gene_rld_sub_mat)
}

GSEA_SIGN_PLOT <- function(gsea_res_all){
    gsea_res_all$FDR <- -log(gsea_res_all$qvalues,10)
    twoord.plot(gsea_res_all$NES,gsea_res_all$p.adjust,gsea_res_all$NES,gsea_res_all$FDR,lcol="black",rcol="red",
            xlab="NES",ylab="p.adjust",rylab="-log10(qvalue)",
            main="Significance vs NES")
    legend("topright",legend=c("p.adjust","-log10(qvalue)"),pch=c(1,2),col=c(1,2))
}

GESA_NES_PLOT <- function(gsea_res_all){
    p1<-ggplot(data=gsea_res_all,aes(x=gsea_res_all$NES,..density..))+
	    geom_histogram(color='white',fill='gray60')+
		geom_line(stat='density')+
		theme_bw() +theme(panel.grid =element_blank())+
		ylab(label = '# of gene sets')+
		xlab(label = 'NES')+
		labs(title = 'NES_Gene_Sets_Number')+theme(plot.title = element_text(hjust = 0.5))
}

GSEA_marker_gene_plot <- function(marker_genes_data){
  marker_genes_data_top <- marker_genes_data[order(marker_genes_data$lfc),]
  marker_genes_data_bottom <- marker_genes_data[order(marker_genes_data$lfc,decreasing = TRUE),]
  
  marker_genes_data_top100 <- marker_genes_data_top[1:100,]
  marker_genes_data_top100$Rank <- 1:100
  marker_genes_data_top100$Type <- "Down"
  marker_genes_data_bottom100 <- marker_genes_data_bottom[1:100,]
  marker_genes_data_bottom100$Rank <- 1:100
  marker_genes_data_bottom100$Type <- "Up"
  #print(marker_genes_data_top100)
  #print(marker_genes_data_bottom100)
  marker_genes_data <- rbind(marker_genes_data_top100,marker_genes_data_bottom100)
  print(marker_genes_data)
  ggplot(marker_genes_data,aes(x = lfc,y = Rank))+
    geom_point(aes(color = Type))
}

RIDGEPLOT <- function(x, showCategory=30, fill="p.adjust",
                                 core_enrichment = TRUE, label_format = 30) {
    # has_package("ggridges")
	default_labeller <- function(n) {
    function(str){
        str <- gsub("_", " ", str)
        ep_str_wrap(str, n)
    }
}

ep_str_wrap <- function(string, width) {
    x <- gregexpr(' ', string)
    vapply(seq_along(x),
           FUN = function(i) {
               y <- x[[i]]
               n <- nchar(string[i])
               len <- (c(y,n) - c(0, y)) ## length + 1
               idx <- len > width
               j <- which(!idx)
               if (length(j) && max(j) == length(len)) {
                   j <- j[-length(j)]
               }
               if (length(j)) {
                   idx[j] <- len[j] + len[j+1] > width
               }
               idx <- idx[-length(idx)] ## length - 1
               start <- c(1, y[idx] + 1)
               end <- c(y[idx] - 1, n)
               words <- substring(string[i], start, end)
               paste0(words, collapse="\n")
           },
           FUN.VALUE = character(1)
    )
}

    if (!is(x, "gseaResult"))
        stop("currently only support gseaResult")

    ## fill <- match.arg(fill, c("pvalue", "p.adjust", "qvalue"))
    if (fill == "qvalue") {
        fill <- "qvalues"
    }
    if (!fill %in% colnames(x@result)) {
        stop("'fill' variable not available ...")
    }

    ## geom_density_ridges <- get_fun_from_pkg('ggridges', 'geom_density_ridges')

    n <- showCategory
    if (core_enrichment) {
        gs2id <- geneInCategory(x)[seq_len(n)]
    } else {
        gs2id <- x@geneSets[x$ID[seq_len(n)]]
    }

    gs2val <- lapply(gs2id, function(id) {
        res <- x@geneList[id]
        res <- res[!is.na(res)]
    })

    nn <- names(gs2val)
    i <- match(nn, x$ID)
    nn <- x$Description[i]

    j <- order(x$NES[i], decreasing=FALSE)

    len <- sapply(gs2val, length)
    gs2val.df <- data.frame(category = rep(nn, times=len),
                            color = rep(x[i, fill], times=len),
                            value = unlist(gs2val))

    colnames(gs2val.df)[2] <- fill
    gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])

    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }
    #print(head(gs2val))
    ggplot(gs2val.df, aes_string(x="value", y="category", fill=fill)) +
        ggridges::geom_density_ridges() +
        ## scale_x_reverse() +
        scale_fill_continuous(low="red", high="blue", name = fill,
            guide=guide_colorbar(reverse=TRUE)) +
        scale_y_discrete(labels = label_func) +
        ## scale_fill_gradientn(name = fill, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
        ## geom_vline(xintercept=0, color='firebrick', linetype='dashed') +
        xlab(NULL) + ylab(NULL) +  theme_dose(8)+
		xlab("log2(Fold change)")
}

gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList
    #print(head(geneList))
    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
	#print(head(geneSet))
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
	#print(df)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
	#print(head(geneList))
    df$geneList <- geneList
	#df$gene <- names(geneList)
	#print(head(df))
    #print(head(object@result))
    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}

gseaScores <- getFromNamespace("gseaScores", "DOSE")

tableGrob2 <- function(d, p = NULL) {
    # has_package("gridExtra")
    d <- d[order(rownames(d)),]
    tp <- gridExtra::tableGrob(d)
    if (is.null(p)) {
        return(tp)
    }

    # Fix bug: The 'group' order of lines and dots/path is different
    p_data <- ggplot_build(p)$data[[1]]
    # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
    p_data <- p_data[order(p_data[["group"]]), ]
    pcol <- unique(p_data[["colour"]])
    ## This is fine too
    ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
    j <- which(tp$layout$name == "rowhead-fg")

    for (i in seq_along(pcol)) {
        tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
    }
    return(tp)
}


gsScore <- function(object, geneSetID) {
    geneList <- object@geneList
    #print(head(geneList))
    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
	#print(head(geneSet))
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
	#print(df)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
	#print(head(geneList))
    df$geneList <- geneList
	df$gene <- names(geneList)
	#print(head(df))
    #print(head(object@result))
    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}
#map_rate和region柱状图绘制
#ggplot(map_rate,aes(x=type,y=reads_number,group=times))+ geom_bar(stat="identity", aes(fill=times)) +theme_bw()+theme(panel.grid = element_blank())

#单样本绘制PCA和层次聚类图
#pca <- prcomp(t(log_gene_fpkm_matrix), scale=FALSE)
#pca.data <- data.frame(Sample=rownames(pca$x),
#+                        X=pca$x[,1],
#+                        Y=pca$x[,2])
#pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
#> pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
#ggplot(data=pca.data,aes(x=X,y=Y,color=Sample))+
#+     geom_point(size=3)+
#+     theme_bw()+theme(panel.grid=element_blank())+
#+     xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
#+     ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))

#sampleDist <- dist(t(log_gene_fpkm_matrix))
#sampleDistMatrix <- as.matrix(sampleDist)
#cluster_plot <- pheatmap(sampleDistMatrix,display_numbers = FALSE,number_color = "black",color = colorRampPalette(c("blue","white"))(256))
