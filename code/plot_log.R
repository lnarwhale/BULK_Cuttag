.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
#------read in------------------
library(getopt)
library(ggbiplot)
library(ggstatsplot)
library(magrittr)
library(RColorBrewer)
spec=matrix(c('log_dir','l',1,'character',
              'name','n',1,'character'),
            byrow = TRUE,ncol = 4)
opt=getopt(spec)
setwd(opt$log_dir)

pic_name <- opt$name

a <- read.csv(paste0(pic_name,"_bamqc.txt"),header = FALSE)
value <- as.numeric(a$V2)
name <- a$V1
png(paste0(pic_name,"_barplot.png"))
barplot(height = value,names.arg = name)
dev.off()

bing <- a[1:2,1:2]
colnames(bing) <- c('names','freq')
bing[2,1] <- c('nodup')
bing[3,1] <- c('dup')
bing[3,2] <- value[1]-value[2]
bing <- bing[-1,]
percent <- round(bing[,2]/sum(bing[,2])*100,1)
label <- paste(bing[,1],"(",percent,"%)")
png(filename = paste0(pic_name,"_pieplot.png"),height = 400,width = 400)
pie(bing[,2],border = "white",col=brewer.pal(5,"Set3"),label = label)
dev.off()


