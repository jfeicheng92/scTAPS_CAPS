##This code is to make figure2 to demonstrate the performance of scTAPS and scCAPS on mESC
##module load R/4.0.3-foss-2020b

.libPaths("/well/ludwig/users/ebu571/R/4.0/skylake")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
library(readr)
options(bitmapType='cairo-png')
setwd("/users/ludwig/ebu571/ebu571/project/single_cell/scTAPS_CAPS_RNA")

raw_dat <- read.table("meth_update/read_summary.txt",header=TRUE)
raw_dat <- raw_dat[raw_dat$sample!="scTAPS203_N702_i7_10_rev_N502_i5_1_rev",]
mapping_jf <- read.table("meth_update/mapping_jf.txt",header=TRUE)


##plot conversion and false positive

smps <- mapping_jf$sample
res1 <- NULL
for(smp in smps){
  input_file <- paste0("/users/ludwig/cfo155/cfo155/scTAPS_CAPS/mESC_update/align/",smp,".spikeins.filter.meth.sta.txt.gz")
  dat <- suppressMessages(fread(paste0('zcat ', input_file), header=TRUE))
  mC <- sum(dat[dat$chr=="J02459.1",3])/sum(dat[dat$chr=="J02459.1",4])*100 %>% round(2)
  hmC <- sum(dat[dat$chr=="144hmC" & dat$pos %in% c(86, 91, 99, 109),3])/sum(dat[dat$chr=="144hmC" & dat$pos %in% c(86, 91, 99, 109),4]) *100 %>% round(2)
  uC <- sum(dat[dat$chr=="unmodified_2kb",3])/sum(dat[dat$chr=="unmodified_2kb",4])*100 %>% round(2)
  
  mC_mC <-   sum(dat[dat$chr=="J02459.1",3])
  mC_aC <- sum(dat[dat$chr=="J02459.1",4])
  hmC_mC <- sum(dat[dat$chr=="144hmC" & dat$pos %in% c(86, 91, 99, 109),3])
  hmC_aC <- sum(dat[dat$chr=="144hmC" & dat$pos %in% c(86, 91, 99, 109),4]) 
  uC_mC <- sum(dat[dat$chr=="unmodified_2kb",3])
  uC_aC <- sum(dat[dat$chr=="unmodified_2kb",4])
  res1 <- rbind(res1, c(smp, mC, hmC, uC, mC_mC, mC_aC, hmC_mC, hmC_aC, uC_mC, uC_aC))
}
res1 <- as.data.frame(res1)
colnames(res1) <- c("smp","mC","hmC","uC", "mC_mC", "mC_aC", "hmC_mC", "hmC_aC", "uC_mC", "uC_aC")
write.table(res1, "plot/sctaps_caps.conversion.txt", col.names = TRUE, quote=F)
res1[!grepl("scTAPS", res1$smp),"lib"] <- "scCAPS"
res1[,-c(1:4,11)] <-  lapply(res1[,-c(1:4,11)], as.numeric)
res1 <- aggregate(. ~ lib, res1[,-c(1:4)], sum)
res1$mC <- res1$mC_mC/res1$mC_aC
res1$hmC <- res1$hmC_mC/res1$hmC_aC
res1$uC <- res1$uC_mC/res1$uC_aC
res <- melt(res1[,-c(2:7)])


res$variable <- factor(res$variable, level=c("hmC", "mC","uC"))
pdf("plot/sctaps_caps_conversion_caps.pdf", width=4,height=4)
ggplot(res[res$lib=="scCAPS",],aes(x=variable,y=value*100,fill=lib))+ geom_bar(stat="identity") +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=2), size=3.5) +
  theme_bw() +
  ylab("modification level%") +
  xlab("lib")+
  theme(legend.position="None") +
  scale_fill_manual(values=c("#c7e1a6")) +scale_y_continuous(limits=c(-1,100),expand = c(0, 0))
dev.off()











####plot mapping and genomic CpG coverage
####meth/read_summary.txt 
mapping_merged <- merge(raw_dat, mapping_jf, by=c("sample"))
mapping_merged$q10_pmap <- as.numeric(mapping_merged$q10_nmap/mapping_merged$nclean_reads)
mapping_merged[grepl("scTAPS", mapping_merged$sample),"lib"] <- "scTAPS"
mapping_merged[!grepl("scTAPS", mapping_merged$sample),"lib"] <- "scCAPS"

cov_dat <- read.table("output/genome_coverage_update.txt")
colnames(cov_dat) <- c("sample", "CpG_num", "per_base")
seldat <- cov_dat[!grepl("bulk_taps",cov_dat$sample) & !grepl("bulk_caps",cov_dat$sample) & !grepl("cells_caps",cov_dat$sample) & !grepl("cells_taps",cov_dat$sample),]
seldat$sample <- gsub(".md.filter.meth.sta.txt.gz", "", seldat$sample)
seldat_mapping <- merge(seldat, mapping_merged, by=c("sample"))


p3 <- ggplot(seldat_mapping[seldat_mapping$lib=="scCAPS",],aes(x=lib,y=q10_pmap*100,fill=lib))+ geom_violin() +
  geom_jitter(size=0.5, width = 0.1) +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=8), size=5)+
  theme_bw() +
  ylab("q10 mapping rate%") +
  xlab("lib")+
  theme(legend.position="None") +scale_fill_manual(values=c("#c7e1a6"))+scale_y_continuous(limits=c(0,100),expand = c(0, 0))

p4 <- ggplot(seldat_mapping[seldat_mapping$lib=="scCAPS",],aes(x=lib,y=round(CpG_num/21342780,3)*100,fill=lib))+ geom_violin() +
  geom_jitter(size=0.5, width = 0.1) +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=8), size=3.5)+
  theme_bw() +
  ylab("genomic CpG coverage%") +
  xlab("lib")+
  theme(legend.position="None") +scale_fill_manual(values=c("#c7e1a6"))+scale_y_continuous(limits=c(0,25),expand = c(0, 0))

level_dat <- read.table("output/meth_level_update.txt", header=TRUE)
p_p3 <- ggplot(level_dat,aes(x="",y=mean_level*100,fill=""))+ geom_violin() +
  geom_jitter(size=0.5, width = 0.1) +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=8), size=5)+
  theme_bw() +
  ylab("aggregated 5hmC level") +
  xlab("scCAPS")+
  theme(legend.position="None") +scale_fill_manual(values=c("#c7e1a6"))+scale_y_continuous(limits=c(0,20),expand = c(0, 0))

library(cowplot)
p <- plot_grid(p3,p4,p_p3, labels = c('A', 'B','C'), ncol=3,label_size = 8)
ggsave(p, width =6, height = 4,filename = "plot/sctaps_caps_mapping_caps.pdf")

####raw reads and CpG coverage

pdf("plot/sctaps_caps_CpG_perbase_cov_sccaps.pdf")
ggplot(seldat_mapping[seldat_mapping$lib=="scCAPS",], aes(x=nclean_reads/1000000, y=as.numeric(CpG_num)/1000000,fill=lib)) +
  geom_smooth(method = "lm", se = FALSE,color="grey40")+
  geom_point(aes(size=1, color=lib))+theme_bw()+scale_color_manual(values=c("#c7e1a6"))+scale_y_continuous(limits=c(0,5),expand = c(0, 0))+scale_x_continuous(limits=c(0,20),expand = c(0, 0))
  dev.off()

pdf("plot/sctaps_caps_CpG_perbase_cov_sccaps_no_legend.pdf",width=4,height=4)
ggplot(seldat_mapping[seldat_mapping$lib=="scCAPS",], aes(x=nclean_reads/1000000, y=as.numeric(CpG_num)/1000000,fill=lib)) +
  geom_smooth(method = "lm", se = FALSE,color="grey40")+
  geom_point(aes(size=1, color=lib))+theme_bw()+scale_color_manual(values=c("#c7e1a6"))+scale_y_continuous(limits=c(0,5),expand = c(0, 0))+scale_x_continuous(limits=c(0,20),expand = c(0, 0))+theme(legend.position="None") 
  dev.off()

####raw signal correlation
window_size <- 100000
cov_min <- 300
cov_max <- 30000


merge_sccaps <- read.table("meth_update/merged_sccaps.bed.gz", header=F, sep="\t")
colnames(merge_sccaps) <- c("chr", "start", "end", "mC", "aC", "ratio")
sd_caps <- read.table("meth_old/standard_caps_merge_CpG.bedGraph.gz", header=T,comment.char = "")
sd_caps$chr <- sd_caps$X.chr
sd_caps <- sd_caps[grepl("chr", sd_caps$chr),] %>% dplyr::select(-X.chr)

merge_sccaps <- merge_sccaps[merge_sccaps$chr%in%paste0("chr",seq(1,19)),]
merge_sccaps$idx <- round(merge_sccaps$end/window_size)*window_size
merge_sccaps_a <- merge_sccaps %>% group_by(chr, idx) %>% summarize(mC=sum(mC), aC=sum(aC))
merge_sccaps_a$rC <- merge_sccaps_a$mC/merge_sccaps_a$aC
sum(merge_sccaps$mC)/sum(merge_sccaps$aC) ##0.05896677

sd_caps <- sd_caps[sd_caps$chr%in%paste0("chr",seq(1,19)),]
sd_caps$idx <- round(sd_caps$end/window_size)*window_size
sd_caps_a <- sd_caps %>% group_by(chr, idx) %>% summarize(mC=sum(mC), aC=sum(total))
sd_caps_a$rC <- sd_caps_a$mC/sd_caps_a$aC
sum(sd_caps$mC)/sum(sd_caps$total) ##0.0437825

dat <- merge(merge_sccaps_a, sd_caps_a, by=c("chr","idx"))
pdf("test.pdf")
hist(dat$aC.x, xlim=c(0,30000),breaks=4000)
dev.off()

dat <- dat[complete.cases(dat) & dat$aC.x<cov_max&dat$aC.x>cov_min,]

pdf("plot/merge_sccaps_raw_correlation.pdf")
  smoothScatter(x=dat$rC.x*100, y=dat$rC.y*100, 
                xlab=paste0("merged sccaps\n cor:", round(cor(dat$rC.x, dat$rC.y),3)), 
                ylab="standard caps", xlim = c(0, 100), ylim = c(0,100))
dev.off()



######correlate with standard CAPS, pick three cells
window_size <- 100000
cov_min <- 50
cov_max <- 30000
cor_caps <- NULL
smps <- c("meth_old/10k_cells_caps.md.filter.meth.sta.txt.gz","meth_update/merged_sccaps.bed.gz", "meth_update/C203_N702_i7_10_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz", "meth_update/C203_N702_i7_7_rev_N502_i5_1_rev.md.filter.meth.sta.txt.gz", "meth_update/C203_N702_i7_8_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz")
#smps <- c("meth_old/standard_caps_merge_CpG.bedGraph.gz")
#j <-  c("meth_update/C203_N702_i7_10_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz")
for(smp in smps){
  for (j in smps){
  library_type <- "caps"
  if(grepl("merged", smp)){
    input <- read.table(smp,header=F, sep="\t")
    colnames(input) <- c("chr", "start", "end", "mC", "aC", "ratio")
  }else if(grepl("standard", smp)){
    input <- read.table(smp,header=F, sep="\t")
    colnames(input) <- c("chr", "start", "end", "ratio","mC", "aC")    
    }else{
    input <- read.table(smp,header=T, sep="\t")
    colnames(input) <- c("chr", "end", "mC", "aC", "ratio")
  }
  input <- input[input$chr%in%paste0("chr",seq(1,19)),]
  input$idx <- round(input$end/window_size)*window_size
  input_a <- input %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC)) 
  input_a$rC <- input_a$mC/input_a$aC
  sum(input$mC)/sum(input$aC)
  
  if(grepl("merged", j)){
    input2 <- read.table(j,header=F, sep="\t")
    colnames(input2) <- c("chr", "start", "end", "mC", "aC", "ratio")
  }else if(grepl("standard", j)) {
    input2 <- read.table(j,header=F, sep="\t")
    colnames(input2) <- c("chr", "start", "end", "ratio","mC", "aC")    
    }else{
    input2 <- read.table(j,header=T, sep="\t")
    colnames(input2) <- c("chr", "end", "mC", "aC", "ratio")
  }
  input2 <- input2[input2$chr%in%paste0("chr",seq(1,19)),]
  input2$idx <- round(input2$end/window_size)*window_size
  input2_a <- input2 %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC)) 
  input2_a$rC <- input2_a$mC/input2_a$aC
  sum(input2$mC)/sum(input2$aC)
  
  dat <- merge(input_a, input2_a, by=c("chr","idx"))
  dat <- dat[complete.cases(dat) & dat$aC.x<cov_max&dat$aC.x>cov_min &dat$aC.y<cov_max&dat$aC.y>cov_min,]
  cor_caps <- rbind(cor_caps, c(smp, j,round(cor(dat$rC.x, dat$rC.y),2)))

}
}

cor_caps <- as.data.frame(cor_caps)
colnames(cor_caps) <- c("smp1", "smp2", "cor_value")
cor_caps$smp1 <- factor(cor_caps$smp1, level=c("meth_old/10k_cells_caps.md.filter.meth.sta.txt.gz", "meth_update/merged_sccaps.bed.gz", "meth_update/C203_N702_i7_10_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz","meth_update/C203_N702_i7_7_rev_N502_i5_1_rev.md.filter.meth.sta.txt.gz", "meth_update/C203_N702_i7_8_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz"))
cor_caps$smp2 <- factor(cor_caps$smp2, level=c("meth_old/10k_cells_caps.md.filter.meth.sta.txt.gz", "meth_update/merged_sccaps.bed.gz", "meth_update/C203_N702_i7_10_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz","meth_update/C203_N702_i7_7_rev_N502_i5_1_rev.md.filter.meth.sta.txt.gz", "meth_update/C203_N702_i7_8_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz"))

write.table(cor_caps, "output/sccaps.correlation.table", quote = FALSE, col.names = F, row.names = F)


get_upper_tri <- function(CorMat){
  CorMat[upper.tri(CorMat)]<- NA
  return(CorMat)
}

library(corrplot)
library(tidyr)
library(reshape2)
cot_value <- spread(cor_caps,smp2,cor_value) %>% as.data.frame()
cot_value <- data.frame(cot_value[,-1], row.names=as.character(cot_value[,1])) %>% as.matrix(rownames = TRUE,colnames=TRUE)
cor_caps1 <- get_upper_tri(cot_value)
meltNum <- melt(cor_caps1, na.rm = T)
meltNum$value <- as.numeric(meltNum$value)
pdf("plot/heatmap.correlation_sc_caps.pdf", width=10, height=10)
ggplot(meltNum, aes(x = Var1, y = Var2, fill = value)) + geom_tile()+coord_fixed()+
  geom_text(aes(label = round(as.numeric(value),2)), color = "black")+  scale_fill_gradient(low = "white", high = "#B2182B",
                      limit = c(0,1), name = "Pearson\nCorrelation") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
dev.off()






####heatmap for scCAPS ###chr12:3,000,000-120,000,000
bin_size <- 1000000 ##1MB

merged_sccaps <- read.table("meth_update/merged_sccaps.bed.gz", header=F, sep="\t")
colnames(merged_sccaps) <- c("chr", "start", "end", "mC", "aC", "ratio")

merged_sccaps <- merged_sccaps[merged_sccaps$chr%in%paste0("chr12") & merged_sccaps$end>=3000000 &merged_sccaps$end<=120000000 ,]
merged_sccaps$idx <- round(merged_sccaps$end/bin_size)*bin_size
merged_sccaps_a <- merged_sccaps %>% 
  group_by(chr, idx) %>% 
  summarize(mC = sum(mC), aC=sum(aC)) 
merged_sccaps_a$rC <- merged_sccaps_a$mC/merged_sccaps_a$aC

dat_mer <- merged_sccaps_a[,c("chr", "idx","rC")]
for(smp in smps){
  input_file <- paste0("meth_update/",smp,".md.filter.meth.sta.txt.gz")
  input <- read.table(input_file,header=T, sep="\t")
  colnames(input) <- c("chr", "end", "mC", "aC", "ratio")
  input <- input[input$chr%in%paste0("chr12") & input$end>=3000000 &input$end<=120000000 ,]
  input$idx <- round(input$end/bin_size)*bin_size
  input_a <- input %>% group_by(chr, idx) %>% summarize(mC=sum(mC), aC=sum(aC))
  input_a$rC <- input_a$mC/input_a$aC
  dat_mer <- merge(dat_mer, input_a[,c("chr", "idx","rC")], by=c("chr", "idx"))
}

nrow(dat_mer)
write.table(dat_mer, "output/sccaps.table", quote = FALSE, col.names = F, row.names = F)
#dat_mer <- read.table("output/sccaps.table")
#dat_mer$uni_id <- seq.int(nrow(dat_mer))
dat_mer$uni_id <- dat_mer$V2/1000000
dat_mer$uni_id <- factor(dat_mer$uni_id, level=seq.int(3,120))
dat_mer <- dat_mer[,-3] %>% dplyr::select(-V1, -V2)
dat_mer <- dat_mer[order(dat_mer$uni_id,decreasing=FALSE),]

final_matrix <- data.frame(dat_mer[,-97], row.names=as.character(dat_mer[,97])) %>% as.matrix(rownames = TRUE)
library(pheatmap)
library(RColorBrewer)
pdf("plot/sccaps_heatmap.pdf",width=24, height=12)
pheatmap(t(final_matrix),cluster_cols = FALSE, cluster_rows = FALSE,border_color=NA,show_rownames=FALSE, show_colnames=FALSE,color=colorRampPalette(c("navy", "white", "red"))(50))
dev.off()




####plot average level ---violin
level_dat <- read.table("output/meth_level_update.txt", header=TRUE)
p <- ggplot(level_dat,aes(x="",y=mean_level*100,fill=""))+ geom_violin() +
  geom_jitter(size=0.5, width = 0.1) +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=8), size=5)+
  theme_bw() +
  ylab("aggregated 5hmC level") +
  xlab("scCAPS")+
  theme(legend.position="None") +scale_fill_manual(values=c("#c7e1a6"))+scale_y_continuous(limits=c(0,100),expand = c(0, 0))
ggsave(p, "plot/mESC_sccaps.level.pdf")



#######Plot coverage over CGI for standard CAPS and merged bulk 
sd_caps_cgi <- read_delim("meth_old/standard_caps_merge_CpG_CGI.bed",delim="\t",col_names=FALSE)
colnames(sd_caps_cgi) <- c("gchr", "gstart", "gend", "index", "chr","start","end","mC","aC","ratio")
sd_caps_cgi <- sd_caps_cgi[,c("index", "aC")]

sc_caps<- read_delim("meth_update/merged_sccaps_CGI.bed",delim="\t",col_names=FALSE)
colnames(sc_caps) <- c("gchr", "gstart", "gend", "index", "chr","start","end","mC","aC","ratio")
sc_caps <- sc_caps[,c("index", "aC")]


###check coverage
sum(sd_caps_cgi$aC)/nrow(sd_caps_cgi) ##22.28995
sum(sc_caps$aC)/nrow(sc_caps) ##11.76267
sf <- (sum(sc_caps$aC)/nrow(sc_caps) )/(sum(sd_caps_cgi$aC)/nrow(sd_caps_cgi))

dat <- merge(aggregate(aC ~ index,sd_caps_cgi,mean),aggregate(aC ~ index,sc_caps,mean), by=c("index"))
dat$aC.sf <- dat$aC.x*sf

mel_CpG <- melt(dat[,c("index", "aC.sf", "aC.y")], id.vars =c("index"))

p1 <- ggplot(mel_CpG, aes(x=index, y=value, color=variable)) +
  geom_line()+
  theme(legend.position = "bottom") +
  #ylim(5.5,10.5) +
  ylab("coverage") +
  xlab("around CGI") +
  scale_x_continuous(breaks=c(0,20,30,50),
                     labels=c("CpG-4k", "                CpG island","", "CpG+4k"))+
  scale_color_manual(values=c("#2D6A4E","#c7e1a6"),name="")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
ggsave("plot/sccaps_coverage_around_CGI_0.60.pdf", p1, width=5,height = 3)

####plot CGI coverage for sctaps

sd_taps_cgi <- read_delim("meth_update/cd8_tcells_merge_CGI.bed",delim="\t",col_names=FALSE)
colnames(sd_taps_cgi) <- c("gchr", "gstart", "gend", "index", "chr","start","end","aC")
sd_taps_cgi <- sd_taps_cgi[,c("index", "aC")]

sc_taps<- read_delim("meth_update/C183_bulk_CGI.bed",delim="\t",col_names=FALSE)
colnames(sc_taps) <- c("gchr", "gstart", "gend", "index", "chr","start","end","aC")
sc_taps <- sc_taps[,c("index", "aC")]

###check coverage
sum(sd_taps_cgi$aC)/nrow(sd_taps_cgi) ##7.762967
sum(sc_taps$aC)/nrow(sc_taps) ##8.577623
sf <- (sum(sc_taps$aC)/nrow(sc_taps) )/(sum(sd_taps_cgi$aC)/nrow(sd_taps_cgi))

dat_taps <- merge(aggregate(aC ~ index,sd_taps_cgi,mean),aggregate(aC ~ index,sc_taps,mean), by=c("index"))
dat_taps$aC.sf <- dat_taps$aC.x*sf

mel_CpG_taps <- melt(dat_taps[,c("index", "aC.sf", "aC.y")], id.vars =c("index"))

p1 <- ggplot(mel_CpG_taps, aes(x=index, y=value, color=variable)) +
  geom_line()+
  theme(legend.position = "bottom") +
  ylim(5,10) +
  ylab("coverage") +
  xlab("around CGI") +
  scale_x_continuous(breaks=c(0,20,30,50),
                     labels=c("CpG-4k", "                CpG island","", "CpG+4k"))+
  scale_color_manual(values=c("#0066CC","#A3DEF4"),name="")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
ggsave("plot/sctaps_coverage_around_CGI.pdf", p1, width=5,height = 3)






#####compared 5hmC level around different chroHMM state

smps <- mapping_jf$sample
hmm_dat <- NULL
for(smp in smps){
input_file <- paste0("meth_update/",smp,".md.filter.meth.sta.HMM.bedgraph.gz")
cell_hmm <- read.table(input_file, header=FALSE)
colnames(cell_hmm) <- c("chr", "start", "end", "status", "mC", "aC")
cell_hmm$status <- gsub('[0-9]+_', '', cell_hmm$status)
active_pro <- sum(cell_hmm[cell_hmm$status=="Active_Promoter"  ,]$mC %>% as.numeric() %>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Active_Promoter",]$aC %>% as.numeric()%>%na.omit()) 
poised_pro <- sum(cell_hmm[cell_hmm$status=="Poised_Promoter",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Poised_Promoter",]$aC %>% as.numeric()%>%na.omit()) 
str_enhancer <- sum(cell_hmm[cell_hmm$status=="Strong_Enhancer",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Strong_Enhancer",]$aC %>% as.numeric()%>%na.omit()) 
poised_enhancer <- sum(cell_hmm[cell_hmm$status=="Poised_Enhancer",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Poised_Enhancer",]$aC %>% as.numeric()%>%na.omit()) 
txn_translation <- sum(cell_hmm[cell_hmm$status=="Txn_Transition",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Txn_Transition",]$aC %>% as.numeric()%>%na.omit())
txn_enlong <- sum(cell_hmm[cell_hmm$status=="Txn_Elongation",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Txn_Elongation",]$aC %>% as.numeric()%>%na.omit())
weak_txn <- sum(cell_hmm[cell_hmm$status=="Weak_Txn",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Weak_Txn",]$aC %>% as.numeric()%>%na.omit())
insulate <- sum(cell_hmm[cell_hmm$status=="Insulator",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Insulator",]$aC %>% as.numeric()%>%na.omit())
repres <- sum(cell_hmm[cell_hmm$status=="Repressed",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Repressed",]$aC %>% as.numeric()%>%na.omit())
heter <- sum(cell_hmm[cell_hmm$status=="Heterochrom",]$mC %>% as.numeric()%>%na.omit()) /sum(cell_hmm[cell_hmm$status=="Heterochrom",]$aC %>% as.numeric()%>%na.omit())
hmm_dat <- rbind(hmm_dat, c(smp, active_pro, poised_pro, str_enhancer,poised_enhancer,txn_translation,txn_enlong,weak_txn,insulate,repres,heter))
}

hmm_dat <- as.data.frame(hmm_dat)
colnames(hmm_dat) <- c("smp","Active_Promoter","Poised_Promoter","Strong_Enhancer", "Poised_Enhancer", "Txn_Transition", "Txn_Elongation", "Weak_Txn", "Insulator", "Repressed", "Heterochrom")
melt_anno <- melt(hmm_dat, id.vars=c("smp"))
melt_anno$variable <- factor(melt_anno$variable, levels = c("Active_Promoter","Poised_Promoter","Strong_Enhancer", "Poised_Enhancer", "Txn_Transition", "Txn_Elongation", "Weak_Txn", "Insulator", "Repressed", "Heterochrom"))

p <- ggplot(melt_anno, aes(x = variable, y = as.numeric(value)*100, fill = variable)) + 
  geom_violin() + geom_jitter(size=0.5, width = 0.1) +
  theme_light() +   stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.2f", ..y..)),
               position=position_nudge(y=8), size=5)+
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8, size = 12,colour = "black"),
        axis.text.y = element_text(colour = "black"))  +scale_y_continuous(limits=c(0,20),expand = c(0, 0))+
  #scale_fill_manual(values = c("#2D6A4E","#c7e1a6")) + 
  xlab("") + ylab("level")

ggsave("plot/scCAPS_chromHMM.pdf", p, width = 10, height = 6)
write.table(hmm_dat, "plot/scCAPS_chromHMM.data.txt", row.names = F, col.names = T, sep = "\t", quote = F )
