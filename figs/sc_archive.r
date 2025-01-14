##This code is to make figure2 to demonstrate the performance of scTAPS and scCAPS on mESC
##module load R/4.0.3-foss-2020b

#.libPaths("/well/ludwig/users/ebu571/R/4.0/skylake")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
library(readr)
setwd("/users/ludwig/ebu571/ebu571/project/single_cell/scTAPS_CAPS_RNA")

raw_dat <- read.table("meth_update/read_summary.txt",header=TRUE)
raw_dat <- raw_dat[raw_dat$sample!="scTAPS203_N702_i7_10_rev_N502_i5_1_rev",]
mapping_jf <- read.table("meth_update/mapping_jf.txt",header=TRUE)



##Fig1b

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


###Fig1c
#####plot for per-base coverage and CpG genome coverage in low-input samples
cov_dat <- read.table("output/genome_coverage.txt")
colnames(cov_dat) <- c("sample", "CpG_num", "per_base")
seldat <- cov_dat[grepl("bulk_taps",cov_dat$sample) | grepl("bulk_caps",cov_dat$sample) | grepl("cells_caps",cov_dat$sample) | grepl("cells_taps",cov_dat$sample),]
seldat[grepl("taps",seldat$sample),c("group")]<- "taps"
seldat[grepl("caps",seldat$sample),c("group")] <- "caps"
seldat[grepl("bulk",seldat$sample),c("cell_num")]<- "bulk" 
seldat[grepl("10",seldat$sample),c("cell_num")]<- "10" 
seldat[grepl("100",seldat$sample),c("cell_num")]<- "100" 
seldat[grepl("1k",seldat$sample),c("cell_num")]<- "1k" 
seldat[grepl("10k",seldat$sample),c("cell_num")]<- "10k" 

library(ggplot2)
pdf("plot/low_input_CpG_perbase_cov.pdf")
ggplot(seldat, aes(x=per_base, y=round(CpG_num/21342780,2)*100, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group, shape=cell_num), size=2)+theme_bw()+scale_y_continuous(limits=c(0,100),expand = c(0, 0))+scale_x_continuous(limits = c(0,25),expand = c(0, 0))+scale_color_manual(values=c("#1B9E77", "#D95F02"))
  dev.off()

pdf("plot/low_input_CpG_perbase_cov_no_legend.pdf", width=4, height=4)
ggplot(seldat, aes(x=per_base, y=round(CpG_num/21342780,2)*100, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group, shape=cell_num), size=2)+theme_bw()+scale_y_continuous(limits=c(0,100),expand = c(0, 0))+scale_x_continuous(limits = c(0,25),expand = c(0, 0))+scale_color_manual(values=c("#1B9E77", "#D95F02"))+ theme(legend.position = "none")
  dev.off()




###Fig1d
######correlate with standard CAPS, pick three cells that have highest number of nclean
window_size <- 100000
cov_min <- 50
cov_max <- 30000
cor_caps <- NULL
smps <- c("meth_old/10k_cells_caps.md.filter.meth.sta.txt.gz","meth_update/merged_sccaps.bed.gz", "meth_update/C203_N702_i7_5_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz", "meth_update/C203_N702_i7_7_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz", "meth_update/C203_N702_i7_4_rev_N502_i5_6_rev.md.filter.meth.sta.txt.gz")
#smps <- c("meth_old/standard_caps_merge_CpG.bedGraph.gz")
#j <-  c("meth_update/C203_N702_i7_5_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz")
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
cor_caps$smp1 <- factor(cor_caps$smp1, level=c("meth_old/10k_cells_caps.md.filter.meth.sta.txt.gz", "meth_update/merged_sccaps.bed.gz", "meth_update/C203_N702_i7_5_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz","meth_update/C203_N702_i7_7_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz", "meth_update/C203_N702_i7_4_rev_N502_i5_6_rev.md.filter.meth.sta.txt.gz"))
cor_caps$smp2 <- factor(cor_caps$smp2, level=c("meth_old/10k_cells_caps.md.filter.meth.sta.txt.gz", "meth_update/merged_sccaps.bed.gz", "meth_update/C203_N702_i7_5_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz","meth_update/C203_N702_i7_7_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz", "meth_update/C203_N702_i7_4_rev_N502_i5_6_rev.md.filter.meth.sta.txt.gz"))

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



###Fig1e
####heatmap for scCAPS ###chr12:3,000,000-120,000,000
bin_size <- 200000##1MB

merged_sccaps <- read.table("meth_update/merged_sccaps.bed.gz", header=F, sep="\t")
colnames(merged_sccaps) <- c("chr", "start", "end", "mC", "aC", "ratio")
##chr12:81,000,000-111,000,000
merged_sccaps <- merged_sccaps[merged_sccaps$chr%in%paste0("chr12") & merged_sccaps$end>=81000000 &merged_sccaps$end<=111000000,]
merged_sccaps$idx <- round(merged_sccaps$end/bin_size)*bin_size
merged_sccaps_a <- merged_sccaps %>% 
  group_by(chr, idx) %>% 
  summarize(mC = sum(mC), aC=sum(aC)) 
merged_sccaps_a$rC <- merged_sccaps_a$mC/merged_sccaps_a$aC

mapping_jf <- read.table("meth_update/mapping_jf.txt",header=TRUE)

smps <- mapping_jf$sample

dat_mer <- merged_sccaps_a[,c("chr", "idx","rC")]
for(smp in smps){
  input_file <- paste0("meth_update/",smp,".md.filter.meth.sta.txt.gz")
  input <- read.table(input_file,header=T, sep="\t")
  colnames(input) <- c("chr", "end", "mC", "aC", "ratio")
  input <- input[input$chr%in%paste0("chr12") & input$end>=81000000 &input$end<=111000000 ,]
  input$idx <- round(input$end/bin_size)*bin_size
  input_a <- input %>% group_by(chr, idx) %>% summarize(mC=sum(mC), aC=sum(aC))
  input_a$rC <- input_a$mC/input_a$aC
  dat_mer <- merge(dat_mer, input_a[,c("chr", "idx","rC")], by=c("chr", "idx"))
}

nrow(dat_mer)
write.table(dat_mer, "output/sccaps_30MB.table", quote = FALSE, col.names = F, row.names = F)
#dat_mer <- read.table("output/sccaps_30MB.table")
#dat_mer$uni_id <- seq.int(nrow(dat_mer))
dat_mer$uni_id <- dat_mer$V2/bin_size
dat_mer$uni_id <- factor(dat_mer$uni_id, level=seq.int(405,555))
dat_mer <- dat_mer[,-3] %>% dplyr::select(-V1, -V2)
dat_mer <- dat_mer[order(dat_mer$uni_id,decreasing=FALSE),]

final_matrix <- data.frame(dat_mer[,-97], row.names=as.character(dat_mer[,97])) %>% as.matrix(rownames = TRUE)
library(pheatmap)
library(RColorBrewer)
pdf("plot/sccaps_heatmap_30mb.pdf",width=24, height=12)
pheatmap(t(final_matrix),cluster_cols = FALSE, cluster_rows = FALSE,border_color=NA,show_rownames=FALSE, show_colnames=FALSE,color=colorRampPalette(c("navy", "white", "red"))(50))
dev.off()





###Fig1f
#In Linux
```
computeMatrix scale-regions -S meth_old/10k_cells_caps.md.filter.meth.sta.bed.bw meth_update/merged_sccaps.bed.bw -R /gpfs3/well/ludwig/users/ebu571/resource/mm9.refGene.gtf.gz \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 5000 \
                             -o meth_update/merged_sccaps.updown.5k.10bin.mat.gz --binSize 10 -p 8

        
plotProfile -m meth_update/merged_sccaps.updown.5k.10bin.mat.gz \
              -out plot/sc_caps.updown10k.5k.pdf \
              --perGroup --legendLocation "center-right" --yMax 0.20 --colors "#2D6A4E" "#c7e1a6"

computeMatrix scale-regions -S /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/meth/cd8_tcells_merge_CpG.ratio.bw /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/meth/C183_bulk_CpG.ratio.bw -R /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/resource/gencode.v43.basic.annotation.gtf.gz \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 5000 \
                             -o meth_update/merged_sctaps.tcell.updown.5k.10bin.mat.gz --binSize 10 -p 8

plotProfile -m meth_update/merged_sctaps.tcell.updown.5k.10bin.mat.gz \
              -out plot/sc_taps.tcell.updown10k.5k.pdf \
             --perGroup --legendLocation "center-right" --yMax 100 --colors "#0066CC" "#A3DEF4"
```


###FigS4
#####plot correlation 
##cp /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/stats/per_cell_100k.cor.mat /users/ludwig/ebu571/ebu571/project/single_cell/scTAPS_CAPS_RNA/output
cor_caps <- read.table("output/each.sccaps.correlation.table")
cor_caps$smp <- "scCAPS"
colnames(cor_caps) <- c("smp2","smp1","cor", "smp")
cor_taps <- read.table("output/per_cell_100k.cor.mat")
cor_taps$smp <- "scTAPS"
colnames(cor_taps) <- c("smp1", "smp2","cor","smp")

##exclude bulk and merged corelation for scTAPS, for scCAPS, bulk is 10k cell sample
###revised heatmap

cor_scCAPS_scTAPS <- rbind(cor_caps[grepl("10k",cor_caps$smp2),c("cor", "smp")],cor_taps[cor_taps$smp2!="C183_bulk_CpG.bin_100k.bed" & cor_taps$smp2!="cd8_tcells_merge_CpG.bin_100k.bed",c("cor", "smp")] )
cor_scCAPS_scTAPS$smp <- factor(cor_scCAPS_scTAPS$smp, level=c("scTAPS", "scCAPS"))
p <- ggplot(cor_scCAPS_scTAPS,aes(x=smp,y=cor, fill=smp))+ geom_boxplot(outlier.shape=NA,width=0.5) +
  geom_jitter(size=0.5, width = 0.2) +
  stat_summary(geom="text", fun=mean,
               aes(label=sprintf("%1.3f", ..y..)),
              position=position_nudge(y=0.1), size=10)+
  theme_bw() +
  ylab("correlation") +
  xlab("methods")+
  theme(legend.position="None") +scale_fill_manual(values=c("#A3DEF4","#c7e1a6"))+scale_y_continuous(limits=c(0,1),expand = c(0, 0))
ggsave("plot/cor_sctaps_sccaps.pdf", p,width=4, height=4)


####cor value and raw read

cor_scCAPS_scTAPS <- rbind(cor_caps[grepl("10k",cor_caps$smp2),c("cor", "smp", "smp1")],cor_taps[cor_taps$smp1!="C183_bulk_CpG.bin_100k.bed" & cor_taps$smp1!="cd8_tcells_merge_CpG.bin_100k.bed",c("cor", "smp", "smp1")] )

cor_scCAPS_scTAPS$smp1 <- gsub("_rev.CpG.bin_100k.bed","_rev",cor_scCAPS_scTAPS$smp1)
cor_scCAPS_scTAPS_com <- merge(cor_scCAPS_scTAPS,taps_caps, by.x=("smp1"), by.y=c("smp") )
cor_scCAPS_scTAPS_com$smp <- factor(cor_scCAPS_scTAPS_com$smp, level=c("scTAPS", "scCAPS"))

p_cor_depth <- ggplot(cor_scCAPS_scTAPS_com, aes(x=as.numeric(nclean)/1000000, y=cor,fill=smp)) +
geom_point(aes(size=0.5, fill=smp),color="black",shape=21)+scale_fill_manual(values=c("#A3DEF4","#c7e1a6"))+
  #geom_smooth(aes(fill=lib, color=lib) ,method = "lm", se = FALSE)+theme_bw()+
  scale_y_continuous(limits=c(0,1),expand = c(0, 0))+ylab("correlation")+theme_bw()+
    scale_x_continuous(limits=c(0,20),expand = c(0, 0))+ylab("correlation")
  
  
ggsave("plot/correlation_raw.pdf", p_cor_depth, width=6, height=4)

