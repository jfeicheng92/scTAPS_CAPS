##This code is to make figure2 to demonstrate the performance of low-input TAPS and CAPS on mESC
#ml use -a /apps/eb/2020b/skylake/modules/all
##module load R/4.0.3-foss-2020b

#.libPaths("/well/ludwig/users/ebu571/R/4.0/skylake")
library(data.table)
library(parallel)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library("cowplot")
setwd("/users/ludwig/ebu571/ebu571/project/single_cell/scTAPS_CAPS_RNA")

###comparision between different methods scTAPS
sctaps_qc <- read.table("output/scTAPS_Tcell",header=TRUE)
sctaps_qc$lib <- "scTAPS"
sctaps_qc$perc_cpg <- as.numeric(sctaps_qc$chr_nC)/24314257

sccaps_qc <- read.table("output/scCAPS_mESC",header=TRUE)
sccaps_qc$lib <- "scCAPS"
sccaps_qc$perc_cpg <- as.numeric(sccaps_qc$chr_nC)/21342780

taps_caps <- rbind(sctaps_qc, sccaps_qc)
##number of CpG in human genome: 24314257, after filter blacklist
##number of CpG in mice: 21342780

others <- read.table("output/other_methods", header=TRUE)

p1 <- ggplot(taps_caps, aes(x=as.numeric(nclean)/1000000, y=perc_cpg)) +
#geom_point(aes(size=1, color=lib))+
  geom_smooth(aes(fill=lib, color=lib) ,method = "lm", se = FALSE)+theme_bw()+
  scale_y_continuous(limits=c(0,0.301),expand = c(0, 0))+ylab("% of CpG covered")+
  scale_x_continuous(limits=c(0,25),expand = c(0, 0))+xlab("Reads/cell(Million)")+
  geom_point(data=others, aes(x=as.numeric(raw)/1000000, y=as.numeric(perc_cpg)),size=0.8, color="black", fill="grey")+
  geom_text(data=others, aes(x=as.numeric(raw)/1000000, y=as.numeric(perc_cpg),label=methods), position = position_nudge(y = -0.01,x=0.01))

p2 <- ggplot(taps_caps, aes(x=as.numeric(nclean)/1000000, y=perc_chr)) +
#geom_point(aes(size=1, color=lib))+
  geom_smooth(aes(fill=lib, color=lib) ,method = "lm", se = FALSE)+theme_bw()+
  scale_y_continuous(limits=c(0,0.301),expand = c(0, 0))+ylab("% of genome covered")+
  scale_x_continuous(limits=c(0,25),expand = c(0, 0))+xlab("Reads/cell(Million)")+
  geom_point(data=others, aes(x=as.numeric(raw)/1000000, y=as.numeric(genomic_cov)),size=0.8, color="black", fill="grey")+
  geom_text(data=others, aes(x=as.numeric(raw)/1000000, y=as.numeric(genomic_cov),label=methods), position = position_nudge(y = -0.01,x=0.01))
#ggsave("plot/sc_mC_comparision.pdf", p2)
pdf("plot/sc_mC_comparision.pdf", width=6,height=6)
plot_grid(p1, p2, labels = c('A', 'B'),ncol=1)
dev.off()
######mapping rate
mapping_sta <- others[,c("methods", "mapping_rate")]
mapping_sta <- rbind(mapping_sta, c("scTAPS", 0.930))
mapping_sta <- rbind(mapping_sta, c("scCAPS+", 0.894))
mapping_sta[mapping_sta$methods=="scTAPS",c("group")] <- c("scTAPS")
mapping_sta[mapping_sta$methods=="scCAPS+",c("group")] <- c("scCAPS+")
mapping_sta[mapping_sta$methods!="scTAPS" & mapping_sta$methods!="scCAPS+",c("group")] <- c("others")

mapping_sta$methods <- factor(mapping_sta$methods, level=c("scBS-seq", "snmC-seq2", "scCabernet_E7.5", "scCabernet_k562", "snmC-seq2.split", "snhmC-seq2_split", "simple_seq", "scTAPS", "scCAPS+"))
p_bar <- ggplot(data=mapping_sta, aes(x=methods, y=as.numeric(mapping_rate),fill=group))+
  geom_bar(stat="identity",width=0.8)+
  theme_bw()+scale_fill_manual(values=c("#A0A0A0","#c7e1a6","#A3DEF4"))+scale_y_continuous(limits=c(0,1),expand = c(0, 0))+theme(axis.text.x = element_text(angle = 45))+
  theme(legend.position="None")
ggsave("plot/scTAPS_scCAPS_compared.mapping.pdf", p_bar, width=6,height=4)
######correlate with standard CAPS, pick three cells
window_size <- 100000
cov_min <- 50
cov_max <- 30000
cor_caps <- NULL
smps <- c("meth_old/10k_cells_caps.md.filter.meth.sta.txt.gz","meth_update/merged_sccaps.bed.gz")
sc_smps <- sccaps_qc$smp
#smps <- c("meth_old/standard_caps_merge_CpG.bedGraph.gz")
#j <-  c("meth_update/C203_N702_i7_5_rev_N502_i5_3_rev.md.filter.meth.sta.txt.gz")
for(smp in smps){
  for (j in sc_smps){
  library_type <- "caps"
  if(grepl("merged", smp)){
    input <- read.table(smp,header=F, sep="\t")
    colnames(input) <- c("chr", "start", "end", "mC", "aC", "ratio")
  }else {
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
  

    input2 <- read.table(paste0("meth_update/", j, ".md.filter.meth.sta.txt.gz"),header=T, sep="\t")
    colnames(input2) <- c("chr", "end", "mC", "aC", "ratio")

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
write.table(cor_caps, "output/each.sccaps.correlation.table", quote = FALSE, col.names = F, row.names = F)

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





#####
####Revise heatmap for correlation in Fig1


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


