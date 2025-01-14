##This code is to make figure2 to demonstrate the performance of low-input TAPS and CAPS on mESC
##module load R/4.0.3-foss-2020b

#.libPaths("/well/ludwig/users/ebu571/R/4.0/skylake")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
library(readr)
setwd("/users/ludwig/ebu571/ebu571/project/single_cell/scTAPS_CAPS_RNA")

###Fig.S1b
smps <- c("100_cells_caps", "100_cells_taps", "10_cells_caps", "10_cells_taps", "10k_cells_caps", "10k_cells_taps", "1k_cells_caps", "1k_cells_taps")
res2 <- NULL
for(smp in smps){
  input_file <- paste0("/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/align/",smp,".spikeins.filter.meth.sta.txt.gz")
  dat <- suppressMessages(fread(paste0('zcat ', input_file), header=TRUE))
  mC <- sum(dat[dat$chr=="J02459.1",3])/sum(dat[dat$chr=="J02459.1",4])*100 %>% round(2)
  hmC <- sum(dat[dat$chr=="144hmC" & dat$pos %in% c(86, 91, 99, 109),3])/sum(dat[dat$chr=="144hmC" & dat$pos %in% c(86, 91, 99, 109),4]) *100 %>% round(2)
  uC <- sum(dat[dat$chr=="unmodified_2kb",3])/sum(dat[dat$chr=="unmodified_2kb",4])*100 %>% round(2)
  res2 <- rbind(res2, c(smp, mC, hmC, uC))
}
res2 <- as.data.frame(res2)
colnames(res2) <- c("smp","mC","hmC","uC")
res2 <- melt(res2,id.vars="smp") 
res2$value <- as.character(res2$value) %>% as.numeric()
res2$lib <- gsub(".*cells_","",res2$smp)%>% toupper()
res2$smp <- gsub("_taps|_caps","",res2$smp)
res2$value <- round(res2$value,4)
res2$smp <- factor(res2$smp, levels=c("10_cells","100_cells","1k_cells","10k_cells"))
res2$variable <- factor(res2$variable, levels=c("hmC","mC","uC"))
p2 <- ggplot(res2,aes(x=smp,y=value,fill=lib))+geom_bar(width=0.8,stat="identity", position = "dodge2") + 
  facet_wrap(lib~ variable) +  theme_bw() +theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
)+
  geom_text(aes(label=value), position = "dodge2", vjust=-0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "NA") +
  ylim(0,100)+
  ylab("modification level on spike-ins%") +
  scale_fill_manual(values=c("#1B9E77", "#D95F02")) 

pdf("plot/lowinput_qc.pdf",width=12, height = 5)
cowplot::plot_grid(p1, p2, ncol = 2,align = "h", rel_widths = c(1,2.6))
dev.off()


###Fig.S1c
window_size <- 100000
cov_min <- 1000
cov_max <- 10000

###CpG coverage
cov_dat <- read.table("output/genome_coverage.txt")
colnames(cov_dat) <- c("sample", "CpG_num", "per_base")
seldat <- cov_dat[grepl("cells_caps",cov_dat$sample) | grepl("cells_taps",cov_dat$sample),]
seldat[grepl("taps",seldat$sample),c("group")]<- "taps"
seldat[grepl("caps",seldat$sample),c("group")] <- "caps"
seldat[grepl("10",seldat$sample),c("cell_num")]<- "10_cells" 
seldat[grepl("100",seldat$sample),c("cell_num")]<- "100_cells" 
seldat[grepl("1k",seldat$sample),c("cell_num")]<- "1k_cells" 
seldat[grepl("10k",seldat$sample),c("cell_num")]<- "10k_cells" 
seldat$cell_num <- factor(seldat$cell_num, levels=c("10_cells","100_cells","1k_cells","10k_cells"))

p1 <- ggplot(seldat,aes(x=cell_num,y=round(CpG_num/21342780,3)*100,fill=group))+geom_bar(width=0.8,stat="identity", position = "dodge2") + 
  facet_wrap(~ group) +  theme_bw() +theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
)+
  geom_text(aes(label=round(CpG_num/21342780,3)*100), position = "dodge2", vjust=-0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "NA") +
  scale_y_continuous(limits=c(0,100),expand = c(0, 0))+
  ylab("genomic CpG coverage") +scale_fill_manual(values=c("#1B9E77", "#D95F02")) 


p2 <- ggplot(seldat,aes(x=cell_num,y=per_base,fill=group))+geom_bar(width=0.8,stat="identity", position = "dodge2") + 
  facet_wrap(~ group) +  theme_bw() +theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
)+
  geom_text(aes(label=round(per_base,1)), position = "dodge2", vjust=-0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "NA") +
  scale_y_continuous(limits=c(0,8),expand = c(0, 0))+
  ylab("per-base CpG depth") +scale_fill_manual(values=c("#1B9E77", "#D95F02")) 

pdf("plot/lowinput_qc_coverage1.pdf",width=4, height = 6)
cowplot::plot_grid(p1, p2, ncol = 1,align = "h")
dev.off()



###Fig.S1f

######plot correlation

window_size <- 100000
cov_min <- 100
cov_max <- 3000


cor_caps <- NULL
smps <- c("10_cells","100_cells", "1k_cells", "10k_cells")
for(smp in smps){
  library_type <- "caps"
  input_file1 <- paste0("meth_old/",smp,"_",library_type,".md.filter.meth.sta.txt.gz")
  lowinput <- fread(paste0('gzip -cd ', input_file1), header=TRUE)
  lowinput <- lowinput[lowinput$chr%in%paste0("chr",seq(1,19)),]
  lowinput$idx <- round(lowinput$pos/window_size)*window_size
  lowinput_a <- lowinput %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC)) 
  lowinput_a$rC <- lowinput_a$mC/lowinput_a$aC
  sum(lowinput$mC)/sum(lowinput$aC)
  
  input_file2 <- paste0("meth_old/standard_",library_type,"_merge_CpG.bedGraph.gz")
  standard <- fread(paste0('gzip -cd ', input_file2), header=TRUE)
  colnames(standard)[1:2] <- c("chr","pos")
  standard <- standard[standard$chr%in%paste0("chr",seq(1,19)),]
  standard$idx <- round(standard$pos/window_size)*window_size
  standard_a <- standard %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(total))
  standard_a$rC <- standard_a$mC/standard_a$aC
  sum(standard$mC)/sum(standard$total)
  
  dat <- merge(lowinput_a, standard_a, by=c("chr","idx"))
  dat <- dat[complete.cases(dat) & dat$aC.x<cov_max&dat$aC.x>cov_min,]
  cor_caps <- rbind(cor_caps, c("sd", smp,round(cor(dat$rC.x, dat$rC.y),2)))

}

for(smp in smps){
  for (j in smps){
  library_type <- "caps"
  input_file1 <- paste0("meth/",smp,"_",library_type,".md.filter.meth.sta.txt.gz")
  lowinput <- fread(paste0('gzip -cd ', input_file1), header=TRUE)
  lowinput <- lowinput[lowinput$chr%in%paste0("chr",seq(1,19)),]
  lowinput$idx <- round(lowinput$pos/window_size)*window_size
  lowinput_a <- lowinput %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC)) 
  lowinput_a$rC <- lowinput_a$mC/lowinput_a$aC
  sum(lowinput$mC)/sum(lowinput$aC)
  
  input_file2 <- paste0("meth/",j,"_",library_type,".md.filter.meth.sta.txt.gz")
  lowinput2 <- fread(paste0('gzip -cd ', input_file2), header=TRUE)
  lowinput2 <- lowinput2[lowinput2$chr%in%paste0("chr",seq(1,19)),]
  lowinput2$idx <- round(lowinput2$pos/window_size)*window_size
  lowinput2_a <- lowinput2 %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC)) 
  lowinput2_a$rC <- lowinput2_a$mC/lowinput2_a$aC
  sum(lowinput2$mC)/sum(lowinput2$aC)
  
  dat <- merge(lowinput_a, lowinput2_a, by=c("chr","idx"))
  dat <- dat[complete.cases(dat) & dat$aC.x<cov_max&dat$aC.x>cov_min,]
  cor_caps <- rbind(cor_caps, c(smp, j,round(cor(dat$rC.x, dat$rC.y),2)))

}
}

cor_caps <- as.data.frame(cor_caps)
colnames(cor_caps) <- c("smp1", "smp2", "cor_value")
cor_caps <- rbind(cor_caps, c("10_cells", "sd", 0.90))
cor_caps <- rbind(cor_caps, c("100_cells", "sd", 0.93))
cor_caps <- rbind(cor_caps, c("1k_cells", "sd", 0.95))
cor_caps <- rbind(cor_caps, c("10k_cells", "sd", 0.98))
cor_caps <- rbind(cor_caps, c("sd", "sd", 1))
cor_caps$smp1 <- factor(cor_caps$smp1, level=c("sd", "10k_cells", "1k_cells", "100_cells", "10_cells"))
cor_caps$smp2 <- factor(cor_caps$smp2, level=c("sd", "10k_cells", "1k_cells", "100_cells", "10_cells"))
cor_caps$cor_value <- as.numeric(cor_caps$cor_value)

get_upper_tri <- function(CorMat){
  CorMat[upper.tri(CorMat)]<- NA
  return(CorMat)
}

library(corrplot)
library(tidyr)
cot_value <- spread(cor_caps,smp2,cor_value) %>% as.data.frame()
cot_value <- data.frame(cot_value[,-1], row.names=as.character(cot_value[,1])) %>% as.matrix(rownames = TRUE,colnames=TRUE)
cor_caps1 <- get_upper_tri(cot_value)
meltNum <- melt(cor_caps1, na.rm = T)

pdf("plot/tpm_correlation_lowinput_caps.pdf")
ggplot(meltNum, aes(x = Var1, y = Var2, fill = value)) + geom_tile()+coord_fixed()+
  geom_text(aes(label = round(value,2)), color = "black")+  scale_fill_gradient(low = "white", high = "firebrick4",
                      limit = c(0.8,1), name = "Pearson\nCorrelation") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
dev.off()




#####For low input taps
window_size <- 100000
cov_min <- 1000
cov_max <- 10000


cor_taps <- NULL
smps <- c("10_cells","100_cells", "1k_cells", "10k_cells")
for(smp in smps){
  library_type <- "taps"
  input_file1 <- paste0("meth/",smp,"_",library_type,".md.filter.meth.sta.txt.gz")
  lowinput <- fread(paste0('gzip -cd ', input_file1), header=TRUE)
  lowinput <- lowinput[lowinput$chr%in%paste0("chr",seq(1,19)),]
  lowinput$idx <- round(lowinput$pos/window_size)*window_size
  lowinput_a <- lowinput %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC)) 
  lowinput_a$rC <- lowinput_a$mC/lowinput_a$aC
  sum(lowinput$mC)/sum(lowinput$aC)
  
  input_file2 <- paste0("meth/standard_",library_type,"_merge_CpG.bedGraph.gz")
  standard <- fread(paste0('gzip -cd ', input_file2), header=TRUE)
  colnames(standard)[1:2] <- c("chr","pos")
  standard <- standard[standard$chr%in%paste0("chr",seq(1,19)),]
  standard$idx <- round(standard$pos/window_size)*window_size
  standard_a <- standard %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(total))
  standard_a$rC <- standard_a$mC/standard_a$aC
  sum(standard$mC)/sum(standard$total)
  
  dat <- merge(lowinput_a, standard_a, by=c("chr","idx"))
  dat <- dat[complete.cases(dat) & dat$aC.x<cov_max&dat$aC.x>cov_min,]
  cor_taps <- rbind(cor_taps, c("sd", smp,round(cor(dat$rC.x, dat$rC.y),2)))

}

for(smp in smps){
  for (j in smps){
  library_type <- "taps"
  input_file1 <- paste0("meth/",smp,"_",library_type,".md.filter.meth.sta.txt.gz")
  lowinput <- fread(paste0('gzip -cd ', input_file1), header=TRUE)
  lowinput <- lowinput[lowinput$chr%in%paste0("chr",seq(1,19)),]
  lowinput$idx <- round(lowinput$pos/window_size)*window_size
  lowinput_a <- lowinput %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC)) 
  lowinput_a$rC <- lowinput_a$mC/lowinput_a$aC
  sum(lowinput$mC)/sum(lowinput$aC)
  
  input_file2 <- paste0("meth/",j,"_",library_type,".md.filter.meth.sta.txt.gz")
  lowinput2 <- fread(paste0('gzip -cd ', input_file2), header=TRUE)
  lowinput2 <- lowinput2[lowinput2$chr%in%paste0("chr",seq(1,19)),]
  lowinput2$idx <- round(lowinput2$pos/window_size)*window_size
  lowinput2_a <- lowinput2 %>% 
    group_by(chr, idx) %>% 
    summarize(mC = sum(mC), aC=sum(aC)) 
  lowinput2_a$rC <- lowinput2_a$mC/lowinput2_a$aC
  sum(lowinput2$mC)/sum(lowinput2$aC)
  
  dat <- merge(lowinput_a, lowinput2_a, by=c("chr","idx"))
  dat <- dat[complete.cases(dat) & dat$aC.x<cov_max&dat$aC.x>cov_min,]
  cor_taps <- rbind(cor_taps, c(smp, j,round(cor(dat$rC.x, dat$rC.y),2)))

}
}

cor_taps <- as.data.frame(cor_taps)
colnames(cor_taps) <- c("smp1", "smp2", "cor_value")
cor_taps <- rbind(cor_taps, c("10_cells", "sd", 0.92))
cor_taps <- rbind(cor_taps, c("100_cells", "sd", 0.92))
cor_taps <- rbind(cor_taps, c("1k_cells", "sd", 0.93))
cor_taps <- rbind(cor_taps, c("10k_cells", "sd", 0.93))
cor_taps <- rbind(cor_taps, c("sd", "sd", 1))
cor_taps$smp1 <- factor(cor_taps$smp1, level=c("sd", "10k_cells", "1k_cells", "100_cells", "10_cells"))
cor_taps$smp2 <- factor(cor_taps$smp2, level=c("sd", "10k_cells", "1k_cells", "100_cells", "10_cells"))
cor_taps$cor_value <- as.numeric(cor_taps$cor_value)

cot_value <- spread(cor_taps,smp2,cor_value) %>% as.data.frame()
cot_value <- data.frame(cot_value[,-1], row.names=as.character(cot_value[,1])) %>% as.matrix(rownames = TRUE,colnames=TRUE)
cor_taps1 <- get_upper_tri(cot_value)
cor_taps1 <- melt(cor_taps1, na.rm = T)

pdf("plot/tpm_correlation_lowinput_taps.pdf")
ggplot(cor_taps1, aes(x = Var1, y = Var2, fill = value)) + geom_tile()+coord_fixed()+
  geom_text(aes(label = round(value,2)), color = "black")+  scale_fill_gradient(low = "white", high = "firebrick4",
                      limit = c(0.8,1), name = "Pearson\nCorrelation") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
dev.off()

###Fig.S1g



``` ##IN LINUX
computeMatrix scale-regions -S meth_old/10k_cells_taps.md.filter.meth.sta.bed.bw meth_old/1k_cells_taps.md.filter.meth.sta.bed.bw meth_old/100_cells_taps.md.filter.meth.sta.bed.bw meth_old/10_cells_taps.md.filter.meth.sta.bed.bw -R /gpfs3/well/ludwig/users/ebu571/resource/mm9.refGene.gtf.gz \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 5000 \
                              -o meth_old/lowinput_taps.updown10k.5k.mat.gz --binSize 10 -p 8

plotProfile -m meth_old/lowinput_taps.updown10k.5k.mat.gz \
              -out plot/lowinput_taps.updown10k.5k.pdf \
              --perGroup --legendLocation "center-right" --yMax 1 --colors "#3B75A1" "#698FB4" "#5EA1C6" "#86B1D0" 

computeMatrix scale-regions -S meth_old/10k_cells_caps.md.filter.meth.sta.bed.bw meth_old/1k_cells_caps.md.filter.meth.sta.bed.bw meth_old/100_cells_caps.md.filter.meth.sta.bed.bw meth_old/10_cells_caps.md.filter.meth.sta.bed.bw -R /gpfs3/well/ludwig/users/ebu571/resource/mm9.refGene.gtf.gz \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 5000 \
                             -o meth_old/lowinput_caps.updown.5k.mat.gz --binSize 10 -p 8
        
plotProfile -m meth_old/lowinput_caps.updown.5k.mat.gz \
              -out plot/lowinput_caps.updown.5k.pdf \
              --perGroup --legendLocation "center-right" --yMax 0.25 --colors "#40916C" "#52B788" "#74C69D" "#95D5B2"              
              
              ```
