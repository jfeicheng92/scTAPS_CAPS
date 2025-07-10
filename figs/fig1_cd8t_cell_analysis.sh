########### merge single-cell to bulk ########### 
module purge
module load BEDTools
meth_dir=/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/align

eval $(ls $meth_dir/C183_N*md.filter.meth.sta.txt.gz|tr '\n' ',' |sed 's/,$/|tail -n +2|cut -f1-4/g'|sed 's/^/cat <(zcat /g;s/,/|tail -n +2|cut -f1-4 ) <(zcat /g;s/$/)/g')|\
    awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$2+2,$3,$4}' |bedtools sort -i - |\
    bedtools groupby -i - -g 1,2,3 -c 4,5 -o sum,sum  >C183_bulk.md.filter.meth.sta.txt

zcat meth/C183_bulk.md.filter.meth.sta.txt.gz|grep ^chr|sort -k1,1 -k2,2n|\
    bedtools intersect  -a <(cut -f1-3 resource/genome_bin100k.bed |sort -k1,1 -k2,2n ) -b - -wa -wb -sorted |\
    cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c 4,5|\
    cat <(echo -e "chr\tstart\tend\tmC\ttotal") - >meth/C183_bulk_CpG.bin_100k.bed
zcat meth/C183_bulk.md.filter.meth.sta.txt.gz|grep ^chr|awk 'BEGIN{FS="\t";OFS="\t"}{if($5>0)printf("%s\t%d\t%d\t%.2f\n", $1,$2,$3,$4/$5*100)}'|\
    cut -f1-4|sort -k1,1 -k2,2n >meth/C183_bulk_CpG.ratio.bedGraph

########### create files for standard taps ########### 
module purge
module load BEDTools
zcat meth/cd8_tcells_merge_CpG.bedGraph.gz|grep ^chr|sort -k1,1 -k2,2n|\
    bedtools intersect  -a <(cut -f1-3 resource/genome_bin100k.bed |sort -k1,1 -k2,2n ) -b - -wa -wb -sorted |\
    cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c 5,6|\
    cat <(echo -e "chr\tstart\tend\tmC_standard\ttotal_standard") - >meth/cd8_tcells_merge_CpG.bin_100k.bed
zcat meth/cd8_tcells_merge_CpG.bedGraph.gz|grep ^chr|cut -f1-4|sort -k1,1 -k2,2n >meth/cd8_tcells_merge_CpG.ratio.bedGraph


for i in `ls meth/*ratio.bedGraph|grep C183`
do
/gpfs3/well/ludwig/users/cfo155/tools/ucsc_utils/bedGraphToBigWig $i resource/GRCh38_spike_ins.fasta.fai ${i/.bedGraph}.bw
done
zcat meth/C183_bulk.md.filter.meth.sta.txt.gz|grep ^chr|sort -k1,1 -k2,2n|cut -f1,2,3,5 >meth/C183_bulk_CpG.totalC.bed
zcat meth/cd8_tcells_merge_CpG.bedGraph.gz|grep ^chr|sort -k1,1 -k2,2n|cut -f1,2,3,6 >meth/cd8_tcells_merge_CpG.totalC.bed
########### get stats for all single-cell sample ########### 
for i in `ls align/C*.md.filter.meth.sta.txt.gz`
do
smp=`basename $i |sed 's/.md.filter.meth.sta.txt.gz//g'`
echo -e "$smp,"\
`cat fastq/${smp}_1_fastqc.html|sed 's/<[^>]*>/\n/g'|grep -i total -A2|sed -n 3p`","\
`cat fastq/${smp}_fastp_1_fastqc.html|sed 's/<[^>]*>/\n/g'|grep -i total -A2|sed -n 3p`","\
`cat stats/${smp}.mapping.txt|tail -n +2|cut -f2-|sed 's/\t/,/g'`","\
`zcat align/${smp}.spikeins.filter.meth.sta.txt.gz |grep J02459.1 |awk '$4>0'|awk 'BEGIN{OFS=",";FS="\t"}{mC+=$3;aC+=$4;nC+=1}END{print mC/aC*100,mC,aC,nC}'`","\
`zcat align/${smp}.spikeins.filter.meth.sta.txt.gz |grep unmodified_2kb |awk '$4>0'|awk 'BEGIN{OFS=",";FS="\t"}{mC+=$3;aC+=$4;nC+=1}END{print mC/aC*100,mC,aC,nC}'`","\
`zcat align/${smp}.spikeins.filter.meth.sta.txt.gz |grep 144hmC |awk '$4>0&&($2=="86"||$2=="91"||$2=="99"||$2=="109")'|awk 'BEGIN{OFS=",";FS="\t"}{mC+=$3;aC+=$4;nC+=1}END{print mC/aC*100,mC,aC,nC}'`","\
`zcat align/${smp}.md.filter.meth.sta.txt.gz|grep ^chr[0-9,X,Y] |awk '$4>0'|awk 'BEGIN{OFS=",";FS="\t"}{mC+=$3;aC+=$4;nC+=1}END{print mC/aC*100,mC,aC,nC}'`
done |cat <(echo "smp,nraw,nclean,q1_nmap,q10_nmap,proper_nmap,mean_isize,dup_rate,lambda_rC,lambda_mC,lambda_aC,lambda_nC,unmeth2kb_rC,unmeth2kb_mC,unmeth2kb_aC,unmeth2kb_nC,hmC_rC,hmC_mC,hmC_aC,hmC_nC,chr_rC,chr_mC,chr_aC,chr_nC") - >stats/all_stats.info



for i in `ls align/C*.md.filter.meth.sta.txt.gz`
do
smp=`basename $i |sed 's/.md.filter.meth.sta.txt.gz//g'`
echo -e "$smp,"\
`cat fastq/${smp}_1_fastqc.html|sed 's/<[^>]*>/\n/g'|grep -i total -A2|sed -n 3p`","\
`cat fastq/${smp}_fastp_1_fastqc.html|sed 's/<[^>]*>/\n/g'|grep -i total -A2|sed -n 3p`","\
`cat stats/${smp}.mapping.txt|tail -n +2|cut -f2-|sed 's/\t/,/g'`","\
`zcat align/${smp}.spikeins.filter.meth.sta.txt.gz |grep J02459.1 |awk '$4>0'|awk 'BEGIN{OFS=",";FS="\t"}{mC+=$3;aC+=$4;nC+=1}END{print mC/aC*100,mC,aC,nC}'`","\
`zcat align/${smp}.spikeins.filter.meth.sta.txt.gz |grep unmodified_2kb |awk '$4>0'|awk 'BEGIN{OFS=",";FS="\t"}{mC+=$3;aC+=$4;nC+=1}END{print mC/aC*100,mC,aC,nC}'`","\
`zcat align/${smp}.spikeins.filter.meth.sta.txt.gz |grep 144hmC |awk '$4>0&&($2=="86"||$2=="91"||$2=="99"||$2=="109")'|awk 'BEGIN{OFS=",";FS="\t"}{mC+=$3;aC+=$4;nC+=1}END{print mC/aC*100,mC,aC,nC}'`","\
`zcat meth/${smp}.CpG.meth.bed.gz|grep ^chr[0-9] |awk '$6>0'|awk 'BEGIN{OFS=",";FS="\t"}{mC+=$5;aC+=$6;nC+=1}END{print mC/aC*100,mC,aC,nC}'`
done |cat <(echo "smp,nraw,nclean,q1_nmap,q10_nmap,proper_nmap,mean_isize,dup_rate,lambda_rC,lambda_mC,lambda_aC,lambda_nC,unmeth2kb_rC,unmeth2kb_mC,unmeth2kb_aC,unmeth2kb_nC,hmC_rC,hmC_mC,hmC_aC,hmC_nC,chr_rC,chr_mC,chr_aC,chr_nC") - >stats/all_stats.new.info

########### methylation around gene & coverage around TSS ########### 
/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/meth/C183_bulk_CpG.ratio.bw
/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/meth/cd8_tcells_merge_CpG.ratio.bw
/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/meth/C183_bulk_CpG.totalC.bed
/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/meth/cd8_tcells_merge_CpG.totalC.bed
/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/resource/gencode.v43.basic.annotation.gtf.gz
/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/human_immune_cells/t_cell_update/resource/GRCh38_spike_ins.fasta.fai
/gpfs3/well/ludwig/users/cfo155/cfDNA/072019_cfDNA/resource/cpgIsland.bed


########### methylation within chromhmm ########### 
for i in `ls *bed`
do
awk -v j=${i/.bed/} 'BEGIN{OFS="\t"}{print $1,$2,$3,j}' $i
done |sort -k1,1 -k2,2n >all_state.txt

for i in `ls meth/*chrhmm.meth.txt`
do
awk -v j=${i/_rev.CpG.chrhmm.meth.txt/} 'BEGIN{OFS="\t"}{print j,$0}' $i
done|cut -d '/' -f2 >stats/all_state.meth.txt
