#### get all statistics ####
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

#### get methylation for region ####
module purge
module load BEDTools
for i in `ls meth/*.meth.bed.gz`
do
bedtools intersect  -a <(cut -f1-3 resource/gencode.vM1.annotation.protein_coding.bed |sort -k1,1 -k2,2n ) \
    -b <(zcat $i|tail -n +2) -wa -wb -sorted |\
    cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c 5,6|\
    cat <(zcat $i|head -1|cut -f1-3,5,6) - >${i/.meth.bed.gz/}.gene.bed

bedtools intersect  -a <(cut -f1-3 resource/genome_bin100k.bed |sort -k1,1 -k2,2n ) \
    -b <(zcat $i|tail -n +2) -wa -wb -sorted |\
    cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c 5,6|\
    cat <(zcat $i|head -1|cut -f1-3,5,6) - >${i/.meth.bed.gz/}.genome_bin100k.bed
bedtools intersect  -a <(cut -f1-3 resource/gencode.vM1.annotation.protein_coding_slope2k.bed |sort -k1,1 -k2,2n ) \
    -b <(zcat $i|tail -n +2) -wa -wb -sorted |\
    cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c 5,6|\
    cat <(zcat $i|head -1|cut -f1-3,5,6) - >${i/.meth.bed.gz/}.gene_2k.bed
done

meth_dir=/users/ludwig/cfo155/cfo155/scTAPS_CAPS/mouse_neuron_update/meth
eval $(ls $meth_dir/C*gene.bed|tr '\n' ',' |sed 's/,$/|cut -f4-/g'|sed 's/^/paste <(cat /g;s/,/|cut -f4- ) <(cat /g;s/$/)/g')|\
    paste <(cat `ls $meth_dir/C*gene.bed|head -1 `|cut -f1-3) -  >all_sample.CpG.gene.txt
eval $(ls $meth_dir/C*genome_bin100k.bed|tr '\n' ',' |sed 's/,$/|cut -f4-/g'|sed 's/^/paste <(cat /g;s/,/|cut -f4- ) <(cat /g;s/$/)/g')|\
    paste <(cat `ls $meth_dir/C*genome_bin100k.bed|head -1 `|cut -f1-3) -  >all_sample.CpG.genome_bin100k.txt


#### call CH methylation ####
module purge
module load BEDTools
module load SAMtools/1.18-GCC-12.3.0
for i in C204 C206 C211 C213 
do
samtools merge -o ${i}.merge.bam ${i}*md.bam
samtools sort -o ${i}.merge.sort.bam ${i}.merge.bam 
done

for i in C204 C206 C211 C213 
do
samtools index ${i}.merge.sort.bam
/gpfs3/well/ludwig/users/cfo155/miniconda3/bin/MethylDackel extract \
    /gpfs3/well/ludwig/users/ebu571/14Jul2022_2022_100cycles/resource/caps_mm9_lambda.fa \
    ${i}.merge.sort.bam -o ${i}.merge -p 13 \
    --OT 10,118,10,118 --OB 10,118,10,118 --cytosine_report
done


#### get cpg number in selected genes ####
grep -P "MXD4|CD47|VEGFA|ALKBH5|TMED7|EXOC4|ARHGAP15"  gencode.v43.annotation.protein_coding.bed|\
    sort -k1,1 -k2,2n|\
    bedtools intersect -a - -b <(sort -k1,1 -k2,2n GRCh38_no_alt.cg.pos.excluded.bed) -wa -wb -sorted |\
    bedtools groupby -i - -grp 1,2,3,4,5,6 -c 7 -o count >selgene_cpg_count.txt