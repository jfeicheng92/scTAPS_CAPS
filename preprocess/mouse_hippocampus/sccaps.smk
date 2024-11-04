"""
Workflow for single-cell caps

samples include:
    neuron N-
    neuron N+ aged
    neuron N+
    neuron N- aged

read length: 128 bp

rawdata: /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/P230072

run:
module purge
ml use -a /apps/eb/skylake/modules/all 
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake --use-envmodules --max-status-checks-per-second 0.01 --snakefile code/sccaps.smk --cluster "sbatch -p short --mem-per-cpu 30G" -j 384 -np

dir=/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/P230072/WTCHG_981766
paste <(ls ${dir}/*gz) <(ls ${dir}/*gz|sed 's/WTCHG_981766/WTCHG_981767/g') <(ls ${dir}/*gz|awk -F'/' '{print $NF}'|sed 's,^,>fastq/,g;s/[ATCG+]*_R//g')|sed 's/^/cat /g;' >P230072_rename_fastq.sh

# extract CpG
/gpfs3/well/ludwig/users/cfo155/miniconda2/bin/seqkit locate -i -P -p CG /gpfs3/well/ludwig/users/ebu571/14Jul2022_2022_100cycles/resource/caps_mm9_lambda.fa |\
    grep -wP chr[0-9,X,Y]* |awk 'BEGIN{OFS="\t"}{print $1,$5-1,$6-1}' |sort -k1,1 -k2,2n |\
    bedtools intersect -a - -b <(zcat mm9-blacklist.bed.gz|sort -k1,1 -k2,2n) -v -sorted >mm9_cpg.blacklist_removed.test.bed
bedtools intersect -a <(cat mm9_genome.cg.bed |sort -k1,1 -k2,2n) \
    -b <(zcat mm9-blacklist.bed.gz|sort -k1,1 -k2,2n) -v -sorted |cut -f1-3 >mm9_genome.cg.blacklist_removed.bed

# merge meth call 
names=(*CpG.meth.bed)
paste $(for i in `seq 0 $(( ${#names[@]} -1 ))`; do echo ${names[$i]};done|tr '\n' ' ') |\
    cut -f1-3,$(seq  5 6 $(( ${#names[@]} * 6 )) |tr '\n' ',' |sed s/,$//),$(seq  6 6 $(( ${#names[@]} * 6 )) |tr '\n' ',' |sed s/,$//) > all_sample.CpG.bed

meth_dir=~/cfo155/scTAPS_CAPS/mouse_neuron_update/meth
eval $(ls $meth_dir/*meth.bed.gz|tr '\n' ',' |sed 's/,$/|cut -f5/g'|sed 's/^/paste <(zcat /g;s/,/|cut -f5 ) <(zcat /g;s/$/)/g')|\
    tail -n +2|cat <(ls $meth_dir/*meth.bed.gz|xargs -n 1 basename|sed 's/.meth.bed.gz//g'|tr '\n' '\t'|sed 's/\t$//g' ) -|\
    paste <(zcat `ls $meth_dir/*meth.bed.gz|head -1 `|cut -f1-4) - >all_sample.CpG.txt

names=(*CpG.100k.bed)
paste $(for i in `seq 0 $(( ${#names[@]} -1 ))`; do echo ${names[$i]};done|tr '\n' ' ') |\
    cut -f1-3,$(seq  4 5 $(( ${#names[@]} * 5 )) |tr '\n' ',' |sed s/,$//),$(seq  5 5 $(( ${#names[@]} * 5 )) |tr '\n' ',' |sed s/,$//) > all_sample.CpG.100k.bed
names=(*CpG.gene.bed)
paste $(for i in `seq 0 $(( ${#names[@]} -1 ))`; do echo ${names[$i]};done|tr '\n' ' ') |\
    cut -f1-3,$(seq  4 5 $(( ${#names[@]} * 5 )) |tr '\n' ',' |sed s/,$//),$(seq  5 5 $(( ${#names[@]} * 5 )) |tr '\n' ',' |sed s/,$//) > all_sample.CpG.gene.bed
names=(*CpG.gene_2k.bed)
paste $(for i in `seq 0 $(( ${#names[@]} -1 ))`; do echo ${names[$i]};done|tr '\n' ' ') |\
    cut -f1-3,$(seq  4 5 $(( ${#names[@]} * 5 )) |tr '\n' ',' |sed s/,$//),$(seq  5 5 $(( ${#names[@]} * 5 )) |tr '\n' ',' |sed s/,$//) > all_sample.CpG.gene_2k.bed

# generate features
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz
zcat gencode.vM1.annotation.gtf.gz |awk -F "\t" '$3=="gene" && $9~/protein_coding/'|bedtools slop -i - -g caps_mm9_lambda.fa.fai -l 2000 -r 2000 |awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$10,$18,$7}' | tr -d '";'  >gencode.vM1.annotation.protein_coding_slope2k.bed
bedtools makewindows -g <(cut -f1,2 caps_mm9_lambda.fa.fai|grep -wP chr[0-9,X,Y]* ) -w 100000 >genome_bin100k.bed
bedtools makewindows -g <(cut -f1,2 caps_mm9_lambda.fa.fai|grep -wP chr[0-9,X,Y]* ) -w 1000000 >genome_bin1m.bed
zcat gencode.vM1.annotation.gtf.gz |awk -F "\t" '$3=="gene" && $9~/protein_coding/' |awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$10,$18,$7}' | tr -d '";'  >gencode.vM1.annotation.protein_coding.bed


module load BEDTools
for region in resource/gencode.vM1.annotation.protein_coding_slope2k.bed resource/genome_bin100k.bed 
do
    for meth_file in meth/all_sample.CpG.test.bed
    do
        n1=$(cat ${meth_file}|head -1|awk -F "\t" '{print NF}')
        cat $meth_file|tail -n +2|\
            bedtools intersect  -a <(cut -f1-3 $region|sort -k1,1 -k2,2n ) -b - -wa -wb -sorted |\
            cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c `seq -s ',' 4 1 $n1`|\
            cat <(cat $meth_file|head -1) - |\
            gzip - >${meth_file%.bed}.`basename ${region%.bed}`.mat.gz
    done
done
for i in `ls *CpG.bedGraph`
do
    echo $i `cat $i |grep chr|awk '{mC+=$6;aC+=$5+$6}END{print mC/aC,mC,aC}'` `cat ${i/_CpG.bedGraph/}_CHG.bedGraph |grep chr|awk '{mC+=$6;aC+=$5+$6}END{print mC/aC,mC,aC}'` `cat ${i/_CpG.bedGraph/}_CHH.bedGraph |grep chr|awk '{mC+=$6;aC+=$5+$6}END{print mC/aC,mC,aC}'`
done > meth.sta.txt


for i in `ls *snp.txt.sta`;do echo $i `cat $i |tail -n +2|tr '\n' '\t'`;done >all_sample.CpG.sta
for i in `ls *snp.txt.spikein.sta`;do echo $i `cat $i |tail -n +2|tr '\n' '\t'`;done >all_sample.CpG.spikein.sta
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["SCCAPS_SMPS"]#"C206_N702_i7_2_rev_N504_i5_3_rev" #
REF = config["SCTAPS_SMPS_REF"]
print(SAMPLES)

rule all:
    input:
        # expand("fastq/{sample}_{readDirection}_fastqc.html",sample = SAMPLES, readDirection=['1','2']),
        # expand("align/{sample}.bwa.bam", sample = SAMPLES),
        # expand("align/{sample}.md.bam", sample = SAMPLES),
        # expand("stats/{sample}.mapping.txt", sample = SAMPLES),
        # expand("align/{sample}.md.read.meth.txt.gz", sample=SAMPLES),
        # expand("align/{sample}.md.filter.meth.sta.txt.gz", sample = SAMPLES),
        # expand("align/{sample}.spikeins.bam", sample = SAMPLES),
        # expand("align/{sample}.spikeins.read.meth.txt.gz", sample = SAMPLES),
        # expand("align/{sample}.spikeins.filter.meth.sta.txt.gz", sample = SAMPLES),
        # expand("meth/{sample}.CpG.meth.bed.gz", sample = SAMPLES),
        expand("meth/{sample}.cytosine_report.txt", sample=SAMPLES),
        expand("meth/{sample}.cytosine_report.snp.txt.gz", sample=SAMPLES),
        expand("meth/{sample}.cytosine_report.sta.txt", sample=SAMPLES),
        expand("bam/{sample}.md.clipoverlap.bam", sample=SAMPLES),
        expand("align/{sample}.c.snp.txt.gz", sample=SAMPLES),
        expand("meth/{sample}.cytosine_report.snp.txt.sta", sample=SAMPLES),
        expand("meth/{sample}.cytosine_report.snp.txt.spikein.sta", sample=SAMPLES)
        # expand("meth/{sample}.144hmC_CHH.bedGraph", sample=SAMPLES),
        # expand("meth/{sample}.CpG.gene.bed", sample=SAMPLES),
        # expand("meth/{sample}.CpG.100k.bed", sample = SAMPLES)



################################ PREPROCESS #################################
rule mergefastq: 
    input:
        "code/P230072_rename_fastq.sh"
    output:
        expand("fastq/{{sample}}_{readDirection}.fastq.gz", readDirection=['1','2'])
    params:
        "{sample}"
    shell:
        """
        cmdline=$(grep {params} {input}|tr '\\n' ';')
        echo $cmdline|bash -
        """

rule fastp: 
    input:
        expand("fastq/{{sample}}_{readDirection}.fastq.gz",readDirection=['1','2'])
    output:
        expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    params:
        prefix="fastq/{sample}.fastp",
        fastp=config["fastp"]
    shell:
        """
        {params.fastp} \
            -i {input[0]} \
            -I {input[1]} \
            -o {output[0]} \
            -O {output[1]} \
            -j {params.prefix}.json \
            -h {params.prefix}.html 
        """

rule fastqc: 
    input:
        raw=expand("fastq/{{sample}}_{readDirection}.fastq.gz",readDirection=['1','2']),
        clp=expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    output:
        expand("fastq/{{sample}}_{readDirection}_fastqc.html",readDirection=['1','2'])
    params:
        fastqc=config["fastqc"]
    shell:
        """
        module purge
        module load {params.fastqc}
        fastqc {input.raw[0]}
        fastqc {input.raw[1]}
        fastqc {input.clp[0]}
        fastqc {input.clp[1]}
        """

rule align_fastp: 
    input:
        expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    output:
        "align/{sample}.bwa.bam"
    log:
        "logs/{sample}.bwa.log"
    params:
        ref=REF,
        tmp="{sample}",
        bwa=config["bwa"],
        samtools=config["samtools"]
    threads: 2
    shell:
        """
        module purge
        module load {params.bwa}
        module load {params.samtools}
        (bwa mem -t {threads} {params.ref} {input} |\
        samtools sort -@ 8 -O BAM -T {params.tmp} >{output}) 1>{log} 2>&1
        """


################################ MM9 #################################
rule markdup: 
    input:
        "align/{sample}.bwa.bam"
    output:
        mdbam="align/{sample}.md.bam",
        matrix="align/{sample}.md_report.txt"
    log:
        "logs/{sample}.picard.log"
    params:
        mdtmp=temp("markdup"),
        picard=config["picard"]

    shell:
        """
        module purge
        module load {params.picard}
        (
            java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar MarkDuplicates \
            I={input} \
            O={output.mdbam} \
            M={output.matrix} \
            TMP_DIR={params.mdtmp}) 1>{log} 2>&1
        """


rule mapping: 
    input:
        bam="align/{sample}.md.bam",
        dup="align/{sample}.md_report.txt"
    output:
        "stats/{sample}.mapping.txt"
    params:
        sample="{sample}",
        samtools=config["samtools"]
    shell:
        """
        module purge
        module load {params.samtools}
        mkdir {params.sample}_temp
        q1_nmap=`samtools view -q 1 {input.bam}|cut -f1|sort -u -T {params.sample}_temp |wc -l `
        q10_nmap=`samtools view -q 10 {input.bam}|cut -f1|sort -u -T {params.sample}_temp|wc -l `
        dup_rate=`grep ^LI {input.dup} -A1|tail -n +2|cut -f9`
        proper_nmap=`samtools view -q 10 {input.bam} |awk '$2==99||$2==163'|awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{sum+=$9}}END{{print NR,sum/NR}}'`
        echo -e "sample\\tq1_nmap\\tq10_nmap\\tproper_nmap\\tmean_isize\\tdup_rate" >{output}
        echo -e "{params.sample}\\t$q1_nmap\\t$q10_nmap\\t$proper_nmap\\t$dup_rate" >>{output}
        rm {params.sample}_temp -rf
        """

rule extract_meth: 
    input:
        bam="align/{sample}.md.bam"
    output:
        "align/{sample}.md.read.meth.txt.gz"
    params:
        python=config["python"],
        nd_taps_extract_parallel=config["nd_taps_extract_parallel"]
    shell:
        """
        module purge
        module load {params.python}
        python3 {params.nd_taps_extract_parallel} -b {input.bam} -t 20 -n 100000
        """

rule call_meth: 
    input:
        "bam/{sample}.md.bam"
    output:
        "meth/{sample}.cytosine_report.txt"
    log:
        "logs/{sample}_methcall.log"
    params:
        ref=REF,
        methydackel=config["methydackel"],
        prefix="meth/{sample}",
        samtools=config["samtools"]
    shell:
        """
        (
        {params.methydackel} extract {params.ref} {input} -o {params.prefix} -p 13 \
            --OT 10,118,10,118 --OB 10,118,10,118 --cytosine_report --CHG --CHH 
        module purge
        ) 1>{log} 2>&1
        """

rule call_meth_sta: 
    input:
         "meth/{sample}.cytosine_report.txt"
    output:
        "meth/{sample}.cytosine_report.sta.txt"
    shell:
        """
        grep chr {input} |awk '{{mC[$6]+=$5;uC[$6]+=$4}}END{{for(i in uC)print i,mC[i],uC[i],mC[i]/(mC[i]+uC[i])}}' >{output}
        """

rule bam_clipoverlap: 
    input:
        "bam/{sample}.md.bam"
    output:
        "bam/{sample}.md.clipoverlap.bam"
    shell:
        """
        /users/ludwig/cfo155/cfo155/conda_tools/bin/bam clipOverlap --in {input} --out {output} --poolSize 100000000
        """

rule cg_snp: 
    input:
        "bam/{sample}.md.clipoverlap.bam"
    output:
        "align/{sample}.c.snp.txt.gz"
    params:
        bcftools="BCFtools",
        ref=REF
    shell:
        """
        module purge
        module load {params.bcftools}
        bcftools mpileup -f {params.ref} {input} -a AD|\
            awk '$4=="C"||$4=="G"'|\
            awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{gsub(".*:","",$10);split($10, arr, ",");sum = 0;for (i in arr) sum += arr[i]; print $1,$2,$4,$5,$10,arr[1],sum}}'|gzip - >{output}
        """

rule snp_meth: 
    input:
        snp="align/{sample}.c.snp.txt.gz",
        meth="meth/{sample}.cytosine_report.txt"
    output:
        "meth/{sample}.cytosine_report.snp.txt.gz"
    shell:
        """
        join -1 1 <(awk 'BEGIN{{OFS="\\t"}}{{print $1"_"$2,$0}}' {input.meth}|sort -k1,1) \
             -2 1 <(zcat {input.snp}|awk 'BEGIN{{OFS="\\t"}}{{print $1"_"$2,$0}}' |sort -k1,1) -a1 -t$'\\t'|gzip - >{output}
        """

rule snp_meth_sta:
    input:
        "meth/{sample}.cytosine_report.snp.txt.gz"
    output:
        "meth/{sample}.cytosine_report.snp.txt.sta"
    log:
        "{sample}.log"
    shell:
        """
        (
        echo -e "context\\tmpileup_mC\\tmpileup_aC\\tmpileup_rC\\tmethydackel_mC\\tmethydackel_aC\\tmethydackel_rC" >{output}
        zcat {input}| awk '$9~/chr/'| awk '$12=="<*>"||($11=="C"&&$12=="T,<*>")||($11=="G"&&$12=="A,<*>")'|\
            awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{uC1[$7]+=$14;aC1[$7]+=$15;uC2[$7]+=$5;aC2[$7]+=$6+$5}}END{{for(i in uC1)print i,uC1[i],aC1[i],1-uC1[i]/aC1[i],uC2[i],aC2[i],1-uC2[i]/aC2[i]}}' >> {output} 
        ) >{log} 2>&1
        """

rule snp_meth_sta_spikein:
    input:
        "meth/{sample}.cytosine_report.snp.txt.gz"
    output:
        "meth/{sample}.cytosine_report.snp.txt.spikein.sta"
    log:
        "{sample}.log"
    shell:
        """
        (
        echo -e "context\\tmpileup_mC\\tmpileup_aC\\tmpileup_rC\\tmethydackel_mC\\tmethydackel_aC\\tmethydackel_rC" >{output}
        zcat {input}| awk '$9~/144hmC/&&($10=="87"||$10=="92"||$10=="100"||$10=="110"||$10=="117")'|\
            awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{uC1[$7]+=$14;aC1[$7]+=$15;uC2[$7]+=$5;aC2[$7]+=$6+$5}}END{{for(i in uC1)print "144hmC:"i,uC1[i],aC1[i],1-uC1[i]/aC1[i],uC2[i],aC2[i],1-uC2[i]/aC2[i]}}' >> {output} 
        zcat {input}| awk '$9~/J02459.1/'| awk '$12=="<*>"||($11=="C"&&$12=="T,<*>")||($11=="G"&&$12=="A,<*>")'|\
            awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{uC1[$7]+=$14;aC1[$7]+=$15;uC2[$7]+=$5;aC2[$7]+=$6+$5}}END{{for(i in uC1)print "J02459.1:"i,uC1[i],aC1[i],1-uC1[i]/aC1[i],uC2[i],aC2[i],1-uC2[i]/aC2[i]}}' >> {output} 
        zcat {input}| awk '$9~/unmodified_2kb/'| awk '$12=="<*>"||($11=="C"&&$12=="T,<*>")||($11=="G"&&$12=="A,<*>")'|\
            awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{uC1[$7]+=$14;aC1[$7]+=$15;uC2[$7]+=$5;aC2[$7]+=$6+$5}}END{{for(i in uC1)print "unmodified_2kb:"i,uC1[i],aC1[i],1-uC1[i]/aC1[i],uC2[i],aC2[i],1-uC2[i]/aC2[i]}}' >> {output} 
        ) >{log} 2>&1
        """

rule cpg_filter:
    input:
        methcall="align/{sample}.md.filter.meth.sta.txt.gz",
    output:
        cpg="meth/{sample}.CpG.meth.bed"
    params:
        sample="{sample}",
        cpg=config["cgpos"],
        bedtools="BEDTools"
    shell:
        """
        module purge
        module load {params.bedtools}
        zcat {input.methcall}|tail -n +2|sort -k1,1 -k2,2n|awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$2+1,$3,$4,$5}}'|\
            bedtools intersect -a {params.cpg} -b - -sorted -wao|\
            awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{if($4!=".")printf"%s\\t%d\\t%d\\t%.1f\\t%d\\t%d\\n",$1,$2,$3,$9,$7,$8;else print $1,$2,$3,"*","0","0"}}'|\
            cat <(echo -e "CHROM\\tSTART\\tEND\\t{params.sample}_RATE\\t{params.sample}_MOD\\t{params.sample}_TOTAL") - > {output.cpg}
        """

rule region_meth_1:
    input:
        "meth/{sample}.CpG.meth.bed.gz"
    output:
        bins="meth/{sample}.CpG.100k.bed",
        gene_2k="meth/{sample}.CpG.gene_2k.bed"
    params:
        bins=config["bins"],
        gene_2k=config["gene_2k"],
        bedtools="BEDTools"
    shell:
        """
        module purge
        module load {params.bedtools}
        bedtools intersect  -a <(cut -f1-3 {params.bins} |sort -k1,1 -k2,2n ) -b <(zcat {input}) -wa -wb -sorted |\
            cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c 5,6|\
            cat <(cat {input}|head -1|cut -f1-3,5,6) - >{output.bins}
        bedtools intersect  -a <(cut -f1-3 {params.gene_2k} |sort -k1,1 -k2,2n ) -b <(zcat {input}) -wa -wb -sorted |\
            cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c 5,6|\
            cat <(cat {input}|head -1|cut -f1-3,5,6) - >{output.gene_2k}
        """

rule region_meth_2:
    input:
        "meth/{sample}.CpG.meth.bed.gz"
    output:
        "meth/{sample}.CpG.gene.bed"
    params:
        region=config["gene"],
        bedtools="BEDTools"
    shell:
        """
        module purge
        module load {params.bedtools}
        bedtools intersect  -a <(cut -f1-3 {params.region} |sort -k1,1 -k2,2n ) -b {input} -wa -wb -sorted |\
            cut -f1-3,7-|bedtools groupby -i - -g 1,2,3 -o sum -c 5,6|\
            cat <(cat {input}|head -1|cut -f1-3,5,6) - >{output}
        """

rule gzip:
    input:
        "meth/{sample}.CpG.meth.bed"
    output:
        "meth/{sample}.CpG.meth.bed.gz"
    shell:
        """
        gzip {input}
        """


######################## spike-ins ########################
rule extract_spikein: 
    input:
        "align/{sample}.bwa.bam"
    output:
        "align/{sample}.spikeins.bam",
    params:
        ref=REF,
        samtools=config["samtools"]
    shell:
        """
        module purge
        module load {params.samtools}
        samtools view -bS {input} -L <(grep -v "chr" {params.ref}.fai|awk 'BEGIN{{OFS="\\t"}}{{print $1,"0",$2}}') >{output}
        """

rule spikein_calls: 
    input:
        "align/{sample}.spikeins.bam"
    output:
        "align/{sample}.spikeins.read.meth.txt.gz"
    params:
        ref=REF,
        python=config["python"]
    shell:
        """
        module purge
        module load {params.python}
        python3 code/nd_taps_extract_parallel.py -b {input} -t 10 -n 100000
        """

rule spikein_trim: 
    input:
        "align/{sample}.spikeins.read.meth.txt.gz"
    output:
        "align/{sample}.spikeins.filter.meth.sta.txt.gz"
    params:
        python=config["python"]
    shell:
        """
        module purge
        module load {params.python}
        python3 code/nd_taps_convert_call.py -r {input} -s 10 -e 118 -c True
        """


rule call_meth_spikein: 
    input:
        bam="align/{sample}.spikeins.bam"
    output:
        "meth/{sample}.144hmC_CHH.bedGraph"
    log:
        "logs/{sample}_methcall.log"
    params:
        ref=REF,
        methydackel=config["methydackel"],
        prefix="meth/{sample}.144hmC",
        samtools=config["samtools"]
    shell:
        """
        (
        module load {params.samtools} 
        samtools index {input.bam}
        {params.methydackel} extract {params.ref} {input.bam} -o {params.prefix} -p 13 -r 144hmC \
            --OT 10,118,10,118 --OB 10,118,10,118 --CHG --CHH 
        ) 1>{log} 2>&1
        module purge
        """