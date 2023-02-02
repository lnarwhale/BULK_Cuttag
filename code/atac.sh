#!/bin/bash


function quality_control(){
	#--function:quality control for chip-seq or atac-seq
	#--input: 1.fastq_dir; 2.name ;3.out_dir
	#--output: the fastq after quality_control
	#--need: 
	echo "the fastq is "$1/$2
	echo "the out_dir is "$3
	/data1/shenluemou/biosoft/fastp \
	--in1 $1/$2"_R1.fq" \
	--in2 $1/$2"_R2.fq" \
	--out1 $3/$2"_R1.fq" \
	--out2 $3/$2"_R2.fq" \
	--html $3/$2".html" \
	--trim_poly_g --poly_g_min_len 5 \
	--trim_poly_x --poly_x_min_len 5 \
	--cut_front --cut_tail --cut_window_size 4 \
	--qualified_quality_phred 15 \
	--low_complexity_filter \
	--complexity_threshold 30 \
	--length_required 30 \
	--thread 4
}

function fqtobw(){
	#--function:fastq to bw
	#--input: 1.fastq_dir; 2.name; 3.reference_dir&name(index); 4.out_dir
	#--output: sam bam bw mapping_picture
	#--need: (eval work_dir(include plot_function.R); map_stat.R)
	echo "the fastq is "$1/$2
	echo "the reference is "$3
	echo "the out_dir is "$4
	bowtie2 -p 10 -x $3 -1 $1/$2"_R1.fq" -2 $1/$2"_R2.fq" -S $4/$2".sam" > $4/$2".log" 2>&1
	samtools view -bS $4/$2".sam" > $4/$2".bam"
	samtools sort -@ 8 $4/$2".bam" -o $4/$2"_sorted.bam"
	samtools index $4/$2"_sorted.bam"
	bamCoverage -b $4/$2"_sorted.bam" -o $4/$2".bw"
	mkdir $4/tmp
	cat $4/$2".log" | grep -E ') aligned concordantly' > $4/tmp/tmp_mapstat.txt
	sed -i 's/[a-zA-Z]//g' $4/tmp/tmp_mapstat.txt
	sed -i 's/[\(\)\%]//g' $4/tmp/tmp_mapstat.txt
	sed -i 's/[ ][ ]*/ /g' $4/tmp/tmp_mapstat.txt
	Rscript map_stat.R -i $4/tmp/tmp_mapstat.txt -w $4/tmp -f $work_dir/"plot_function.R" -s $2 > $4/tmp/$2"_forstat.log" 2>&1
	mv $4/tmp/tmp_map.txt $4/$2"_map.txt"
	mv $4/tmp/tmp_mappie.pdf $4/$2"_mappie.pdf"
	rm -rf $4/tmp
}

function bamqc(){
	#--function:bam quality control
	#--input: 1.bam_dir; 2.name; 3.out_dir; 
	#--output: bam after qc
	#--need:
	touch $1/$2"_bamqc.txt"
	/data1/shenluemou/biosoft/sambamba/bin/sambamba markdup -r $1/$2"_sorted.bam" $1/$2"_markdup.bam"
	yuan=$(samtools view -h $1/$2"_sorted.bam" | grep -v '@' | awk '{print $1}' | sort | uniq | wc -l)
	markdup=$(samtools view -h $1/$2"_markdup.bam" | grep -v '@' | awk '{print $1}' | sort | uniq | wc -l)
	samtools view -h $1/$2"_markdup.bam" | grep -v 'MT' | samtools view -bS -o $3/$2".bam"
	de_mt=$(samtools view -h $3/$2".bam" | grep -v '@' | awk '{print $1}' | sort | uniq | wc -l)
	samtools index $3/$2".bam"
	bamCoverage -b $3/$2".bam" -o $3/$2".bw"
	echo "bam_yuan,"$yuan >> $1/$2"_bamqc.txt"
	echo "bam_rmdup,"$markdup >> $1/$2"_bamqc.txt"
	echo "bam_demt,"$de_mt >> $1/$2"_bamqc.txt"
	Rscript plot_log.R -l $1 -n $2
}

function callpeak(){
	#--function:callpeak
	#--input: 1.bam_dir; 2.name; 3.out_dir; 4.type of callpeaking 	
	#--output: macs2 callpeak result
	#--need:
	makeTagDirectory $3 $1/$2".bam" 
	findPeaks $3 -style $4 -o $3/$2"_callpeak_tmp.txt"
	cat $3/$2"_callpeak_tmp.txt" | grep -v "#" | awk '{print $1"\t""chr"$2"\t"$3"\t"$4"\t"$5}' > $3/$2"_callpeak.txt"
}


function motif(){
	#--function:motif find
	#--input:1.callpeak_dir; 2.name; 3.species; 4.out_dir
	#--output:
	#--need:
	if [ $3 = "human" ];then
		ref=hg19
	fi
	if [ $3 = "mouse" ];then
		ref=mm10
	fi
	if [ $3 = "melanogaster" ];then
		ref=dm6
	fi
	findMotifsGenome.pl $1/$2"_callpeak.txt" $ref $4 -size 200 -mask
}



function chipseeker(){
	#--function:chipseeker
	#--input: 1.callpeak_dir; 2.name; 3.out_peak_call; 4.work_dir; 5.species
	#--output:
	#--need:
	cat $1/$2"_callpeak_tmp.txt" | grep -v "#" | awk '{print "chr"$2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$8}' > $3/$2"_tochipseeker.txt"
	#Rscript chipseeker.R -w $4 -p $3 -s $5 > $4/chip.log
}






