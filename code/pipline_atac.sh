source atac.sh

#该流程可用于小鼠和人类chip-seq和atac-seq数据的分析，且默认使用的数据为双端测序
work_dir=$(pwd)
cd ..
all_dir=$(pwd)
cd ${work_dir}

group="Gli2CT1 Gli3CT1 K27m3CT1"
species="mouse" #human or mouse OR melanogaster
reference_dir="/data1/shenluemou/reference/mapping/bowtie2/mouse/grcm38/mouse_grcm38"
gtf_dir="/data1/shenluemou/reference/gtf/mouse/grcm38/mouse_grcm38.gtf"
bed_dir="/data1/shenluemou/reference/bed/mouse/grcm38/mouse_grcm38.bed"
homer_mode="factor"


mkdir $all_dir/result
#--------quality_control---------------------------------
mkdir $all_dir/result/1qc
echo "-------------------begin 1qc----------------------"
for name in ${group};do
	quality_control ${all_dir}/data ${name} ${all_dir}/result/1qc
	echo "finish the 1qc of "${name}
done
echo "-------------------finish all 1qc---------------------"

#-------------map---------------------------------------
mkdir $all_dir/result/2map
echo "-------------------begin map-----------------------"
for name in ${group};do
	fqtobw ${all_dir}/result/1qc ${name} $reference_dir $all_dir/result/2map
	echo "finish the 2map of "${name}
done
echo "-------------------finish all 2map------------------"


#----------------begin 3bamqc----------------------------
mkdir $all_dir/result/3bamqc
echo "--------------------begin bamqc---------------------"
for name in ${group};do
	bamqc $all_dir/result/2map ${name} $all_dir/result/3bamqc
	echo "finish the 3bamqc of "${name}
done
multiBigwigSummary  bins -b $all_dir/result/3bamqc/*bw -o $all_dir/result/3bamqc/matrix.npz -p 8
plotCorrelation -in $all_dir/result/3bamqc/matrix.npz  \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
--plotFileFormat pdf \
-o $all_dir/result/3bamqc/matrix.pdf   \
--outFileCorMatrix $all_dir/result/3bamqc/matrix.tab
echo "--------------------finish all 3bamqc---------------"


#-------------------4callpeak--------------------------------
mkdir $all_dir/result/4callpeak
echo "----------------------begin callpeak-----------------"
for name in ${group};do
	mkdir $all_dir/result/4callpeak/$name
	callpeak $all_dir/result/3bamqc ${name} $all_dir/result/4callpeak/$name $homer_mode
	echo "finish the 4callpeak of "${name}
done
echo "-----------------------finish all 4callpeak-----------"


#-------------------5motif-------------------------------------
mkdir $all_dir/result/5motif
echo "---------------------begin motif----------------------"
for name in ${group};do
	mkdir $all_dir/result/5motif/$name
	motif $all_dir/result/4callpeak/$name ${name} $species $all_dir/result/5motif/$name
	echo "finish the 5motif analysis of "${name}
done
echo "--------------------finish all motif------------------"


#-------------------6chipseeker-------------------------------
mkdir $all_dir/result/5chipseeker
mkdir $all_dir/result/5chipseeker/peak
echo "-------------------begin chipseeker-------------------"
for name in ${group};do
	chipseeker $all_dir/result/4callpeak/$name $name $all_dir/result/5chipseeker/peak $all_dir/result/5chipseeker $species
done
Rscript chipseeker.R -w $all_dir/result/5chipseeker -p $all_dir/result/5chipseeker/peak -s $species > $all_dir/result/5chipseeker/chip.log


#------------------

































