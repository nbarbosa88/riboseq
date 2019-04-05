#!/bin/bash

#author UrMi
#pipeline to run ribocode on .SRR files
#input directory containing SRR files

#Step 0: Load all required modules, files, prereqs. These are
# Modules: sra-toolkit, bowtie2, STAR, RiboCode
# Files: Humangenome.fa corresponding HumanAnnotation.gtf SRR files, rRNA seq database
# Indexes: STAR index, bowtie2 rRNA index

module load bowtie2/2.3.1-py3-ge4lv4s
module load star/2.5.3a-rdzqclb
module load python/3.6.3-u4oaxsb
module load py-pip/9.0.1-py3-dpds55c


#copy needed files to scratch for faster access
rRNAIndex=/work/LAS/mash-lab/usingh/riboSeq/sampleData/humanData/bowtieIndex/rRNAindex
STARindex=/work/LAS/mash-lab/usingh/riboSeq/sampleData/humanData/human_STARindex
RiboCodeAnnot=/work/LAS/mash-lab/usingh/riboSeq/sampleData/humanData/RiboCode_annot
proc=16


#Step 1: Convert .SRA files to FASTQ. Use fastq-dump to control quality
file_dir=$1
#make list of all SRA files in input directory
file_list=($file_dir/*.sra)

#run this many fastq dump in parallel
size=10
i=0
for f in "${file_list[@]}"; do
	echo "$f"
	#run fastqdump
	this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
	echo $this_fname
	echo "fastq-dump --readids --split-files --dumpbase --skip-technical --clip --read-filter pass --outdir $file_dir $f && rm -f "$f" &"
	fastq-dump --readids --split-files --dumpbase --skip-technical --clip --read-filter pass --outdir $file_dir $f && rm -f "$f" & 
	v=$(( $(($i+1)) % $size)) 
	if [ "$v" -eq "0" ]; then
  		echo $i
		echo $v
		echo "waiting..."
		wait
	fi
	i=$(($i+1))
done
echo "finally waiting for fastq-dump to finish..."
wait

########Check if files are RPF not RNA seq############

#Step 2: Remove rRNA reads from FASTQ using bowtie or Sortmerna
echo "Filtering rRNA using bowtie2"

file_list=($file_dir/*pass_1.fastq)

#run for all fastq files
for f in "${file_list[@]}"; do
	this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
	echo $this_fname
	
	echo "bowtie2 -p $proc --norc --un "$file_dir/$this_fname"_norRNA.fastq -q $f -x $rRNAIndex -S "$file_dir/$this_fname"_bt2Out.SAM"
	bowtie2 -p $proc --norc --un "$file_dir/$this_fname"_norRNA.fastq -q $f -x $rRNAIndex -S "$file_dir/$this_fname"_bt2Out.SAM


	if [ $? -ne 0 ]; then
                fail_flag=true
                echo "FAILED BOWTIE2 FOR " "$this_fname"
                echo "$this_fname" >> "$file_dir"/failed_bowtie.log
                failed_salmon+=("$this_fname")
                continue
        fi
	#remove unwanted files
	rm -f "$f"
	rm -f "$file_dir/$this_fname"_bt2Out.SAM									
done



#output of above step is files with no alignment to rRNA ending with name _norRNA.fastq


#Step 3: Align clean reads to human genome using STAR

echo "Running Star..."
file_list=($file_dir/*_norRNA.fastq)

for f in "${file_list[@]}"; do
        this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
        echo $this_fname
	
	echo "STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir $STARindex --readFilesIn $f --outFileNamePrefix "$file_dir/$this_fname" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd"
	STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir $STARindex --readFilesIn $f --outFileNamePrefix "$file_dir/$this_fname" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd

        if [ $? -ne 0 ]; then
                fail_flag=true
                echo "FAILED STAR FOR " "$this_fname"
                echo "$this_fname" >> "$file_dir"/failed_star.log
                failed_salmon+=("$this_fname")
                continue
        fi
        #remove unwanted files
        rm -f "$f"
	#rm -f 
done

#output of above step is BAM files with name "$this_fnamei"  SRR1234_pass_1_norRNAAligned.toTranscriptome.out.bam

#exit

#Step 4: Run RiboCode to identify translated ORFs and save results for each SRR file

echo "Running RiboCode"
file_list=($file_dir/*toTranscriptome.out.bam)

for f in "${file_list[@]}"; do
	echo $f
        #this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
	this_fname=$(echo "$f" | tr "_" "\t" | awk '{print $1}')
        echo $this_fname
	
	#step 4a run metaplots
	echo "metaplots -a $RiboCodeAnnot -r "$f" -o "$this_fname""
	metaplots -a $RiboCodeAnnot -r "$f" -o "$this_fname"

	#output of metaplots is config file used in next step. Output name is "$this_fname"_pre_config.txt

	#step 4b run RiboCode
	echo "RiboCode -a $RiboCodeAnnot -c "$this_fname"_pre_config.txt -l no -g -o "$this_fname"_RiboCoderesult"
	RiboCode -a $RiboCodeAnnot -c "$this_fname"_pre_config.txt -l no -g -o "$this_fname"_RiboCoderesult
	

        if [ $? -ne 0 ]; then
                fail_flag=true
                echo "FAILED RC FOR " "$this_fname"
                echo "$this_fname" >> "$file_dir"/failed_ribocode.log
                failed_salmon+=("$this_fname")
                continue
        fi

	#cleanup
	rm -f "$this_fname"_pass_1_norRNAAligned.toTranscriptome.out.bam
	rm -f "$this_fname"_pass_1_norRNAAligned.sortedByCoord.out.bam
done





