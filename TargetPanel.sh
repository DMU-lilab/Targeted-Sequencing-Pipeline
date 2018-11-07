#!/bin/bash
###### Time:2018.06.01 
###### This is a pipeline for the targetsequencing of 48 cancer-ralated gene

###### What u need
#       ||||||||
###### Softwore: FastQC,PrimerTrim[coding by xlzh],Samtools,BWA,Varscan,Annovar
###### file:HG38 human reference genome, 2158 pairs of primer information

###### Workflow
#      ||||||||
###### Quality Control at first by FastQC software whatever data you have, u must to know the quality of data and it is better for follow-up work
###### Remove primer ---> Mapping[BWA mem, Ref hg38] ---> Remove multiple-mapping reads ---> Calculate mutation frequency of each site[VasScan, Q30] ---> Annotation[Annovar]

############### In: The input path which including all the rowdata of fastq.Formation: XXX_cfDNA_R1.fq.gz, XXX_cfDNA_R2.fq.gz
############### Index: The mapping index of reference genome
############### Ref: The reference genome
############### Out: The output path which including Trimfq, bampath, filterbampath, snvfile and so on.

In=$1 
Index=$2
Ref=$3
Out=$4

###############Step1:Remove primer
###############Step2:Mapping
###############Step3:Remove mutliple-reads
###############Step4:Caculate AF

echo "Step1:Remove primer==========================================================================================================="
Trimfq=$Out'/Trimfq/'
if [ ! -d $Trimfq ];then mkdir -p $Trimfq;fi   ###Create a new folder
Newfq=$Trimfq"Newfq"
Ampliconf=$Trimfq"Ampliconf"
if [ ! -d $Newfq ];then mkdir -p $Newfq;fi
if [ ! -d $Ampliconf ];then mkdir -p $Ampliconf;fi

cd $In;fa=`ls *.gz|cut -d "_" -f 1,2|uniq`
for i in $fa;do
        echo "-----remove primer of $i--------"
	R1="_R1.fq.gz"
        R2="_R2.fq.gz"
        /home/jctian/software/Test/pTrimmer-master/pTrimmer-1.3.1 \
		-s pair \
                -a /home/jctian/Documents/Panelprimer/hg38amplicon.txt \
                -f $i$R1 \
                -r $i$R2 \
                -o $Trimfq
	mv $Trimfq"Summary.ampcount" $Trimfq$i".count"
        echo "-----------$i is finished---------------"
done
##
mv $Trimfq*.count $Ampliconf
mv $Trimfq*.fq $Newfq


echo "Step2:Mapping================================================================================================================="
bampath=$Out'/bampath/'
if [ ! -d $bampath ];then mkdir -p $bampath;fi

fr1="_trim_R1.fq"
fr2="_trim_R2.fq"
cd $Newfq;fa=`ls|cut -d '_' -f 1,2|uniq`
for i in $fa;do
	echo "--------map file $i---------"
	bwa mem -t 10 -M $Index \
		$i$fr1 $i$fr2 \
		|samtools view -bS - >  ${bampath}"tmp.bam"
		samtools sort ${bampath}"tmp.bam" -o ${bampath}$i".sort.bam"
		samtools index ${bampath}$i".sort.bam"
		rm ${bampath}"tmp.bam"
done

echo "Step3:Filter XA tag,remove multiple hit======================================================================================="
filterbampath=$Out'/filterbampath'
if [ ! -d $filterbampath ];then mkdir -p $filterbampath;fi

for f in $(ls $bampath);do
	if [ ${f##*.} == "bam" ];then
		outbname=${f%%.*}".bam"
		echo "filtering" $f
		samtools view -b -x XA -F 0x100 $bampath"/"$f > $filterbampath"/"$outbname
		samtools sort -@ 20 $filterbampath"/"$outbname -o $filterbampath"/"${f%%.*}".sort.bam"
		samtools index $filterbampath"/"${f%%.*}".sort.bam"
		rm $filterbampath"/"$outbname
	fi
done


echo "Step4:Caculate the AF of every site============================================================================================"
AFfile=$Out'/AFfile/'
if [ ! -d $AFfile ];then mkdir -p $AFfile;fi

cd $filterbampath;fa=`ls *.bam`
for i in $fa;do
        echo "============$i=============="
	m=`echo $i| cut -d '.' -f 1`
        java -Xmx40G -jar /home/jctian/software/Varscan/VarScan.v2.3.9.jar mpileup2cns <(samtools mpileup -d 500000 -f $Ref -Q 0 $i) \
                --min-coverage 1 \
                --min-reads2 0 \
                --p-value 1 \
                --min-avg-qual 30 \
                --min-var-freq 0 \
                --output > $AFfile$m
	echo "==========$i finished========="
done

echo "Step5:Integrate all data========================================================================================================"
table=$Out"/table"
if [ ! -d $table ];then mkdir -p $table;fi

cd $AFfile;fa=`ls`
        m=`echo $fa`
	##### use a python script to combine all files in one
        python /home/jctian/20180819test_for_newpipeline/Dep_AF.py $m > $table"/Allsites.txt"
	##### replace NM(NM means no data)
	sed -i 's/NM/0/g' $table"/Allsites.txt"
	##### remove "%" for better calculation in R
	sed -i 's/%//g' $table"/Allsites.txt" 
echo "----------finished------------"

#######this is for integrating all sites in R and drawing any pictures
TXT=$table"/TXT"
Routcome=$Out"/Rdata"
if [ ! -d $Routcome ];then mkdir -p $Routcome;fi
#Picture=$table"/Picture"

echo "Step6:Creat a Rdata of sites and amplicon, and draw basic pictures==============================================================="
	Rscript Integrat.r -i $table"/Allsite.txt" -p /home/jctian/Documents/Panelprimer/hg38ampinfor.txt -a $Ampliconf -o $Routcome 











###-------------------------------------------------------------- CUTTING LINE-----------------------------------------------------------------------------
##########################################################################################################################################################
:<<!                                                                                                                                                     
notExce(){
########################################################################################################
#######This steps for the SNV calling
echo "Step4.2:Select the SNV========="
snvfile=$Out'/snvfile/'
if [ ! -d $snvfile ];then mkdir -p $snvfile;fi

cd $filterbampath;fa=`ls *.bam`
for i in $fa;do
	echo "============$i=============="
	####mutiple bam running
	m=`ls *.bam|grep $i"_"`
	java -Xmx40G -jar /home/jctian/software/Varscan/VarScan.v2.3.9.jar mpileup2cns <(samtools mpileup -d 500000 -f $Ref -Q 0 $m)
		--min-coverage 200 \
		--min-reads2 10 \
		--p-value 1 \
		--min-avg-qual 30 \
		--min-var-freq 0 \
		--output-vcf 1 > $snvfile$i".vcf"
done


echo "Step5.2:Annotation by Annovar======"
annout=$Out'/annout'
if [ ! -d $annout ];then mkdir -p $annout;fi

cd ${snvfile};fa=`ls`
for i in $fa;do
        m=`echo $i|cut -d '.' -f 1`
        perl ~/software/annovar/table_annovar.pl -vcfinput ${i} ~/software/annovar/hg38 \
                -buildver hg38 -out ${annout}"/"${m} \
                -remove -protocol refGene,avsnp147,exac03,1000g2015aug_eas \
                -operation g,f,f,f -nastring .
done

## move annovar files
if [ ! -d $annout"/txt" ];then mkdir $annout"/txt";fi
if [ ! -d $annout"/vcf" ];then mkdir $annout"/vcf";fi
if [ ! -d $annout"/temp" ];then mkdir $annout"/temp";fi

mv $annout"/"*'.hg38_multianno.txt' $annout"/txt"
mv $annout"/"*'.hg38_multianno.vcf' $annout"/vcf"
mv $annout"/"*'.'* $annout"/temp"
}
!
############################################################################################################################################################
