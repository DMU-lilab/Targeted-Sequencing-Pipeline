# Targeted-Sequencing-Pipeline
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

###Usage: bash TargetPanel.sh In Index Ref Out
###For example: bash TargetPanel.sh /home/rawfastq /home/Document/hg38/hg38.fastq /home/Document/hg38/hg38.fastq  /home/Result 

