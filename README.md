# Targeted-sequencing Pipeline 
## This is a pipeline for the targetsequencing of 48 cancer-ralated gene
### Time:2018.11.07
##
## What u need ?
#### Softwore: FastQC,P-trimmer[coding by xlzh],Samtools,BWA,Varscan...
#### Script: Dep_AF.py, Integrat.r
#### Files : HG38 human reference genome, hg38amplicon.txt, hg38ampinfor.txt
##
## Workflow
#### Quality Control at first by FastQC software whatever data you have, u must to know the quality of data and it is better for follow-up work
### Steps:
#### 1. Remove primer
#### 2. Mapping[BWA, Ref hg38]
#### 3. Remove multiple-mapping reads 
#### 4. Calculate mutation frequency of each site[VasScan, Q30]
#### 5. Make AFdata
### Result: 
#### 1. Basic Pictures about AlleleFrequency and AmpliconDepth, 
#### 2. Basic Rdata about AF and Amplicon.
#### 3. If you want to call SNV, you can input the outputing final BAM files as an input for PlasamMutationDetector.
##
## Usage
## bash TargetPanel.sh In Index Ref Out
## 
#### In: The input path which including all the rowdata of fastq.Formation: XXX_cfDNA_R1.fq.gz, XXX_cfDNA_R2.fq.gz
#### Index: The mapping index of reference genome
#### Ref: The reference genome
#### Out: The output path which including Trimfq, bampath, filterbampath, snvfile and so on.
## 
## Example of Usage
bash TargetPanel.sh /home/rowdata /home/Documents/hg38/hg38.fastq /home/Documents/hg38/hg38.fastq /home/Result


















