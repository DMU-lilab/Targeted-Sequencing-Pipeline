#!/usr/bin/env Rscript
########################################################################################
## LeeLab                                                                             ##
## 2018-08-18                                                                         ##
##                                                                                    ##
## usage: Rscript Integrat.r -i input.file /                                          ##
##                           -p primer.file /                                         ##
##                           -a amp_path /                                            ##
##                           -o outpath                                               ##
########################################################################################


library(data.table)
library(reshape2)
library(ggplot2)
library(optparse) 
library(gridExtra)


Transtable <- function(pathfile, pri_path){
 ### 读入文件
  table_all <- fread(pathfile, sep="\t", header=T)
 ### 转变格式
  table_m <- melt(table_all,id=1:2) 
 ### 分割列
  list_v <- strsplit(as.character(table_m$variable), "_")
  list_v <- data.frame(do.call(rbind, list_v))
  table_m[,':='(Patients=list_v$X1, Type=list_v$X2, df=list_v$X3)]
 ### 转变格式
  table_f <- dcast(table_m, Chrom+Pos+Patients+Type~df, value.var = "value")
 ### 重新排序
  table_f <- data.table(table_f)[,c(1,2,3,5,6,4)]
  colnames(table_f)[1] <- "Chr"
### 输入primer信息
  primer <- fread(pri_path,sep="\t",header=T)
  pn <- primer[,.(Chr=Chr, Pos=seq(Insert_Start, Insert_Stop)),by=c("Gene","AMPN")]
### 合并
  newt <- merge(table_f, pn, by=c("Chr", "Pos"))
  return(newt)
}

#################################
## Amplicon Information  ##
## Creat a summary table ##
###########################	
Amp_Table <- function(amp_path, Primer){
    ## 读入primer文件，并转换格式
	primer <- fread(Primer, sep="\t", header=T)
	pm <- primer[,.(AuxInfo=paste(Gene, Chr, Amplicon_Start, 
	                              Insert_Start, Insert_Stop, 
	                              Amplicon_Stop,sep=",")),by=c("AMPN","Gene")]
######> pm
##	AMPN    Gene                                           AuxInfo
##	1:     MPL-A01     MPL      MPL,chr1,43337816,43337835,43337933,43337953
##	2:     MPL-A02     MPL      MPL,chr1,43338041,43338062,43338145,43338168
##	3:     MPL-A03     MPL      MPL,chr1,43338088,43338112,43338201,43338223
##	4:     MPL-A04     MPL      MPL,chr1,43338510,43338532,43338623,43338646
##	5:     MPL-A05     MPL      MPL,chr1,43338595,43338615,43338711,43338734
##	---                                                                      
##	2154: SMARCB1-A19 SMARCB1 SMARCB1,chr22,23825301,23825322,23825417,23825439
##	2155: SMARCB1-A20 SMARCB1 SMARCB1,chr22,23825385,23825407,23825491,23825516
##	2156: SMARCB1-A21 SMARCB1 SMARCB1,chr22,23833531,23833555,23833647,23833668
##  2157: SMARCB1-A22 SMARCB1 SMARCB1,chr22,23833560,23833578,23833674,23833692
##  2158: SMARCB1-A23 SMARCB1 SMARCB1,chr22,23834064,23834084,23834176,23834196
	####
  Amp_count <- list.files(amp_path)
  Amp <- data.frame()
     ## 整合所有Amplicon文件
        for ( i in Amp_count){
            amp_e <- fread(paste(amp_path,i,sep="/"), sep="\t", header=T)
            amp_e$Sample <- sub(".count", "", i) 
            Amp <- rbind(Amp, amp_e)
        }
     ##
  Amp[,':='(Log2=round(log(AmpCount+1,2),2), Log10=round(log(AmpCount+1,10),2))]
     ##
	Amall2 <- data.frame(do.call(rbind,strsplit(as.character(Amp$Sample), "_")))
	colnames(Amall2) <- c("Patients", "Type")
	Amp <- cbind(Amp, Amall2)
	Amp[, Sample := NULL]
     ##
	m <- merge(Amp,pm,by=c("AuxInfo"))
  return(m)
} 
 
############################################################
## For drawing some pictures                     ##
## Amplicon pictures :Boxplot, Violin, Histogram ##
###################################################
Amp_Draw <- function(AMP){
	## Boxplot
	Amp_boxplot <- ggplot(AMP,aes(Gene,Log2,fill=Gene)) + geom_boxplot() +
		#scale_x_discrete(limits=c("WBC","cfDNA1","cfDNA2","cfDNA3","tissue")) + 
		ggtitle("Amplicon Depth") + xlab("Gene") + ylab("Log2 value") +
	  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1)) +
	  guides(fill = FALSE) + ylim(5,15)
	## Violoin
	Amp_violin <- ggplot(AMP,aes(Type,Log2,fill=Type)) + geom_violin() +
		#scale_x_discrete(limits=c("WBC","cfDNA1","cfDNA2","cfDNA3","tissue")) +
		ggtitle("Amplicon Depth") + xlab("Type") + ylab("Amplicon Depth") +
		guides(fill=guide_legend(title="Legend_Title")) +
		theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
	## Histogram
	Amp_histogram <- ggplot(AMP,aes(Log2,fill=Type)) + geom_histogram(bins=90) +
		ggtitle("Amplicon Depth") + xlab("AmpliconDepth") + ylab("Count") + 
		guides(fill=guide_legend(title="Legend_Title")) + 
		theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) + 
		facet_grid(~Patients)
	return(grid.arrange(Amp_boxplot,Amp_violin,Amp_histogram,nrow=2,ncol=2))
}

####################################################
### Allele Fraction : Scatter plot; ECDF ##
###########################################
AF_Draw <- function(AF){
	  ## Scatterplot
	AF_scatter <- ggplot(unique(AF[,1:7]), aes(Type, Freq, colour=Type)) + geom_jitter() +
		ggtitle("AllSample Allele Frequency") + xlab("Type") + ylab("Allele Frequency") +
		#scale_color_manual(values = c("red","blue","orange","green","black"))+
	  theme_bw() +
		theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
	
    ## ECDF plot
	AF_ecdf <- ggplot(AF,aes(Freq,colour=Type)) + stat_ecdf() + ylim(0.999,1) +
		ggtitle("AllSample ECDF") + xlab("allele frequency") + ylab("") +
		scale_x_continuous(breaks=seq(0, 100, by=25)) +
	  theme_bw()+
	  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
		##
	return(grid.arrange(AF_scatter,AF_ecdf,nrow=1))
}
###################################
## Main Function        ##
## options in a list    ##
##########################
option_list <- list(
    make_option(c("-i","--input-file"), help="input file(used to read the file)", default = ""),
    make_option(c("-p","--primer-file"), help="input file(used to add primer information)", default = ""),
    make_option(c("-a","--amplicon-path"), help="input file(used to make a table)", default = "./"),
    make_option(c("-o","--output-path"), help="use to save the result", default = "./")
)

## Get command line options
      arguments <- parse_args(OptionParser(usage = "%prog [options] ", option_list = option_list), positional_arguments = 0)
      opt <- arguments$options

        ##
      Input <- opt$`input-file`
      Primer <- opt$`primer-file`
	    amp_path <- opt$`amplicon-path`
      Output <- opt$`output-path`
        ## Amp_path <- opt$`amp_path`
       
	## save data:allsites and amplicon
      m <- Transtable((Input), Primer)
	    n <- Amp_Table(amp_path, Primer)
      save(m, file=paste0(Output,"/Site.R"))
	    save(n, file=paste0(Output,"/Amplicon.R"))
	
	### draw picture of Allele frequency
     pdf(paste0(Output,"/AF.pdf"), width = 12, height = 8)
	   AF_Draw(m)
	   dev.off()
	### draw picture of AmpliconCount
	   pdf(paste0(Output,"/Amplicon.pdf"), width = 10, height = 5)
	   Amp_Draw(n)
	   dev.off()

