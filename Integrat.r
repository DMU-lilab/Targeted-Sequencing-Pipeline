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
library(ggplot2)
library(optparse) 
library(gridExtra)


Transtable <- function(pathfile){
 ### 读入文件
  table_all <- fread(pathfile, sep="\t", header=T)
 ### 转变格式
  table_m <- melt(table_all,id=1:2) 
 ### 分割列
  list_v <- strsplit(as.character(table_m$variable), "_")
  list_v <- data.frame(do.call(rbind, list_v))
  table_m[,':='(Sample=list_v$X1, Type=list_v$X2, df=list_v$X3)]
 ### 转变格式
  table_f <- dcast(table_m, Chrom+Pos+Sample+Type~df, value.var = "value")
 ### 重新排序
  table_f <- data.table(table_f)[,c(1,2,3,5,6,4)]	
  return(table_f)
}


Data_combine <- function(table, pri_path){
 ## 读入primer文件
  primer <- fread(pri_path,sep="\t",header=T)
 ## Chr Pos Sample Depth Freq Type
 ## 先把Chrom Pos 取出，unique，比对到primer文件，生成需要的注释
  Depth_Freq_All_Samp <- table
  Chr_Pos <- unique(Depth_Freq_All_Samp[,.(Chr,Pos)])
  setkey(primer,"Chr")
  setkey(Chr_Pos,"Chr","Pos")
  dt_lst <- list()
  for (chrI in unique(primer[["Chr"]])){
    dt_pos_chr <- Chr_Pos[chrI,]
    dt_insert_chr <- primer[chrI,]
    dt_primer_chr <-  rbind(dt_insert_chr[,.(Chr,Amplicon_Start,Insert_Start-1,Ampl,Gene,PrDir="F")],dt_insert_chr[,.(Chr,Insert_Stop+1,Amplicon_Stop,Ampl,Gene,PrDir="R")],use.names=FALSE)
    setnames(dt_primer_chr,c("Chr","Start","End","Ampl","Gene","Prdir"))
    dt_primer_chr <- rbind(dt_primer_chr,data.table(chrI,0,0,"NA","NA","NA"),use.names=FALSE) #将yid中的NA，替换自定义列的列号，以便下一步的运行
    setkey(dt_primer_chr,"Chr","Start","End")
    rg_pos <- dt_pos_chr[,.(Start=Pos,End=Pos)]
    rg_insert <- dt_insert_chr[,.(Insert_Start,Insert_Stop)]
    rg_primer <- dt_primer_chr[,.(Start,End)]
    setkey(rg_pos,"Start","End")
    setkey(rg_insert,"Insert_Start","Insert_Stop")
    setnames(rg_primer,c("Start","End"))
    setkey(rg_primer,"Start","End")
    dt_pos_insert_overl <- foverlaps(rg_pos, rg_insert, which=TRUE,nomatch=0)
    dt_pos_insert <- unique(data.table(dt_pos_chr[dt_pos_insert_overl[["xid"]],],Gene=dt_insert_chr[["Gene"]][dt_pos_insert_overl[["yid"]]],Ampl1=dt_insert_chr[["Ampl"]][dt_pos_insert_overl[["yid"]]],Ampl2=as.character(NA)))
    torm <- list()
    for (i in 1:(nrow(dt_pos_insert)-1)){if(dt_pos_insert[i,"Pos"]==dt_pos_insert[i+1,"Pos"]){set(dt_pos_insert,i,ncol(dt_pos_insert),dt_pos_insert[i+1,"Ampl1"]);torm <- c(torm,list(i+1))}}
    torm <- unlist(torm)
    dt_pos_insert <- dt_pos_insert[-torm,]
    setkey(dt_pos_insert,"Chr","Pos")
    rg_pos_insert <- dt_pos_insert[,.(Start=Pos,End=Pos)]
    setkey(rg_pos_insert)
    dt_pos_insert_primer_overl <- foverlaps(rg_pos_insert, rg_primer, which=TRUE,nomatch=NA)
    dt_pos_insert_primer_overl[is.na(yid),yid:=as.integer(row.names(dt_primer_chr[Ampl=="NA"]))]
    dt_pos_insert_primer <- unique(data.table(dt_pos_insert[dt_pos_insert_primer_overl[["xid"]]],Prdir1=dt_primer_chr[["Prdir"]][dt_pos_insert_primer_overl[["yid"]]],Prdir2=as.character(NA),Primer1=dt_primer_chr[["Ampl"]][dt_pos_insert_primer_overl[["yid"]]],Primer2=as.character(NA)))
    for (i in 1:(nrow(dt_pos_insert_primer)-1)){if(dt_pos_insert_primer[i,"Pos"]==dt_pos_insert_primer[i+1,"Pos"] & dt_pos_insert_primer[i,"Primer1"] != 1){set(dt_pos_insert_primer,i,ncol(dt_pos_insert_primer),dt_pos_insert_primer[i+1,"Primer1"]);set(dt_pos_insert_primer,i,7L,dt_pos_insert_primer[i+1,"Prdir1"]);torm <- c(torm,list(i+1))}}
  # add PrStart and Pr End
    tmp1 <- merge(dt_pos_insert_primer, dt_primer_chr, all.x=T, by.x=c("Primer1","Prdir1"), by.y=c("Ampl","Prdir"), allow.cartesian=TRUE)
    setnames(tmp1,c("Start","End","Chr.x","Gene.x"),c("PrStart1","PrEnd1","Chr","Gene"))
    tmp1$Chr.y <- NULL
    tmp1$Gene.y <- NULL
    tmp2 <- merge(tmp1, dt_primer_chr, all.x=T, by.x=c("Primer2","Prdir2"), by.y=c("Ampl","Prdir"), allow.cartesian=TRUE)
    setnames(tmp2, c("Start", "End", "Chr.x", "Gene.x"), c("PrStart2", "PrEnd2", "Chr", "Gene"))
    tmp2$Chr.y <- NULL
    tmp2$Gene.y <- NULL
  # add Primer 5’ Postion
    tmp2[Prdir1=="F", Pr5Pos1:=Pos-PrStart1+1]
    tmp2[Prdir1=="R", Pr5Pos1:=PrEnd1-Pos+1]
    tmp2[Prdir2=="F", Pr5Pos2:=Pos-PrStart2+1]
    tmp2[Prdir2=="R", Pr5Pos2:=PrEnd2-Pos+1]
  
    dt_pos_insert_primer <- tmp2[,.(Chr=Chr, Pos=Pos, Gene=Gene, Ampl1=Ampl1, Ampl2=Ampl2, Primer1=Primer1, Prdir1=Prdir1, PrStart1=PrStart1, PrEnd1=PrEnd1, Pr5Pos1=Pr5Pos1, Primer2=Primer2, Prdir2=Prdir2, PrStart2=PrStart2, PrEnd2=PrEnd2, Pr5Pos2=Pr5Pos2)]
    dt_lst[[chrI]] <- dt_pos_insert_primer
  }
  dt_all_pos <- rbindlist(dt_lst)
  #再将添加Ampl和Primer的Chr,Pos merge回去原来的大数据框
  Depth_Freq_All_Samp_Primer <- merge(Depth_Freq_All_Samp, dt_all_pos,by=c("Chr","Pos"), all.y=TRUE)
  return(Depth_Freq_All_Samp_Primer)
}


#################################
## Amplicon Information  ##
## Creat a summary table ##
###########################	
Amp_Table <- function(amp_path){
        Amp_count <- list.file(amp_path)
        Amp <- data.frame()
     ## 整合所有Amplicon文件
        for ( i in Amp_count){
            amp_e <- fread(amp_path"/"i, sep="\t", header=T)
            amp_e$Sample <- sub(".count", "", i) 
            Amp <- rbind(Amp, amp_e)
        }
     ## sp  <- strsplit(Amp$AuxInfo, "_")
     ## sp  <- data.frame(do.call(rbind,sp))	
        Amall <- data.frame(do.call(rbind, strsplit(Amp$AuxInfo,",")))
        colnames(Amall) <- c("Gene", "Chr", "FP", "IS", "IE", "RP")
        Amp <- cbind(Amp, Amall)
        Amp[, AuxInfo := NULL]
        Amp[,':='(Log2=round(log(AmpCount+1,2),2), Log10=round(log(AmpCount+1,10),2))]
        return(Amp)
} 
 
############################################################
## For drawing some pictures                     ##
## Amplicon pictures :Boxplot, Violin, Histogram ##
###################################################
Amp_Draw <- function(AMP){
	## Boxplot
	Amp_boxplot <- ggplot(AMP,aes(Type,Log2,fill=Type)) + geom_boxplot() +
		#scale_x_discrete(limits=c("WBC","cfDNA1","cfDNA2","cfDNA3","tissue")) + 
		ggtitle("Amplicon Depth") + xlab("x-axis") + ylab("y-axis") +
		guides(fill=guide_legend(title="Legend_Title")) + ylim(5,15)
	## Violoin
	Amp_violin <- ggplot(AMP,aes(Type,log2,fill=Type)) + geom_violin() +
		#scale_x_discrete(limits=c("WBC","cfDNA1","cfDNA2","cfDNA3","tissue")) +
		ggtitle("Amplicon Depth") + xlab("Type") + ylab("Amplicon Depth") +
		guides(fill=guide_legend(title="Legend_Title")) +
		theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
	## Histogram
	Amp_histogram <- ggplot(AMP,aes(log2,fill=Type)) + geom_histogram(bins=90) +
		ggtitle("Amplicon Depth Distribution") + xlab("AmpliconDepth") + ylab("Count") + 
		guides(fill=guide_legend(title="Legend_Title")) + 
		theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) + 
		facet_grid(~Type)
	return(grid.arrange(Amp_boxplot,Amp_violin,Amp_histogram,nrow=2,ncol=2))
}

####################################################
### Allele Fraction : Scatter plot; ECDF ##
###########################################
AF_Draw <- function(AF){
	## Scatterplot
	AF_scatter <- ggplot(AF, aes(Type, Freq, colour=Type)) + geom_jitter() +
		ggtitle("AllSample Allele Frequency") + xlab("Type") + ylab("Allele Frequency") +
		#scale_color_manual(values = c("red","blue","orange","green","black"))+
		theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
	## ECDF plot
	AF_ecdf <- ggplot(AF,aes(Freq,colour=Type)) + stat_ecdf() + ylim(0.999,1) +
		facet_wrap(~Classification) +
		ggtitle("AllSample ECDF") + xlab("allele frequency") + ylab("") +
		scale_x_continuous(breaks=seq(0, 100, by=25)) + 
		#scale_color_manual(values = c("red","blue","orange","green","black"))
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
        m <- Data_combine(Transtable(Input), Primer)
	n <- Amp_Table(amp_path)
        save(m, file=Output"/Site.R")
	save(n, file=Output"/Amplicon.R")
	
	### draw picture
        pdf(Output"/AF.pdf"), width = 10, height = 5)
	AF_Draw(m)
	dev.off()
	###
	pdf(Output"/Amplicon.pdf"), width = 10, height = 5)
	Amp_Draw(n)
	dev.off()
