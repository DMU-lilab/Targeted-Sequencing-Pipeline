#!/usr/bin/python2.7
## Date: 2017-10-10
## JiChao Tian
## Inputfile: VCF with VarScan formation including all sites; Format like that
######################################################################################################################
#Chrom   Position        Ref     Var     Cons:Cov:Reads1:Reads2:Freq:P-value     StrandFilter:R1+:R1-:R2+:R2-:pval   #
#chr1    27207291        A       .       A:1023:1023:0:0%:1E0    Pass:0:0:0:0:1E0        1       0       0       0   #
#chr1    27207292        A       .       A:1021:1021:0:0%:1E0    Pass:0:0:0:0:1E0        1       0       0       0   #
#chr1    27207293        G       .       G:1025:1025:0:0%:1E0    Pass:0:0:0:0:1E0        1       0       0       0   #
#chr1    27207294        A       .       A:1025:1025:0:0%:1E0    Pass:0:0:0:0:1E0        1       0       0       0   #
#chr1    27207295        T       .       T:1026:1026:0:0%:1E0    Pass:0:0:0:0:1E0        1       0       0       0   #
#chr1    27207296        G       .       G:1025:1025:0:0%:1E0    Pass:0:0:0:0:1E0        1       0       0       0   #
######################################################################################################################
## 
## This is for making one table including allsites  from  mutiple samples 
##
## Usage: python Dep_AF.py file1 file2 ....fileN > allsites.txt 
##
##    or  python Dep_AF.py * > allsites.txt 

import sys

def main():
    SNPdict, argc, argv = {}, len(sys.argv), sys.argv
    for index in range(1, argc):
        fp = open(argv[index], 'r')
        line = fp.readline()
        for line in fp.readlines():
            llist = line.split()
            conlist = llist[4].split(':')
            key, cons = (llist[0], int(llist[1])), [conlist[1], conlist[4]]
            if not conlist[2].isdigit():
                continue
            if key not in SNPdict:
                SNPdict[key] = [[] for i in xrange(argc-1)]
            SNPdict[key][index-1] = cons
        fp.close()
    ##output hander:" Chrom Pos Sample1_Depth Sample1_Freq Sample2_Depth Sample2_Freq ........."
    sys.stdout.write("Chrom\tPos\t")
    for i in range(1, argc):
	if i == argc-1:
		sys.stdout.write("%s_Depth\t%s_Freq\n" %(argv[i],argv[i]))
	else:
        	sys.stdout.write("%s_Depth\t%s_Freq\t" %(argv[i],argv[i]))

    for k in sorted(SNPdict.keys()):
        sys.stdout.write("%s\t%s\t" %(k[0],k[1]))
        for i in range(argc-1):
            m = SNPdict[k][i]
	    if i == argc -2:
                if not m:
                    sys.stdout.write("NM\tNM\n")
                else:
                    sys.stdout.write("%s\t%s\n" %(m[0],m[1]))
            else:
                if not m:
                    sys.stdout.write("NM\tNM\t")
                else:
                    sys.stdout.write("%s\t%s\t" %(m[0],m[1]))



if __name__ == '__main__':
    main()

