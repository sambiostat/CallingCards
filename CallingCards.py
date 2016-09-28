#usr/bin/python


##########################################################
#Porject: Calling_Cards
#Data:2015/03/17
#Author:JINCHUN ZHANG
##########################################################
import subprocess
import argparse

parser = argparse.ArgumentParser(prog='CallingCards.py', description="Calling Cards Wrap Script, please copy the whole 'Code' into the parent directory of the folder 'Input' ")
parser.add_argument("--p",dest="path", type=str, default=["./"],nargs=1, help="Type the directory contains Input and Code,ended with /")
parser.add_argument("--c",dest="cutoff", type=str,default=["0"], nargs=1, help="Set the cutoff(Default value is 0) for the read.counts, only keep Insertion when it's count>=cutoff")
parser.add_argument("--d",dest="dist", type=str, default=["2500"], nargs=1, help="Set the maximum within cluster distance(Default value is 2500bp)")
parser.add_argument("--f",dest="pseduo", type=str, default=["0.05"], nargs=1, help="Set the pseudocount_value(Default value is 0.05) for the statistic testing")
parser.add_argument("--a",dest="Alpha", type=str, default=['0.05'], nargs=1, help="Set the alpha value (Default value is 0.05) for statistical testing, when P_value smaller than this alpha value we say its significant.") 
parser.add_argument("--e",dest="exp", type=str, nargs=1, help="Type the experimental file name,such as Sp1.Allow repeated measurement for experimental data( if repeated data called Sp1_Rep1, Sp1_Rep2, the input option can still be SP1.")
parser.add_argument("--w",dest="wt", type=str, nargs=1, help="Set the name of wild type dataset, similar as --e")
parser.add_argument("--o",dest="out", type=str, nargs=1, help="Set the name of output, like Female_Output")

args=parser.parse_args()
pwd=args.path[0]
cutoff=args.cutoff[0]
d=args.dist[0]
f=args.pseduo[0]
a=args.Alpha[0]
exp=args.exp[0]
wt=args.wt[0]
out=args.out[0]
output=str(pwd)+str(out)+'/'

print 'Calling_Cards Pipeline Started!'

subprocess.call("python "+pwd+"Code/Precutchr.py --a "+pwd+" --b " + cutoff+" --e "+exp+" --w "+wt+" --o "+out, shell=True)
subprocess.call("python "+pwd+"Code/QTcluster.py --a "+output +" --b "+ d+" --e "+exp, shell=True )
subprocess.call("python "+pwd+"Code/Combine.py --a "+output+" --e "+exp+' --w '+wt,shell=True)
subprocess.call("python "+pwd+"Code/Simpletest.py --a " +f+" --b "+a+" --c "+ output+' --w '+wt, shell=True)
subprocess.call("python "+pwd+"Code/Annotated.py --a "+pwd+" --o "+out,shell=True)
subprocess.call("mv "+output+"Annotation/Total_Insertions_P_value_Full_Annotated.txt "+output+"Annotation/"+exp+"_Insertions_Annotated.txt",shell=True)
print 'Calling_Cards Pipeline Completed!'
print 'Please check the Output Folder, the final annotated data with P_value is saved at the Annotation folder:)'
