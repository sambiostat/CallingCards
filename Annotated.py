#usr/bin/python
##########################################################
#Porject: Calling_Cards
#Data:2015/03/17
#Author:JINCHUN ZHANG
##########################################################
'''
Used for annotated Significant output
Input: the Significant insertion file, output of the Simpleteste.py
	the refGene data with bed format which can be found at the Code file.
Ouput: the txt file annotated clustered file.
	header:
           Chromosome,
	   Start, 
	   Stop, 
	   BkgHop,
	   BkgHop_withBonus, 
	   ExpHop, 
	   ExpHop_withBonus, 
	   BkgFraction, 
	   ExpFraction, 
	   P_Hyper,
	   P_Poisson,     
           CUTStart,
           CUTEnd,
	   Closest_Upstream_Gene
	   Closest_Upstream_Gene_Common_Name
           CUTstrand,
           CDTStart,
           CDTEnd,
	   Closest_Downstream_Gene
	   Closest_Downstream_Gene_Common_Name	
           CDTstrand
'''


import argparse
import subprocess
import os

parser = argparse.ArgumentParser(prog='Annotated.py', description='Find the closest down-stream and up-stream gene for significant insertions')
parser.add_argument("--a",dest="path", type=str, nargs=1,default=".", help="Type the parent directory of Input and Output")
parser.add_argument("--o",dest="out", type=str, nargs=1,default=".", help="Type the name of Output folder")

args=parser.parse_args()
pwd=args.path[0]
out=args.out[0]
Local=pwd+out+'/Total_Insertions_P_value.txt'
Ref=pwd+'Code/refGene.Sorted.bed'
outfolder=pwd+out

print "Annotated.py Starts!"

if not os.path.exists(outfolder+'/Annotation'):
    os.makedirs(outfolder+'/Annotation')
subprocess.call("sort -k1,1 -k2,2n -k3,3n "+ Local+" | sed '$d' | awk '{print $1,$2,$3}' OFS='\t' > "+outfolder+"/Annotation/Total_Insertions.sorted.bed", shell=True)
subprocess.call("sort -k1,1 -k2,2n -k3,3n "+ Local+" | sed '$d' > "+outfolder+"/Annotation/Total_Insertions.Sorted.txt", shell=True)
subprocess.call("closest-features "+outfolder+"/Annotation/Total_Insertions.sorted.bed " +Ref+" > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated.txt" ,shell=True)
subprocess.call("sed -e 's:NA:NA\tNA\tNA\tNA\tNA\tNA\t:g' -e 's/|/\t/g' "+outfolder+"/Annotation/Total_Insertions.sorted.annotated.txt > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated2.txt", shell=True)
subprocess.call("cat "+outfolder+"/Annotation/Total_Insertions.sorted.annotated2.txt |awk '{print $5,$6,$7,$8,$9,$11,$12,$13,$14,$15}' OFS='\t' > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated3.txt",shell=True)
subprocess.call("paste "+outfolder+"/Annotation/Total_Insertions.Sorted.txt "+outfolder+"/Annotation/Total_Insertions.sorted.annotated3.txt > "+outfolder+"/Annotation/Total_Insertions_P_value_Full_Annotated.txt",shell=True )
subprocess.call("sed -i '1iChromosome	Start	End	BkgHop	BkgHop_withBonus	ExpHop	ExpHop_withBonus	BkgFraction	ExpFraction	P_Hyper	P_Poisson	CUTStart	CUTEnd	Closest_Upstream_Gene	Closest_Upstream_Gene_Common_Name	CUTstrand	CDTStart	CDTEnd	Closest_Downstream_Gene	Closest_Downstream_Gene_Common_Name	CDTstrand' "+outfolder+"/Annotation/Total_Insertions_P_value_Full_Annotated.txt",shell=True )
subprocess.call("rm "+outfolder+"/Annotation/Total_Insertions.sorted.annotated*.txt "+outfolder+"/Annotation/Total_Insertions.sorted.bed", shell=True)

print "Annotated.py Completed!"

