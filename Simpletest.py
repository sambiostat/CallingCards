#!/usr/bin/python

##########################################################
#Porject: Calling_Cards
#Data:2015/03/17
#Author:JINCHUN ZHANG
##########################################################
'''
Simpletest.py is used to find subcluster has significant insertion due to SP1 transcripson.
Input: SP1 data which is the output of QTcluster.py and WT data which is the output of Precutchr.py
Output: One txt file with all Significant Clusters. 
	Header is  
		Chromsome
		Total_WT_insertion 
		Total_Exp_Insertion 
		Start   
		Stop    
		P_value_from_Poission 
		P_value_from_Hypergeomeric     
		Exphop  
		WThop
'''

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import poisson
from scipy.stats import hypergeom
from tabulate import tabulate
import argparse
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(add_help=True)
parser = argparse.ArgumentParser(prog='Simpletest.py', description='Find the Experimental cluster that have significantly more independent insertions than the Wild type ')

parser.add_argument("--a",dest="fake", type=str, default='0.05', nargs=1, help="Set the pseudocount_value")
parser.add_argument("--b",dest="Alpha", type=str, default='0.05', nargs=1, help="Set the alpha value")
parser.add_argument("--c",dest="path", type=str, nargs=1, help="The directory of OUTPUT/")
parser.add_argument("--w",dest="wt", type=str, nargs=1, help="The name of wild type data/")

args=parser.parse_args()
pseudo_count=float(args.fake[0])
alpha=float(args.Alpha[0])
pwd=args.path[0]
wt=args.wt[0]

if not os.path.exists(pwd+'QTout'):
	print "Wrong Directory, please put the one contains a folder QTout which has a file named Combine.csv."
	sys.exit("Ending Script")
print "Simpletest.py start!"


'''
   Function _TEST_: Test the Null Hypothesis that the SP1 == WT.
	#1.Poisson distribution.
	# Test statistics: Pr(X>=SP1hop|WT)
	#2.Hypergeometirc districution
        # scistat.hypergeom.cdf(x,M,n,N)
        # where x is observed number of type I events (white balls in draw) (experiment hops at locus)
        # M is total number of balls (total number of hops,wt+exp total)
        # n is total number of white balls (total number of experimental hops)
        # N is the number of balls drawn (total hops at a locus)

'''
def _TEST_ (block,exptotal,bkgtotal,pseudo_count):
	sigPoi=0
	sigHyp=0
	ind=block.index
	for i in block.index:
	        lamda=block.loc[i,'BkgHop_withBonus']
        	obs=block.loc[i,'ExpHop_withBonus']
		BkgHop=block.loc[i,'BkgHop']
        	ExpHop=block.loc[i,'ExpHop']
        	P_poisson=1-poisson.cdf(int(obs)-1,lamda+pseudo_count)
	        if int(BkgHop) < 100000000 and int(ExpHop) < 100000000:
        	        P_hyper=1-hypergeom.cdf(ExpHop-1,(bkgtotal+exptotal),exptotal, (ExpHop+BkgHop))
        	else:
                	P_hyper='***'
		block.loc[i,'BkgFraction']=float(BkgHop/bkgtotal)
		block.loc[i,'ExpFraction']=float(ExpHop/exptotal)
		block.loc[i,'P_Hyper']= P_hyper
                block.loc[i,'P_Poisson']=P_poisson
		if P_poisson < alpha and lamda*obs != 0:
		        sigPoi=sigPoi+1
	        if P_hyper < alpha:
        	        sigHyp=sigHyp+1
	return block, sigPoi, sigHyp
		

summary=[]
outfile_name=pwd+'Total_Insertions_P_value.txt'
infile_name=pwd+'QTout/'+wt+'_Combine.csv'
bkgtotal=pwd+'PrecutOut/'+wt+'/combine.txt'
bkgtotal=np.loadtxt(bkgtotal,dtype=int)
bkgtotal=sum(bkgtotal[:,3])
ALL=pd.read_csv(infile_name)
exptotal=sum(ALL['ExpHop'])
counter=0

for c in range(1,22):
	block=ALL[ALL['Chromosome']==c]
	if len(block)>0:
		counter=counter+1
		out=_TEST_(block,exptotal,bkgtotal,pseudo_count)	
		block=out[0]
		x='chr'+str(c)
		block['Chromosome']=x
		singleton=block[block['Cluster']==99999]
		singleton=len(singleton[singleton['ExpHop']==1])
		summary.append([x,(len(block)-singleton),singleton,out[1],out[2]])
		del block['Cluster']
		if counter==1:
			block.to_csv(outfile_name,sep='\t',header=True,dtype='S17',index=False)
		else:
			block.to_csv(outfile_name,sep='\t',header=False,dtype='S17',index=False,mode='a')
total=['Total',0,0,0,0]
for item in summary:
	total[1]=total[1]+item[1]
	total[2]=total[2]+item[2]
	total[3]=total[3]+item[3]
	total[4]=total[4]+item[4]

summary.append(total)
print "Simpletest Result Summary:"
print "Wild Type total Hops:", bkgtotal
print "Experiment Total Hops:", exptotal
print tabulate(summary,headers=['Chr','Cluster','Singleton','Sigificant_By_Poisson_CDF','Significant_By_Hypergeometric_CDF'],tablefmt='orgtbl')
print "Simpletest.py Completed!"

