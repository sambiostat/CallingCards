#!/usr/bin/python

##########################################################
#Porject: Calling_Cards
#Data:2015/03/17
#Author:JINCHUN ZHANG
##########################################################

'''
Used to calculate the number of cluster and singleton of the Experimental data after QT clustered.
Input: The directory contains QT output
Output: Screen output with the number of cluster and number of singleton

'''
import os
import sys
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(add_help=True)
parser = argparse.ArgumentParser(prog='Recategory.py', description='Summary the QToutput')

parser.add_argument("--a",dest="path", type=str, nargs=1, help="The directory of QToutput")
args=parser.parse_args()
pwd=args.path[0]

if not os.path.exists(pwd):
	print "Wrong Directory, please put the one ended with folder QTout which is the output from QTcluster.py"
	sys.exit("Ending Script")

cluster=0
singleton=0
for c in range(1,25):
	inf=pwd+'chr'+str(c)+'clustered.csv'
	if os.path.isfile(inf):
		chr=pd.read_csv(inf)
		w=chr['Cluster'].tolist()
		w=np.unique(np.array(w))
		w=w[-2]
		block=chr[chr['Cluster']==99999]
		block=block[block['Independent_Insertion']>1]
		s=chr[chr['Independent_Insertion']==1]
		s=len(s)
		singleton=singleton+s
		b=len(block)
		t=b+w
		cluster=cluster+t

print 'The number of cluster is:', cluster
print 'The number of singleton is:',singleton

