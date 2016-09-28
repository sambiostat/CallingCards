#!/usr/bin/python

##########################################################
#Porject: Calling_Cards
#Data:2015/03/17
#Author:JINCHUN ZHANG
##########################################################

'''
Input: txt data with format'chromosome	position count'
Output: 21 txt format data for 21 chromosomes, sorted by position, 1st column is Chromosome and the 2nd column is the position, 3rd column is their number of independent insertion with 5^N bonus for N>1
*cutoff is used for filter out the independent insertions based on their own count which is the 3rd column of the input dataset
*pwd is the directory that has WT.txt and Exprimental.txt.
'''
import argparse
import numpy as np
import os
import subprocess

print "Precut.py start!"

parser = argparse.ArgumentParser(prog='Precut.py', description='Count the Independent Insertions for each position and cut the file into 21 files.')
parser.add_argument("--e",dest="exp",type=str,nargs=1,help="Type the name of experimental group,such as'sp1'. which all saved in Input Folder")
parser.add_argument("--a",dest="path", type=str, nargs=1,default=".", help="Type the directory of Input(ended with Input/)")
parser.add_argument("--b",dest="cutoff", type=int, nargs=1, default=0,help="Set the cutoff for the read.counts, only keep Insertion when it's count > cutoff")
parser.add_argument("--w",dest="wt", type=str, nargs=1, default=0,help="Set the name of wild type dataset, similar as --e")
parser.add_argument("--o",dest="out", type=str, nargs=1, default=0,help="Set the name of output, like Female_Output")

args=parser.parse_args()
pwd=args.path[0]
cutoff=args.cutoff[0]
exp=args.exp[0]
wt=args.wt[0]
out=args.out[0]
print "Cutoff for filterout reads is: ", cutoff

subprocess.call("cat "+pwd+"Input/*"+wt+"* > "+ pwd+"Input/"+wt+".txt", shell=True)
subprocess.call("cat "+pwd+"Input/*"+exp+"*> " + pwd+"Input/"+exp+".txt",shell=True)

for x in ([str(wt),str(exp)]):
	if not os.path.exists(pwd+out+'/PrecutOut/'+x):
		os.makedirs(pwd+out+'/PrecutOut/'+x)
	sdin=open(pwd+'Input/'+x+'.txt','r')
	inf=sdin.readlines()
	sdin.close()
	chr={
            1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],
            13:[],14:[],15:[],16:[],17:[],18:[],19:[],20:[],21:[]
  	  }
	for i in inf:
        	a=i.split('\t',2)
		    if len(a[0])<3:
        		a[0]=int(a[0])
			if a[0]==23:
				a[0]=20
			if a[0]==24:
				a[0]=21
        		a[1]=int(a[1])
	        	a[2]=int(a[2].replace('\n',''))
        		if a[2]>cutoff:
            			chr[a[0]].append([a[1],a[2]])
    for i in chr.keys():
		if len(chr[i])>0:
	        	chr[i]=np.array(chr[i])
        		uniq=np.unique(chr[i][:,0])
        		new=np.zeros((len(uniq),4),dtype='int32')
	        	new[:,1]=uniq
                new[:,0]=i
        		for n in range(len(uniq)):
            			freq=chr[i][chr[i][:,0]==uniq[n]]
	            		freq=freq[:,0]
				new[n,3]=len(freq)
				if len(freq)==1:
                    			new[n,2]=1
	                	else:
        	            		new[n,2]=len(freq)+5**len(freq)		
        		outf=pwd+out+'/PrecutOut/'+x+'/chr'+str(i)+'.txt'
			with open(outf,'w') as of:
				of.write('Chromosome\tPosition\tIndependent_Insertions_withBonus\tIndependent_Insertions\n')
        			np.savetxt(of, new,fmt="%d")			

subprocess.call("tail -n +2 -q "+pwd+out+"/PrecutOut/"+wt+"/chr*.txt> "+ pwd+out+"/PrecutOut/"+wt+"/combine.txt", shell=True)
print "Precutchr.py completed!"
