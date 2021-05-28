#!/usr/bin/env python3

import sys,re

### Open file of crispr systems found by CRISPRCasTyper
crisprList=[]
crisprFile = sys.stdin
for crisprLine in crisprFile:
	crisprList.append(crisprLine)


### Open cluster-contig list to create list of clusters
clusterFile=open(sys.argv[2],'r')

### Making function to make list of clusters in tsv file
seenCluster=[]
#Making a function that returns the name of the clusters found
def findClusters(filename):
	seenCluster=[]
	#Opening tsv file with clusters in
	with open(filename,'r') as infile:
		for clusterLine in infile:
			#Add cluster name to list 
			clusterName=re.search(r'^\S+_(\d+)\s',clusterLine).group(1)
			if clusterName not in seenCluster:
				seenCluster.append(clusterName)
		return(seenCluster)
#Using the function
clusterList = findClusters(sys.argv[2])


#Setting name/ID of sample
sampleName=0
sampleName=re.search(r'results/(\S+)',sys.argv[1]).group(1)

#Making NA list
sampleRow=[]
index=0


### Making the main function
def searchClusterFile(clusterIndexNumber, clusterFilename):
	with open(clusterFilename,'r') as infile:
		#Make found result X and Y to be able to check any results after (if there are X and Ys then we messed up)
		clusterResult=["X","Y"]
		#for every line in the clusters file
		for line in infile:
            #if the cluster number in the line is the same as the cluster number we're looking into
			if re.search(r'^\S+_(\d+)\s',line).group(1)==clusterList[clusterIndexNumber]:
				#if cluster not within the sample
				if sampleName!=re.search(r'^(\S+)_\d+',line).group(1):
					#if we have not already changed X and Y to NA then we do that
					if clusterResult==["X","Y"]:
						clusterResult=["NA","NA"]
				
				#otherwise look through CRISPRCasTyper list
				else:
					for c in range(len(crisprList)):
						contig=re.search(r'(NODE\S+)', line).group(1)
						#if cluster is in sample and there is a CRISPRCas system to find
						if contig in crisprList[c]:
							clusterResult[0]="1"
						#if cluster does not have CRISPRCas:
						else:
							#if we have not already made "with CRISPR" as 1
							#(In case no CRISPR is to find in cluster)
							if clusterResult[0]!="1":
								clusterResult[1]="1"
								clusterResult[0]="0"
							#otherwise we have "with crispr" a 0
							#(because there are no crisprs to find)
							else:
								clusterResult[1]="0"
								
		#Append result for that cluster to a list
		sampleRow.extend(clusterResult)	

### Using the function for looking for each type
for i in range(len(clusterList)):
	searchClusterFile(i, sys.argv[2])
	i+=1

### Print result
print("S"+sampleName+"\t"+'\t'.join(sampleRow))
