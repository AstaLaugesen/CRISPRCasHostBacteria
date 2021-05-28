#!/usr/bin/env python3

import sys,re


#Function to print template of all clusters
seenCluster=[]
for line in sys.stdin:
	clusterName=re.search(r'\S+_(\d+)',line).group(1)
	if clusterName not in seenCluster:
		seenCluster.append(clusterName)

#Print frst col.
print("SampleName:", end="\t")
#Print the rest
for i in range(len(seenCluster)):
	print("Cluster_"+seenCluster[i]+"_w_CRISPR\t"+"Cluster_"+seenCluster[i]+"_w/o_CRISPR", end="\t")
