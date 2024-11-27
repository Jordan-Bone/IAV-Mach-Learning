#!/bin/bash

pro=(ha m1 m2 na nep np ns1 pa pb1 pb2)

for f in ${pro[@]}; 
do 
	
	mmseqs easy-linclust "$f".fasta "$f"_cluster_75 tmp --min-seq-id 0.75 -c 0.8 --cov-mode 1
	mmseqs easy-linclust "$f".fasta "$f"_cluster_85 tmp --min-seq-id 0.85 -c 0.8 --cov-mode 1
	mmseqs easy-linclust "$f".fasta "$f"_cluster_95 tmp --min-seq-id 0.95 -c 0.8 --cov-mode 1
done