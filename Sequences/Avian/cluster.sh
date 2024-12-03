#!/bin/bash

pro=(ha m1 m2 na nep np ns1 pa pb1 pb2)
# pro=(ha na np pb2)

for f in ${pro[@]}; 
do 
	echo "$f"
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_75_c80_cov1 tmp --min-seq-id 0.75 -c 0.8 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_85_c80_cov1 tmp --min-seq-id 0.85 -c 0.8 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_90_c80_cov1 tmp --min-seq-id 0.90 -c 0.8 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_95_c80_cov1 tmp --min-seq-id 0.95 -c 0.8 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_99_c80_cov1 tmp --min-seq-id 0.99 -c 0.8 --cov-mode 1
# 
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_75_c70_cov1 tmp --min-seq-id 0.75 -c 0.7 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_85_c70_cov1 tmp --min-seq-id 0.85 -c 0.7 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_90_c70_cov1 tmp --min-seq-id 0.90 -c 0.7 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_95_c70_cov1 tmp --min-seq-id 0.95 -c 0.7 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_99_c70_cov1 tmp --min-seq-id 0.99 -c 0.7 --cov-mode 1
# 
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_75_c50_cov1 tmp --min-seq-id 0.75 -c 0.5 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_85_c50_cov1 tmp --min-seq-id 0.85 -c 0.5 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_90_c50_cov1 tmp --min-seq-id 0.90 -c 0.5 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_95_c50_cov1 tmp --min-seq-id 0.95 -c 0.5 --cov-mode 1
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_99_c50_cov1 tmp --min-seq-id 0.99 -c 0.5 --cov-mode 1
# 	
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_75_c80_cov0 tmp --min-seq-id 0.75 -c 0.8
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_85_c80_cov0 tmp --min-seq-id 0.85 -c 0.8
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_90_c80_cov0 tmp --min-seq-id 0.90 -c 0.8
 	mmseqs easy-linclust "$f".fasta "$f"_cluster_95_c80_cov0 tmp --min-seq-id 0.95 -c 0.8
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_99_c80_cov0 tmp --min-seq-id 0.99 -c 0.8
# 	
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_75_c70_cov0 tmp --min-seq-id 0.75 -c 0.7
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_85_c70_cov0 tmp --min-seq-id 0.85 -c 0.7
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_90_c70_cov0 tmp --min-seq-id 0.90 -c 0.7
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_95_c70_cov0 tmp --min-seq-id 0.95 -c 0.7
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_99_c70_cov0 tmp --min-seq-id 0.99 -c 0.7
# 
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_75_c50_cov0 tmp --min-seq-id 0.75 -c 0.5
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_85_c50_cov0 tmp --min-seq-id 0.85 -c 0.5
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_90_c50_cov0 tmp --min-seq-id 0.90 -c 0.5
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_95_c50_cov0 tmp --min-seq-id 0.95 -c 0.5
# 	mmseqs easy-linclust "$f".fasta "$f"_cluster_99_c50_cov0 tmp --min-seq-id 0.99 -c 0.5
	
done

rm *all_seqs.fasta
rm *.tsv

for h in *.fasta;        
do 
	echo $h "," >> sizes.txt
	grep -o '>' $h | wc -l >> sizes.txt;
done

for j in *rep_seq.fasta;
do cp $j ../../feats/avian/.
done

