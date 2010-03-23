#!/bin/bash

echo "select the test from below"
select TEST in `seq 1 5`;
do
case $TEST in
1)
#cluster_utils.py --format=dag --log_evalue data/athaliana_lyrata.dag data/athaliana_lyrata.qa
quota_align.py --merge --quota=1:1 data/athaliana_lyrata.qa
cluster_utils.py --print_grimm data/athaliana_lyrata.qa.filtered 
;;

2)
cluster_utils.py --format=raw --scale=1000 data/maize_sorghum.raw data/maize_sorghum.qa
quota_align.py --merge --Dm=20 --min_size 5 --quota=2:1 data/maize_sorghum.qa
;;

3)
cluster_utils.py --format=raw --scale=1000 data/grape_grape.raw data/grape_grape.qa
quota_align.py --merge --Dm=20 --self --quota=2:2 data/grape_grape.qa
;;

4)
cluster_utils.py --format=raw --scale=1000 data/brachy_brachy.raw data/brachy_brachy.qa
quota_align.py --merge --Dm=20 --self --quota=1:1 data/brachy_brachy.qa
;;

5)
cluster_utils.py --format=maf data/al_scaffold_1_vs_at_chr_1.maf data/al_scaffold_1_vs_at_chr_1.qa
quota_align.py --merge --Dm=20000 --quota=1:1 --Nm=40000 data/al_scaffold_1_vs_at_chr_1.qa
maf_utils.py data/al_scaffold_1_vs_at_chr_1.qa.filtered data/al_scaffold_1_vs_at_chr_1.maf
;;

esac
done
