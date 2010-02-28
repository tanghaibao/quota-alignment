#!/bin/bash

TEST=3
case $TEST in
1)
cluster_utils.py --format=dag data/athaliana_alyrata.dag data/athaliana_alyrata.qa
quota_align.py --merge --quota=1:1 data/athaliana_alyrata.qa
cluster_utils.py --print_grimm data/athaliana_alyrata.qa.filtered 
;;

2)
cluster_utils.py --precision=1000 data/grape_grape data/grape_grape.qa
quota_align.py --merge --self --quota=2:2 data/grape_grape.qa
;;

3)
cluster_utils.py --format=raw --precision=1000 data/maize_sorghum data/maize_sorghum.qa
quota_align.py --merge --Dm=30 --quota=2:1 data/maize_sorghum.qa
;;

4)
cluster_utils.py --precision=1000 data/brachy_brachy data/brachy_brachy.qa
quota_align.py --merge --self --quota=1:1 data/brachy_brachy.qa
;;

5)
cluster_utils.py --format=maf data/al_scaffold_1_vs_at_chr_1.maf data/al_scaffold_1_vs_at_chr_1.qa
quota_align.py --merge --Dm=10000 --quota=1:1 --Nm=40000 data/al_scaffold_1_vs_at_chr_1.qa
maf_utils.py data/al_scaffold_1_vs_at_chr_1.qa.filtered data/al_scaffold_1_vs_at_chr_1.maf
;;

esac
