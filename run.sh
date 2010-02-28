#!/bin/bash

# test case 1
#cluster_utils.py --dag data/athaliana_alyrata.dag data/athaliana_alyrata.qa
#quota_align.py --merge --quota 1:1 data/athaliana_alyrata.qa
#cluster_utils.py --print_grimm data/athaliana_alyrata.qa.filtered 

# test case 2
#cluster_utils.py --precision 1000 data/grape_grape data/grape_grape.qa
#quota_align.py --merge --self --quota 2:2 data/grape_grape.qa

# test case 3 
#cluster_utils.py --precision 1000 data/maize_sorghum data/maize_sorghum.qa
#quota_align.py --merge --quota 2:1 data/maize_sorghum.qa

# test case 4
#cluster_utils.py --precision 1000 data/brachy_brachy data/brachy_brachy.qa
#quota_align.py --merge --self --quota 1:1 data/brachy_brachy.qa

# test case 5
cluster_utils.py --maf data/at_chr1_vs_vv_chr1.maf data/at_chr1_vs_vv_chr1.qa
quota_align.py --merge --quota 1:1 --Nm 5000 data/at_chr1_vs_vv_chr1.qa

