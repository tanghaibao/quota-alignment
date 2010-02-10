#!/bin/bash

# test case 1
#cluster_utils.py --dag athaliana_alyrata.dag athaliana_alyrata.cluster 
#quota_align.py --quota 1:1 athaliana_alyrata.cluster

# test case 2
cluster_utils.py --precision 1000 athaliana_athaliana athaliana_athaliana.cluster
quota_align.py --merge --self --quota 1:1 athaliana_athaliana.cluster

# test case 3 
#cluster_utils.py --precision 1000 maize_sorghum maize_sorghum.cluster
#quota_align.py --merge --quota 2:1 maize_sorghum.cluster

# test case 4
#cluster_utils.py --precision 1000 brachy_brachy brachy_brachy.cluster
#quota_align.py --merge --self --quota 1:1 brachy_brachy.cluster

