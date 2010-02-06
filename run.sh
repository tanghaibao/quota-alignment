#!/bin/bash

# test case 1
#cluster_utils.py --dag athaliana_alyrata.dag athaliana_alyrata.cluster 
#quota_align.py --quota 1:1 athaliana_alyrata.cluster

# test case 2
cluster_utils.py --precision 1000 athaliana_grape athaliana_grape.cluster
quota_align.py --quota 4:1 athaliana_grape.cluster
