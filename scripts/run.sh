#!/bin/bash

#awk '{print $1"\t"$3"\t"$4"\t"$2}' athaliana.genes >athaliana.bed
#awk '{print $1"\t"$3"\t"$4"\t"$2}' grape.genes >grape.bed
blast_to_raw.py ~/blast/results/athaliana_grape.blastp --Nmax=20 --qbed=athaliana.bed --sbed=grape.bed
../cluster_utils.py --format=raw athaliana_grape.filtered.raw athaliana_grape.qa
../quota_align.py --merge --Dm=20 --min_size=5 --quota=4:1 athaliana_grape.qa
