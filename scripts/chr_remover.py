#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re

"""
this script aims to remove the specified chromosome(s) in bed and pep file.
you can specify a exact name of a chromosome or use a regular expression.
"""

def handle_bed(infile, chrname):
    genelist, out_list = [], []
    with open(infile, 'r') as f:
        for line in f:
            ls = line.split('\t')
            if re.findall(chrname, ls[0]):
                genelist.append(ls[3])
            else:
                out_list.append(line)
    with open(infile+'.'+chrname+'_removed', 'w') as f:
        f.writelines(out_list)
    return genelist

def handle_pep(infile, genelist, chrname):
    out_list = []
    with open(infile, 'r') as f:
        for line in f:
            if line[0] == '>' and line.split(' ')[0][1:] in genelist:
                flag = 0
            elif line[0] == '>':
                out_list.append(line)
                flag = 1
            elif flag:
                out_list.append(line)
    with open(infile+'.'+chrname+'_removed', 'w') as f:
        f.writelines(out_list)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bedfile', 
            help='input bed file name')
    parser.add_argument('pepfile',
            help='input pep file name')
    parser.add_argument('chrname',
            help='a regular expression which matches the chromosome you want to remove')
    args = parser.parse_args()

    gene_list = handle_bed(args.bedfile, args.chrname)
    handle_pep(args.pepfile, gene_list, args.chrname)

if __name__ == '__main__':
    main()
