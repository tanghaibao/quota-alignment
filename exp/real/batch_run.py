#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Runs the four real examples in the quota-align paper. this script will assemble a table in the paper that contains the statistics to benchmark

example_name,depth(before),depth(after),ks(before),ks(after),run_time
"""

cmd1 = "~/projects/quota-alignment/cluster_utils.py --format=dag --log_evalue data/%s.dag data/%s.qa"
cmd2 = "~/projects/quota-alignment/quota_align.py --merge --Dm=0 --quota=%s --Nm=40 data/%s.qa"

import sys
import time
import numpy as np
import scipy
from subprocess import Popen

from cluster_utils import read_clusters

examples = "athaliana_lyrata athaliana_poplar athaliana_athaliana grape_grape".split()
quota = "1:1 4:2 1:1 2:2".split()

def sh(cmd):
    print >>sys.stderr, "[CMD]", cmd
    p = Popen(cmd, shell=True)
    p.communicate()

def time_quota_align(example, q):
    sh(cmd1 % (example, example))
    start_time = time.time()
    sh(cmd2 % (q, example))
    end_time = time.time()

    return end_time - start_time 

def batch_timer(examples, quota):
    fw = file("run_time.csv", "w")
    for example, q in zip(examples, quota):
        print >>fw, "%s,%.3f" % (example, time_quota_align(example, q))

def get_ks(fn):
    fp = file(fn)
    ks_list = []
    for row in fp:
        if row[0]=="#": continue
        try: 
            ks = float(row.split()[0])
            if 0 <= ks < 3: ks_list.append(ks)
        except: 
            continue
    return scipy.mean(ks_list), scipy.std(ks_list), len(ks_list) 

def batch_ks(examples):
    fw = file("ks.csv", "w")
    for example in examples:
        print >>fw, example + ",%.3f(%.3f)-%d" % \
                get_ks("ks/%s.ks.1" % example) + \
                ",%.3f(%.3f)-%d" % get_ks("ks/%s.ks.filtered.1" % example)

def get_depth(fn, q, axis=0):
    intervals = []
    clusters = read_clusters(fn)
    intervals = []
    length = 0
    for cluster in clusters:
        interval = [x[axis][1] for x in cluster]
        intervals.append(interval)
        length = max(length, max(interval))
    length += 1
    print >>sys.stderr, length, "total"
    
    depths = np.zeros(length, np.uint8)
    for interval in intervals:
        start, stop = min(interval), max(interval)
        depths[start:stop+1] += 1

    cutoff = int(q.split(":")[abs(axis-1)])
    exceed_quota = depths[depths>cutoff].shape[0]

    anchor_num = sum(len(c) for c in clusters) 
    cluster_num = len(clusters)
    return anchor_num, cluster_num, exceed_quota, length


def get_both_depth(fn, q):
    anchor_num, cluster_num, exceed_quota_x, length_x = get_depth(fn, q, axis=0) 
    anchor_num, cluster_num, exceed_quota_y, length_y = get_depth(fn, q, axis=1) 
    return anchor_num, cluster_num, (exceed_quota_x + exceed_quota_y) * 100. / (length_x + length_y)


def batch_depth(examples, quota):
    fw = file("depth.csv", "w")
    for example, q in zip(examples, quota):
        a1, b1, c1 = get_both_depth("data/%s.qa" % example, q)
        a2, b2, c2 = get_both_depth("data/%s.qa.filtered" % example, q)
        print >>fw, "%s,%d(%d),%.1f,%d(%d),%.1f" % (example, \
                a1, b1, c1, a2, b2, c2)

if __name__ == '__main__':
    batch_timer(examples, quota)
    #batch_ks(examples)
    #batch_depth(examples, quota)

