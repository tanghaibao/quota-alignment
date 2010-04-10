#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import itertools
import sys
from subprocess import Popen

cmd0 = "simulation.py --inv_prob %.3f --loss_prob %.3f"
cmd1 = "~/projects/quota-alignment/scripts/blast_to_raw.py blast/abblast%s --qbed bed/abed%s --sbed bed/bbed%s >/dev/null"
cmd2 = "~/projects/quota-alignment/quota_align.py --format raw blast/abblast%s.raw --merge --Dm 20 --quota 2:3 --min_size 3"

def sh(cmd):
    print >>sys.stderr, "[CMD]", cmd
    p = Popen(cmd, shell=True)
    p.communicate()

def get_lines(fn, ignore_def=True):
    fp = file(fn)
    if ignore_def:
        lines = sum(1 for x in fp if x[0]!="#")
    return lines

fw = file("summary", "w")
fw.write("inv_prob,loss_prob,res,total,recovery_rate\n")
# fix gene loss, do inversions only
I = [i/1000. for i in xrange(0, 50, 5)]
L = [i/10. for i in xrange(10)]
for inv_prob, loss_prob in itertools.product(I, L): 
    print >>sys.stderr, inv_prob, loss_prob
    tag = ("_%.3f_%.3f" % (inv_prob, loss_prob)).replace(".", "_")
    sh(cmd0 % (inv_prob, loss_prob))
    sh(cmd1 % (tag, tag, tag))
    sh(cmd2 % tag)

    res = get_lines("blast/abblast%s.raw.filtered" % tag)
    total = get_lines("blast/abblast%s" % tag)
    fw.write("%.3f,%.3f,%d,%d,%.3f\n" % (inv_prob, loss_prob, \
            res, total, res*100./total))
    fw.flush() 
fw.close()
