#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Implement a few clique solvers, currently the Bron-Kerbosch clique enumeration is implemented

The clq_data is formatted like following, N number of nodes, >K-size cliques to report
nodes are assumed within [0, N)
N K
edges

Example:

>>> clq_data = '''
... 5 0 
... 0 1 
... 0 2
... 0 3
... 1 2
... 2 3
... End'''
>>> print BKSolver(clq_data).results
[[0, 1, 2], [0, 2, 3], [4]]

"""

import os
import os.path as op
import sys
from subprocess import Popen, call

class AbstractCLQSolver(object):

    # Base class
    def __init__(self, clq_data,
                 work_dir=op.join(op.dirname(__file__),"work")):

        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        clqfile = work_dir + "/data.clq" # problem instance
        print >>sys.stderr, "write problem spec to ", clqfile

        fw = file(clqfile, "w")
        fw.write(clq_data)
        fw.close()

        outfile = self.run(clqfile, work_dir)
        self.results = self.parse_output(outfile)


    def run(self, lp_data, work_dir):
        pass


    def parse_output():
        pass



class BKSolver(AbstractCLQSolver):

    def run(self, clqfile, work_dir="work"):

        outfile = work_dir + "/data.clq.out" 
        if os.path.exists(outfile):
            os.remove(outfile)

        path = os.path.dirname(__file__)
        retcode = call(os.path.join(path, "bk_cliques <%s >%s") % \
                (clqfile, outfile), shell=True)

        if retcode==127:
            print >>sys.stderr, "\nError:"
            print >>sys.stderr, "You need to install program `bk_cliques' on your path"
            print >>sys.stderr, "it should come with this package, " \
                    "type `make' and then copy the binary to your path"
            sys.exit(1)

        elif retcode < 0:
            print >>sys.stderr, "\nError:"
            print >>sys.stderr, "bk_cliques terminated unexpectedly"
            sys.exit(1)

        return outfile


    def parse_output(self, outfile):

        cliques = []
        fp = file(outfile)
        for row in fp:
            clique = [int(x) for x in row.split()]
            clique.sort()
            cliques.append(clique)
        
        cliques.sort()
        return cliques



if __name__ == '__main__':
    import doctest
    doctest.testmod()

