#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Implement a few MIP solvers, based on benchmark found on <http://scip.zib.de/>
SCIP solver is ~16x faster than GLPK solver.
However, I found in rare cases it will segfault. 
Therefore the default is SCIP, the program will switch to GLPK solver for crashed cases.

The input lp_data is assumed in .lp format, see below

>>> lp_data = '''
... Maximize
...  5 x1 + 3 x2 + 2 x3
... Subject to
...  x2 + x3 <= 1
... Binary
...  x1
...  x2
...  x3
... End'''
>>> print GLPKSolver(lp_data).results
[0, 1]
>>> print SCIPSolver(lp_data).results
[0, 1]
"""

import os
import os.path as op
import sys
from subprocess import call

class AbstractMIPSolver(object):
    """
    Base class for LP solvers
    """
    def __init__(self, lp_data,
                 work_dir=op.join(op.dirname(__file__),"work"),
                 verbose=False):

        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        lpfile = work_dir + "/data.lp" # problem instance
        print >>sys.stderr, "write MIP instance to '%s'" % lpfile

        fw = file(lpfile, "w")
        fw.write(lp_data)
        fw.close()

        retcode, outfile = self.run(lpfile, work_dir, verbose=verbose)
        if retcode < 0:
            self.results = [] 
        else:
            self.results = self.parse_output(outfile)
        
        if self.results:
            print >>sys.stderr, "optimized objective value (%d)" % self.obj_val


    def run(self, lp_data, work_dir):
        pass


    def parse_output():
        pass


class GLPKSolver(AbstractMIPSolver):
    """
    GNU Linear Programming Kit (GLPK) solver, wrapper for calling GLPSOL executable
    """
    def run(self, lpfile, work_dir="work", verbose=False):

        outfile = work_dir + "/data.lp.out" # verbose output
        listfile = work_dir +"/data.lp.list" # simple output
        # cleanup in case something wrong happens
        for f in (outfile, listfile):
            if os.path.exists(f): 
                os.remove(f)

        cmd = "glpsol --cuts --fpump --lp %s -o %s -w %s"
        if not verbose: cmd += " >/dev/null"
        retcode = call(cmd % (lpfile, outfile, listfile), shell=True)

        if retcode==127:
            print >>sys.stderr, "\nError:"
            print >>sys.stderr, "You need to install program `glpsol' on your path"
            print >>sys.stderr, "[http://www.gnu.org/software/glpk/]"
            return -1, None

        return retcode, listfile


    def parse_output(self, listfile):

        filtered_list = []

        fp = file(listfile)
        header = fp.readline()
        columns, rows = header.split()
        rows = int(rows)
        data = fp.readlines()
        self.obj_val = int(data[0].split()[-1])
        # the info are contained in the last several lines
        results = [int(x) for x in data[-rows:]]
        results = [i for i, x in enumerate(results) if x==1]

        return results


class SCIPSolver(AbstractMIPSolver):
    """
    SCIP solver, wrapper for calling SCIP executable
    """
    def run(self, lpfile, work_dir="work", verbose=False):

        outfile = work_dir + "/data.lp.out" # verbose output
        if os.path.exists(outfile): 
            os.remove(outfile)

        cmd = "scip -f %s -l %s"
        if not verbose:
            cmd += " >/dev/null"

        retcode = call(cmd % (lpfile, outfile), shell=True)

        if retcode==127:
            print >>sys.stderr, "\nError:"
            print >>sys.stderr, "You need to install program `scip' on your path"
            print >>sys.stderr, "[http://scip.zib.de/]"
            return -1, None

        return retcode, outfile


    def parse_output(self, outfile):

        fp = file(outfile)
        for row in fp:
            if row.startswith("objective value"): 
                obj_row = row
                break

        results = []
        for row in fp:
            #objective value:               8
            #x1                             1   (obj:5)
            #x2                             1   (obj:3)
            if row.strip()=="": break # blank line ends the section
            x = row.split()[0]
            results.append(int(x[1:])-1) # 0-based indexing

        if results:
            self.obj_val = int(obj_row.split(":")[1])

        return results


if __name__ == '__main__':

    import doctest
    doctest.testmod()

