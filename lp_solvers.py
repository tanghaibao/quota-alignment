#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Implement a few MIP solvers,
based on benchmark @ http://scip.zib.de/
SCIP solver is faster than GLPK solver

Example:

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
import sys
from subprocess import Popen

class AbstractMIPSolver(object):

    # Base class
    def __init__(self, lp_data, work_dir="work"):

        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        lpfile = work_dir + "/data.lp" # problem instance
        print >>sys.stderr, "Write problem spec to ", lpfile

        fw = file(lpfile, "w")
        fw.write(lp_data)
        fw.close()

        outfile = self.run(lpfile, work_dir)
        self.results = self.parse_output(outfile)


    def run(self, lp_data, work_dir):
        pass


    def parse_output():
        pass



class GLPKSolver(AbstractMIPSolver):

    def run(self, lpfile, work_dir="work"):

        outfile = work_dir + "/data.out" # verbose output
        listfile = work_dir +"/data.list" # simple output
        os.remove(outfile)
        os.remove(listfile)

        try:
            proc = Popen("glpsol --cuts --fpump --lp %s -o %s -w %s" % \
                    (lpfile, outfile, listfile), shell=True)
        except OSError as detail:
            print >>sys.stderr, "Error:", detail
            print >>sys.stderr, "You need to install program `glpsol' on your path"
            print >>sys.stderr, "[http://www.gnu.org/software/glpk/]"
            sys.exit(1)

        proc.communicate()

        return listfile


    def parse_output(self, listfile):
        filtered_list = []

        fp = file(listfile)
        header = fp.readline()
        columns, rows = header.split()
        rows = int(rows)
        data = fp.readlines()
        # the info are contained in the last several lines
        results = [int(x) for x in data[-rows:]]
        results = [i for i, x in enumerate(results) if x==1]

        return results



class SCIPSolver(AbstractMIPSolver):
    
    def run(self, lpfile, work_dir="work"):

        outfile = work_dir + "/data.out" # verbose output
        os.remove(outfile)

        try:
            proc = Popen("scip -f %s -l %s" % \
                    (lpfile, outfile), shell=True)
        except OSError as detail:
            print >>sys.stderr, "Error:", detail
            print >>sys.stderr, "You need to install program `scip' on your path"
            print >>sys.stderr, "[http://scip.zib.de/]"
            sys.exit(1)

        proc.communicate()

        return outfile


    def parse_output(self, outfile):

        fp = file(outfile)
        for row in fp:
            if row.startswith("objective value"): break

        results = []
        for row in fp:
            #objective value:               8
            #x1                             1   (obj:5)
            #x2                             1   (obj:3)
            if row.strip()=="": break # blank line ends the section
            x = row.split()[0]
            results.append(int(x[1:])-1)

        return results


if __name__ == '__main__':
    import doctest
    doctest.testmod()

