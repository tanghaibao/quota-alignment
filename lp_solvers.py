#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Implement a few MIP solvers 

Example:

>>> lp_data = '''
... Maximize
...  5 x(1) + 3 x(2) + 2 x(3)
... Subject to
...  x(2) + x(3) <= 1
... Binary
...  x(1)
...  x(2)
...  x(3)
... End'''
>>> print GLPKSolver(lp_data).results
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
        filtered_list = [int(x) for x in data[-rows:]]
        filtered_list = [i for i, x in enumerate(filtered_list) if x==1]

        return filtered_list



class SCIPSolver(AbstractMIPSolver):
    
    def run():
        pass

    def parse_output():
        pass


if __name__ == '__main__':
    import doctest
    doctest.testmod()

