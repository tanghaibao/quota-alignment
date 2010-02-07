Bounded synteny alignment
=========================

Introduction
------------

Installation
------------
Dependencies:

- python version >=2.6

- GNU linear programming kit GLPK [http://www.gnu.org/software/glpk/]

- SCIP mixed integer programming solver [http://scip.zib.de/]

Usage
-----
# Inside this folder, type `make` and then put the compiled program `bk_cliques` on your PATH.

# Look at a sample input (.dag or .cluster file), and change your file accordingly. Mostly I recommend the cluster file, with the following format:
    # cluster1
    chr1 pos1 chr2 pos2 score
    ...
    # cluster2

# Run `quota_align.py` and read the possible options.

also see `run.sh` for usage.

Reference
---------
Tang et al. Guided synteny alignment between duplicated genomes. (in preparation)
