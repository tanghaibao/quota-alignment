Quota synteny alignment
=========================

.. contents ::

Introduction
------------

Typically in comparative genomics, we can identify anchors, chain them into syntenic blocks and interpret these blocks as derived from a common descent. However, when comparing two genomes undergone ancient genome duplications (plant genomes in particular), we have large number of blocks that are not orthologous, but are paralogous. This has forced us sometimes to use ad-hoc rules to screen these blocks. 

This program tries to screen the clusters based on the coverage constraints imposed by the user. For example, between rice-sorghum comparison, we can enforce 1:1 ratio to get all the orthologous blocks; or maybe 4:2 to grab orthologous blocks between athaliana-poplar. But the quota has to be given by the user. The program than tries to optimize the scores of these blocks globally.

To see the algorithm in action without installing, please go to `CoGe SynMap tool <http://synteny.cnr.berkeley.edu/CoGe/SynMap.pl>`. Select "Analysis Options".

Installation
------------

- Download the most recent codes at::

    git clone http://github.com/tanghaibao/quota-alignment.git 

Dependencies:

- python version >=2.6

- GNU linear programming kit `GLPK <http://www.gnu.org/software/glpk/>`_::

    cd quota-alignment/
    make
    mkdir tools
    cd tools/
    wget http://ftp.gnu.org/gnu/glpk/glpk-4.42.tar.gz
    tar xzf glpk-4.42.tar.gz
    cd glpk-4.42/
    ./configure
    make
    sudo make install
    glpsol


- (*optional*) SCIP mixed integer programming solver linked with `CLP <http://scip.zib.de/download.shtml>`_, choose the binary that fits your machine::

    unzip scip-1.2.0.linux.x86_64.gnu.opt.clp.zip
    ./scip-1.2.0.linux.x86_64.gnu.opt.clp
    sudo apt-get install liblapack-dev
    sudo ldconfig
    ./scip-1.2.0.linux.x86_64.gnu.opt.clp
    sudo ln -s /usr/lib/liblapack.so{,.3}
    sudo cp scip-1.2.0.linux.x86_64.gnu.opt.clp /usr/local/bin/scip
    sudo chmod +x !$
    scip
    cd ../
    ./lp_solvers.py

- compile program ``bk_cliques`` shipped within this package::

    make
    sudo cp bk_cliques /usr/local/bin


Usage
-----
- Look at a sample input (``.dag`` or ``.cluster`` file), and change your file accordingly. Mostly I recommend the ``.cluster`` format::

    # cluster1
    chr1 pos1 chr2 pos2 score
    ...
    # cluster2

The utility script ``cluster_utils.py`` can be used for converting the Freeling lab ``.dag`` format to the ``.cluster`` format, it can also print out the block sequences for downstream `GRIMM <http://grimm.ucsd.edu/GRIMM/>`_ rearrangment analysis (use ``--print-grimm`` option).

- Run ``quota_align.py`` and read the possible options.

also see ``run.sh`` for usage.


Cookbook
--------
To be written.


Reference
---------
Tang et al. Guided synteny alignment between duplicated genomes.
