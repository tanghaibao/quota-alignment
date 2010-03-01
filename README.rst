Quota synteny alignment
=========================

.. contents ::

Introduction
------------

Typically in comparative genomics, we can identify anchors, chain them into syntenic blocks and interpret these blocks as derived from a common descent. However, when comparing two genomes undergone ancient genome duplications (plant genomes in particular), we have large number of blocks that are not orthologous, but are paralogous. This has forced us sometimes to use ad-hoc rules to screen these blocks. 

This program tries to screen the clusters based on the coverage constraints enforced by the user. For example, between rice-sorghum comparison, we can enforce ``1:1`` ratio to get all the orthologous blocks; or maybe ``4:2`` to grab orthologous blocks between athaliana-poplar. But the quota has to be given by the user. The program than tries to optimize the scores of these blocks globally.

To see the algorithm in action without installation, please go to `CoGe SynMap tool <http://toxic.berkeley.edu/CoGe/SynMap.pl>`_. Select "Analysis Options", select algorithm options for "Merge Syntenic Blocks" (``quota_align.py --merge``) and/or "Syntenic Depth" (``quota_align.py --quota``).

Installation
------------

- Download the most recent codes at::

    git clone http://github.com/tanghaibao/quota-alignment.git 

Required dependencies:

- Python version >=2.6

- GNU linear programming kit `GLPK <http://www.gnu.org/software/glpk/>`_. Please put the executable ``glpsol`` on the ``PATH``::

    wget http://ftp.gnu.org/gnu/glpk/glpk-4.42.tar.gz
    tar xzf glpk-4.42.tar.gz
    cd glpk-4.42/
    ./configure
    make
    sudo make install

Optional dependencies:

- `SCIP <http://scip.zib.de/download.shtml>`_ faster integer programming solver, choose the binary (32-bit, 64-bit) that fits your machine and select the one linked with CLP (for fast speed), note that in order to run SCIP, `LAPACK <http://www.netlib.org/lapack/>`_ needs to be installed too. Please put executable ``scip`` on the ``PATH``::

    unzip scip-1.2.0.linux.x86_64.gnu.opt.clp.zip
    sudo apt-get install liblapack-dev
    sudo ldconfig
    ./scip-1.2.0.linux.x86_64.gnu.opt.clp
    sudo ln -s /usr/lib/liblapack.so{,.3}
    sudo cp scip-1.2.0.linux.x86_64.gnu.opt.clp /usr/local/bin/scip
    sudo chmod +x !$

- ``bx-python`` package, this is only required when user wants to analyze ``.maf`` formatted data::

    easy_install bx-python


Usage
-----
``quota_align.py`` works only on ``.qa`` format, but the script ``cluster_utils.py`` can convert a few formats (including ``.dag`` and ``.maf``) to ``.qa`` format. Look at sample input (in the ``data/`` folder), and change your file accordingly. Mostly I recommend the ``.qa`` format::

    ###
    chr1 pos1 chr3 pos2 score
    ###
    chr1 pos1 chr2 pos2 score
    chr1 pos1 chr2 pos2 score

Note that the symbol ``#`` separates each cluster, which contains one or more *anchor points*. Anchor points within the same cluster must be on the same chromosome pair. If you don't plan to use a chainer, you can put each anchor point in its own cluster, and then rely on ``quota-align.py --merge``.

The utility script ``cluster_utils.py`` can be used for converting the Freeling lab ``.dag`` format to the ``.qa`` format, it can also print out the block sequences for downstream `GRIMM <http://grimm.ucsd.edu/GRIMM/>`_ rearrangment analysis (use ``--print-grimm`` option).

Run ``quota_align.py`` or ``cluster_utils.py`` for all possible options. 


Cookbook
--------
The default package comes with the test data for case 1 and 2. More test data set can be downloaded `here <http://chibba.agtec.uga.edu/duplication/data/quota-align-test.tar.gz>`_. Unpack into the folder, and execute ``run.sh``, also change ``TEST`` variable in ``run.sh`` for selecting different test cases.

For finding paralogous blocks, ``--self`` option must be turned on. The reason for that is in the self-matching case, the constraints on the union of the constraints on **both** axis, rather than on each axis separately.

Reference
---------
Tang et al. Guided synteny alignment between duplicated genomes.
