Quota synteny alignment
=========================

:Author: Haibao Tang (tanghaibao)
:Email: tanghaibao@gmail.com
:License: BSD

.. contents ::

Introduction
------------

Typically in comparative genomics, we can identify anchors, chain them into syntenic blocks and interpret these blocks as derived from a common descent. However, when comparing two genomes undergone ancient genome duplications (plant genomes in particular), we have large number of blocks that are not orthologous, but are paralogous. This has forced us sometimes to use *ad-hoc* rules to screen these blocks. So the question is: **given the expected coverage (quota) along both x- and y-axis, select a subset of the anchors with maximized total score**.

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

- `SCIP <http://scip.zib.de/download.shtml>`_ faster integer programming solver, choose the binary (32-bit, 64-bit) that fits your machine and choose the one linked with CLP (for fast speed), note that in order to run SCIP, `LAPACK <http://www.netlib.org/lapack/>`_ needs to be installed too. Please put executable ``scip`` on the ``PATH``::

    unzip scip-1.2.0.linux.x86_64.gnu.opt.clp.zip
    sudo apt-get install liblapack-dev
    sudo ldconfig
    ./scip-1.2.0.linux.x86_64.gnu.opt.clp
    sudo ln -s /usr/lib/liblapack.so{,.3}
    sudo cp scip-1.2.0.linux.x86_64.gnu.opt.clp /usr/local/bin/scip
    sudo chmod +x !$

- `bx-python <http://bitbucket.org/james_taylor/bx-python/wiki/Home>`_ package, this is only required when user wants to analyze ``.maf`` formatted data::

    easy_install bx-python


Usage
-----
``quota_align.py`` works only on ``.qa`` format, but the script ``cluster_utils.py`` can convert a few formats (including ``.raw``, ``.dag`` and ``.maf``) to ``.qa`` format. Look at sample input (in the ``data/`` folder), and change your file accordingly. Mostly I recommend the ``.qa`` format, which include a list of *anchor points*::

    ###
    chr1 pos1 chr3 pos2 score
    ###
    chr1 pos1 chr2 pos2 score
    chr1 pos1 chr2 pos2 score

Each anchor point correspond to a BLAST match. Note that the symbol ``#`` separates each cluster, which contains one or more anchor points. Anchor points within the same cluster **must** be on the same chromosome pair. If you have no idea what cluster each anchor point belongs and don't plan to use a third-party chainer, you can specify the format to be ``.raw``.

The utility script ``cluster_utils.py`` can be used for converting various formats to the ``.qa`` format, it can also print out the integer sequences (representing block ids) for downstream `GRIMM <http://grimm.ucsd.edu/GRIMM/>`_ rearrangment analysis (use ``--print_grimm`` option).

Run ``quota_align.py`` or ``cluster_utils.py`` for all possible options. 

.. image:: doc/ball1.gif

Cookbook
--------
The default package comes with the test data for case 1 and 2 in ``run.sh``. More test data set can be downloaded `here <http://chibba.agtec.uga.edu/duplication/data/quota-align-test.tar.gz>`_. Unpack into the folder, and execute ``run.sh``.

BLAST anchors chaining and quota-based screening
::::::::::::::::::::::::::::::::::::::::::::::::::::
First you need to figure out a way to convert the BLAST result into the following format (called ``.raw`` format)::

    1       6       1       4848    1e-12 
    1       7       1       4847    2e-10 
    1       8       1       4847    0 
    1       9       1       4846    3e-14 

Where the five columns correspond to chr1, pos1, chr2, pos2, and E-value. Then you can convert the format (called ``.raw`` format) to the ``.qa`` format as required::

    cluster_utils.py --format=raw --log_evalue maize_sorghum.raw maize_sorghum.qa

``--log_evalue`` changes the E-value to ``max(int(-log10(E-value)),50)`` for score, so that all the BLAST anchors score in the range 0-50.

Then we can do something like::

    quota_align.py --merge --Dm=20 --quota=2:1 maize_sorghum.qa 

``--merge`` asks for chaining, and the distance cutoff ``--Dm=20`` for extending the chain; ``--quota=2:1`` turns on the quota-based screening (and asks for two-to-one match, in this case, lineage specific WGD in maize genome, make every **2** maize region matching **1** sorghum region).

BLASTZ anchors chaining and quota-based screening
:::::::::::::::::::::::::::::::::::::::::::::::::::::
Most often you will have the ``.maf`` file. First convert it to ``.qa`` format::

    cluster_utils.py --format=maf athaliana_lyrata.maf athaliana_lyrata.qa 

Then you want to do the chaining and the screening in one step::

    quota_align.py --merge --Dm=20000 --quota=1:1 --Nm=40000 athaliana_lyrata.qa 

``--merge`` asks for chaining, and the distance cutoff ``--Dm=20000`` for extending the chain; ``--quota=1:1`` turns on the quota-based screening (and asks for one-to-one match), and the overlap cutoff ``--Nm=40000``. The reason to specify an overlap cutoff is because the quota-based screening is based on 1D block overlap. Sometimes due to the over-chaining, two blocks will only *slightly* overlap. Therefore the distance ``40000`` is how much *slight* overlap we tolerate.

Finally you can get the screened ``.maf`` file by doing::

    maf_utils athaliana_lyrata.qa athaliana_lyrata.maf

Your final screened ``.maf`` file is called ``athaliana_lyrata.maf.filtered``.

Find quota-screened paralogous blocks
:::::::::::::::::::::::::::::::::::::::::
First we need to figure out how to get the input data. See the last two sections for preparing data from BLAST and BLASTZ. Then we can do something like the following::

    cluster_utils.py --format=raw grape_grape.raw grape_grape.qa
    quota_align.py --merge --Dm=20 --self --quota=2:2 grape_grape.qa

The reason for setting up ``--quota=2:2`` is because grape has `pale-hexaploidy event <http://www.nature.com/nature/journal/v449/n7161/full/nature06148.html>`_. Therefore many regions will have 3 copies, but we need to remove the self match. Therefore we should do ``2:2`` instead. ``--self`` option must be turned on for finding paralogous blocks. The reason for that is in the self-matching case, the constraints on the union of the constraints on **both** axis, rather than on each axis separately. 

For a lineage that has tetraploidy event (genome doubling), using the example of brachypodium (which has undergone an ancient tetraploidy), we can do::

    cluster_utils.py --format=raw brachy_brachy.raw brachy_brachy.qa
    quota_align.py --merge --Dm=20 --self --quota=1:1 brachy_brachy.qa

Note in this case, ``--quota=1:1`` since we have most regions in 2 copies, but we need to ignore the self match. Therefore the rule is when searching paralogous blocks (always do ``--quota=x:x``, where ``x`` is the multiplicity minus 1).

Format block order for GRIMM analysis
:::::::::::::::::::::::::::::::::::::
This is only supported when ``--quota=1:1``. For example::

    quota_align.py --merge --quota=1:1 athaliana_lyrata.qa
    cluster_utils.py --print_grimm athaliana_lyrata.qa.filtered

The script will print this::

    >genome X
    1 2 3 4 5 6 7 8 9 10 11$
    12 13 14 15 16 17 18 19$
    20 21 22 23 24 25 26 27 28 29 30 31$
    32 33 34 35 36$
    37 38 39 40 41$
    42 43 44 45 46 47 48 49 50$
    51 52 53 54 55 56 57 58$
    59 60 61$
    62 63$
    >genome Y
    -1 2 -3 4 -6 -7 5 8 10 9 11 -14 13 -12 15 16 17 18 -19$
    37 38 24 -25 26 29 28 -30 -27 31 32 33 -34 35 36$
    -21 -20 22 23 39 40 41$
    -50 49 -48 44 46 -45 47 63 -62 -55 -54 53 -52 51$
    -42 43 56 57 -58 -59 60 -61$

This is the input format for Glenn Tesler's `GRIMM <http://grimm.ucsd.edu/GRIMM/>`_ software. You can either run it locally or on their `website <http://nbcr.sdsc.edu/GRIMM/grimm.cgi>`_.


Reference
---------
Tang et al. Guided synteny alignment between duplicated genomes through integer programming.
