Quota synteny alignment
=========================

Introduction
------------

Installation
------------

- download the most recent codes at
  git clone git@github.com:tanghaibao/quota-alignment.git

Dependencies:

- python version >=2.6

- GNU linear programming kit GLPK [http://www.gnu.org/software/glpk/]::

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
  cd ../../


- (*optional*) SCIP mixed integer programming solver linked with CLP [http://scip.zib.de/download.shtml], choose your platform::

  unzip scip-1.2.0.linux.x86_64.gnu.opt.clp.zip
  ./scip-1.2.0.linux.x86_64.gnu.opt.clp
  sudo apt-get install liblapack-dev
  sudo ldconfig
  ./scip-1.2.0.linux.x86_64.gnu.opt.clp
  sudo ln -s /usr/lib/liblapack.so{,.3}
  sudo cp scip-1.2.0.linux.x86_64.gnu.opt.clp /usr/local/bin/
  scip
  sudo chmod +x !$
  cd ../
  ./lp_solvers.py

- compile program `bk_cliques` shipped within this package::
  make
  sudo cp bk_cliques /usr/local/bin


Usage
-----
- Look at a sample input (.dag or .cluster file), and change your file accordingly. Mostly I recommend the cluster file, with the following format::

    # cluster1
    chr1 pos1 chr2 pos2 score
    ...
    # cluster2

- Run `quota_align.py` and read the possible options.

also see `run.sh` for usage.

Reference
---------
Tang et al. Guided synteny alignment between duplicated genomes. (in preparation)
