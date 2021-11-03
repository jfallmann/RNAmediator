============
Installation
============

For fast resolve of conda package dependencies, we recommend `mamba`_  as drop in replacement of ``conda``,
simply install mamba into your base conda environment and run commands replacing ``conda`` with ``mamba``.

.. _mamba: https://github.com/TheSnakePit/mamba


Conda
-----

For details on ``conda`` please refer to the `conda manual`_ .

.. _conda manual: https://docs.conda.io/en/latest/

Installation should be as easy as running ``conda install -c conda-forge -c bioconda rissmed``.
Best practice would be to install into an encapsulated environment via ``conda create -n rissmed -c conda-forge -c bioconda rissmed``


Pip
---

To create a working environment for this repository please install the
``rissmed.yaml`` environment as found in the ``envs`` directory
like so:

``conda env create -n rissmed -f envs/rissmed.yaml``

Followed by installation of ``RIssmed``, therefore activate the created environment

``conda activate rissmed``

and run 

``pip install rissmed`` to install the latest version


Parallelization
---------------

For distribution of jobs one can either rely on local hardware or use
scheduling software like
Slurm_ or the
SGE_.

.. _Slurm: https://slurm.schedmd.com/documentation.html
.. _SGE: https://docs.oracle.com/cd/E19957-01/820-0699/chp1-1/index.html
This manual will only show examples on local and SLURM usage.