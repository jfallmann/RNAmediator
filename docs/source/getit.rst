Installation
============

For fast resolve of conda package dependencies,
we recommend `mamba`_  as drop in replacement of ``conda``,
simply install mamba into your base conda environment
and run commands replacing ``conda`` with ``mamba``.

.. _mamba: https://github.com/TheSnakePit/mamba


Conda
-----

For details on ``conda`` please refer to the `conda manual`_ .

.. _conda manual: https://docs.conda.io/en/latest/

Installation should be as easy as running:

.. code-block ::

    conda install -c conda-forge -c bioconda rnamediator

.. note::

    Best practice would be to install into an encapsulated environment via:

    .. code-block ::

        conda create -n rnamediator -c conda-forge -c bioconda rnamediator


Pip
---

Installation via Pip is simple. Just type:

.. code-block ::

    pip install rnamediator

.. note::

    Most of the required packages will be installed. However, the ViennaRNA_ package needs to be in PATH. Also, make sure that the Python API of the
    ViennaRNA package is installed.

    .. _ViennaRNA: https://www.tbi.univie.ac.at/RNA/


Parallelization
---------------

For distribution of jobs one can either rely on local hardware or use
scheduling software like
Slurm_ or the
SGE_.

.. _Slurm: https://slurm.schedmd.com/documentation.html
.. _SGE: https://docs.oracle.com/cd/E19957-01/820-0699/chp1-1/index.html

This manual will only show examples on local and SLURM usage.
