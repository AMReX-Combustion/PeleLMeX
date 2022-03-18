.. PeleLMeX documentation master file, created by
      sphinx-quickstart on Sat Oct 20 12:17:48 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the documentation for `PeleLMeX`
=========================================

`PeleLMeX` is the non-subcycling version of `PeleLM <https://amrex-combustion.github.io/PeleLM/>`_, an adaptive-mesh low Mach number hydrodynamics 
code for reacting flows. If you need help or have questions, please join the users `forum <https://groups.google.com/forum/#!forum/pelelmusers>`_.
The documentation pages appearing here are distributed with the code in the ``Docs`` folder as "restructured text" files.  The html is built
automatically with certain pushes to the `PeleLM` GibHub repository and are maintained online by `ReadTheDocs <https://pelelmex.readthedocs.io/en/latest>`_.
A local version can also be built as follows ::

    cd ${PELELMEX_HOME}/Docs
    make html

where ``PELELMEX_HOME`` is the location of your clone of the `PeleLMeX` repository.  To view the local pages,
point your web browser at the file ``${PELELMEX_HOME}/Docs/build/html/index.html``.

**Current docs build status on ReadTheDocs:**

.. image:: https://readthedocs.org/projects/pelelmex/badge


.. toctree::
   :maxdepth: 2
   :caption: Documentation contents:

   LMeXControls.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
