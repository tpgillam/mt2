===
mt2
===

.. image:: https://img.shields.io/pypi/v/mt2.svg
        :target: https://pypi.python.org/pypi/mt2

.. image:: https://img.shields.io/pypi/pyversions/mt2.svg
        :target: https://pypi.python.org/pypi/mt2

.. image:: https://github.com/tpgillam/mt2/workflows/Build/badge.svg?branch=master
        :target: https://github.com/tpgillam/mt2/actions?query=workflow%3ABuild


This package may be used to evaluate MT2 in all its variants.
This includes both symmetric and asymmetric MT2.
MT2 is also known as the "stransverse mass".

This package provides an interface to the bisection algorithm of http://arxiv.org/abs/1411.4312, via an implementation detailed below.
The variable MT2 itself is described `here <http://arxiv.org/abs/hep-ph/9906349>`__.
Related information may be found in papers relating to MT2 linked from `here <https://www.hep.phy.cam.ac.uk/~lester/mt2/index.html>`__.

Getting started
---------------

Install from `PyPI <https://pypi.org/project/mt2/>`__ with e.g. pip:

.. code-block:: bash

    pip install mt2

One can then compute MT2 as follows; here for the "symmetric" case, where both invisible particles have the same mass:

.. code-block:: python

    from mt2 import mt2

    # The units of all quantities are the same, e.g. GeV
    val = mt2(
        100, 410, 20,  # Visible 1: mass, px, py
        150, -210, -300,  # Visible 2: mass, px, py
        -200, 280,  # Missing transverse momentum: x, y
        100, 100)  # Invisible 1 mass, invisible 2 mass
    print("Expected mT2 = 412.628.  Computed mT2 = ", val)

mt2 is also `available on conda-forge <https://github.com/conda-forge/mt2-feedstock>`__. 

Examples
--------

Vectorisation
*************

The ``mt2`` function supports broadcasting over its arguments if they are array-like.
For example, one could scan over a grid of invisible particle masses like so:

.. code-block:: python

    n1 = 20
    n2 = 20
    mass_1 = numpy.linspace(10, 200, n1).reshape((-1, 1))
    mass_2 = numpy.linspace(10, 200, n2).reshape((1, -1))

    # `val` has shape (n1, n2)
    val = mt2(
        100, 410, 20,  # Visible 1: mass, px, py
        150, -210, -300,  # Visible 2: mass, px, py
        -200, 280,  # Missing transverse momentum: x, y
        mass_1, mass_2)  # Invisible 1 mass, invisible 2 mass

Note on performance
^^^^^^^^^^^^^^^^^^^

With `full precision`, the main reason to use vectorisation as above is convenience.
The time spent in the C++ MT2 calculation is somewhat larger than the overhead introduced by a Python ``for`` loop.
Vectorising can give a runtime reduction of ⪅30% in this case.

`However`, the benefit can be more significant when using a lower precision.
This corresponds to a larger value for the ``desired_precision_on_mt2`` argument.
This is because less time is spent in C++, so proportionally the Python overhead of a ``for`` loop is more significant.

Toy MC
******

A fun example using a toy Monte-Carlo simulation can be viewed in `this notebook <https://github.com/tpgillam/mt2/blob/master/examples/mc.ipynb>`__

Other notes
-----------

For further information, see the documentation:

.. code-block:: python

    help(mt2)

Also exported is ``mt2_ufunc``.
This is the raw implementation as a `numpy ufunc <https://numpy.org/doc/stable/reference/ufuncs.html>`_.
Usage is the same as for ``mt2``, but it supports some additional arguments, like ``where``.
The reader should refer to the numpy documentation for a description of these.

Implementation
**************

The underlying implementation of the Lester-Nachman algorithm used in this package is by Rupert Tombs, found in ``src/mt2_bisect.h``.
It provides results consistent with the implementation provided with http://arxiv.org/abs/1411.4312, but is 3x to 4x faster.
Note that this does *not* implement the "deci-sectioning" described in the paper, since it is found to provide a more significant performance penalty in the majority of cases.
Our version is also scale invariant, and is suitable for large ranges of input magnitude.

The legacy implementation, as it appears on arXiv, is also wrapped and exposed as ``mt2_arxiv`` for those that wish to independently cross-check the re-implementation.
If you find any discrepancies, please file a bug report!
**We strongly encourage all users to use the primary** ``mt2`` **method, due to the higher performance and scale invariance.**

Performance
***********

The default installation method via pip uses a precompiled wheel for your platform.
If you wish to compile from source for your platform, you could instead install like so:

.. code-block:: bash

    pip install mt2 --no-binary :all:

Since this can allow use of newer compilers, and code more optimised for your architecture, this can give a `small` speedup.
On the author's computer, there was 1% runtime reduction as measured with ``examples/benchmark.py``.


License
-------

Please cite:

* http://arxiv.org/abs/hep-ph/9906349, if you use MT2 in an academic paper, and
* http://arxiv.org/abs/1411.4312 if you use this particular calculator.

All files in this repository are released under the MIT license.


Other implementations
---------------------

A list of alternative implementations of the MT2 calculation can be found here:

https://www.hep.phy.cam.ac.uk/~lester/mt2/#Alternatives

In Python, the other wrapper of the same algorithm known to the authors is by Nikolai Hartmann, here: https://gitlab.cern.ch/nihartma/pymt2


Authors
-------
* @kesterlester: Original C++ implementation of mT2.
* @rupt: Current C++ implementation used in this package.
* @tpgillam: Python packaging
