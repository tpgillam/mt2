===
mt2
===

.. image:: https://img.shields.io/pypi/v/mt2.svg
        :target: https://pypi.python.org/pypi/mt2

.. image:: https://github.com/tpgillam/mt2/workflows/Build/badge.svg?branch=master
        :target: https://github.com/tpgillam/mt2/actions?query=workflow%3ABuild

MT2, Asymmetric MT2, and the Stransverse Mass
---------------------------------------------

This package may be used to evaluate MT2 in all its variants.

Specifically, it provides a `numpy ufunc <https://numpy.org/doc/stable/reference/ufuncs.html>`_ interface to the bisection algorithm and C++ implementation of http://arxiv.org/abs/1411.4312.
The variable MT2 itself is described `here <http://arxiv.org/abs/hep-ph/9906349>`__.
Related information may be found in papers relating to MT2 linked from `here <https://www.hep.phy.cam.ac.uk/~lester/mt2/index.html>`__.

Example Usage
-------------

Single calculation of MT2:

.. code-block:: python

    from mt2 import mt2

    # The units of all quantities are the same, e.g. GeV
    val = mt2(
        100, 410, 20,  # Visible 1: mass, px, py
        150, -210, -300,  # Visible 2: mass, px, py
        -200, 280,  # Missing transverse momentum: x, y
        100, 100)  # Invisible 1 mass, invisible 2 mass
    print("Expected mT2 = 412.628.  Computed mT2 = ", val)

An example using broadcasting:

.. code-block:: python

    # to appear soon ....

A more interesting example using a toy Monte-Carlo simulation can be viewed in `this notebook <https://github.com/tpgillam/mt2/blob/master/examples/mc.ipynb>`__

License
-------

Please cite:

* http://arxiv.org/abs/hep-ph/9906349, if you use MT2 in an academic paper, and
* http://arxiv.org/abs/1411.4312 if you use this particular calculator.

All files other than ``src/lester_mt2_bisect_v7.h`` are released under the MIT license.

