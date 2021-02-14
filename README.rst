===
mt2
===

.. image:: https://img.shields.io/pypi/v/mt2.svg
        :target: https://pypi.python.org/pypi/mt2

.. image:: https://github.com/tpgillam/mt2/workflows/Build/badge.svg?branch=master
        :target: https://github.com/tpgillam/mt2/actions?query=workflow%3ABuild

This package may be used to evaluate MT2 in all its variants.
This includes both symmetric and asymmetric MT2, and is also known as the "stransverse mass".

This package provides an interface to the bisection algorithm and C++ implementation of http://arxiv.org/abs/1411.4312.
The variable MT2 itself is described `here <http://arxiv.org/abs/hep-ph/9906349>`__.
Related information may be found in papers relating to MT2 linked from `here <https://www.hep.phy.cam.ac.uk/~lester/mt2/index.html>`__.

Getting started
---------------

Install from pip:

.. code-block:: bash

    pip install mt2

This will, by default, use a precompiled wheel for your platform.
If you wish to compile from source for your platform, you can pass the ``--no-binary`` option to ``pip``.
This may give

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


Note that the main reason to use vectorisation as above is convenience.
The speed penalty for writing as a python ``for`` loop is often low, due to the amount of time spent inside the C++ MT2 calculation.

Toy MC
******

A fun example using a toy Monte-Carlo simulation can be viewed in `this notebook <https://github.com/tpgillam/mt2/blob/master/examples/mc.ipynb>`__

More information
****************

For further information, see the documentation:

.. code-block:: python

    help(mt2)

Under the hood, an implementation is provided as a `numpy ufunc <https://numpy.org/doc/stable/reference/ufuncs.html>`_.
This is exported as ``mt2_ufunc``, and may be useful if wanting to use the ``where`` argument.


License
-------

Please cite:

* http://arxiv.org/abs/hep-ph/9906349, if you use MT2 in an academic paper, and
* http://arxiv.org/abs/1411.4312 if you use this particular calculator.

All files other than ``src/lester_mt2_bisect_v7.h`` are released under the MIT license.

