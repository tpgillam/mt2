===
mt2
===


.. image:: https://img.shields.io/pypi/v/mt2.svg
        :target: https://pypi.python.org/pypi/mt2

.. image:: https://travis-ci.org/tpgillam/mt2.svg?branch=master
        :target: https://travis-ci.org/github/tpgillam/mt2

MT2, Asymmetric MT2, and the Stransverse Mass
---------------------------------------------

This package may be used to evaluate MT2 in all its variants.

Specifically, it provides a numpy ufunc interface to the bisection algorithm and C++ implementation of http://arxiv.org/abs/1411.4312.
The variable MT2 itself is described in http://arxiv.org/abs/hep-ph/9906349.
Related information may be found in papers relating to MT2 linked from https://www.hep.phy.cam.ac.uk/~lester/mt2/index.html.

License
-------

Please cite:

* http://arxiv.org/abs/hep-ph/9906349, if you use MT2 in an academic paper, and
* http://arxiv.org/abs/1411.4312 if you use this particular calculator.


Features
--------

Example usage ::

    from mt2 import mt2

    pxA =   410
    pyA =    20
    mVisA = 100
    chiA =  100

    pxB =  -210
    pyB =  -300
    mVisB = 150
    chiB =  100

    pxMiss = -200
    pyMiss =  280

    val = mt2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB)

    print("Expected mT2 = 412.628.  Computed mT2 = ", val)

