===
mt2
===


.. image:: https://img.shields.io/pypi/v/mt2.svg
        :target: https://pypi.python.org/pypi/mt2

.. image:: https://travis-ci.org/tpgillam/mt2.svg?branch=master
        :target: https://travis-ci.org/github/tpgillam/mt2


Stransverse mass computation as a numpy ufunc.


* Free software: MIT license
* Documentation: https://mt2.readthedocs.io.


Features
--------

Example usage ::

    import mt2
  
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

    val = mt2.get_mT2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB)
    
    print("Expected mT2 = 412.628.  Computed mT2 = "+str(val))
    

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
