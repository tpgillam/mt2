=======
History
=======

1.2.2 (2024-08-28)
------------------

* Bugfix release -- wheel were not built correctly in v1.2.1

1.2.1 (2024-08-05)
------------------

* Move support to Python 3.9-3.12. Support numpy 2. Thanks to @lgray
* Various build system & package modernisation.

1.2.0 (2021-05-05)
------------------

* Add `mt2_arxiv` public API for using (slower) implementation as published on arXiv, which may be wanted for independent testing.

1.1.0 (2021-04-18)
------------------

* Re-implementation of the core algorithm, giving 3x-4x speedup overall.

1.0.0 (2021-02-14)
------------------

* First non-beta release
* Minor formatting alterations
* Improved README

0.2.0 (2021-02-13)
------------------

* Fix numpy build dependencies to use oldest possible version, so that wheels are maximally compatible.
* Don't support Python 3.5

0.1.6 (2021-02-11)
------------------

* Fix build for numpy versions prior to 1.19

0.1.5 (2021-02-10)
------------------

* Migrate to preliminary v7 of C file
* Disable copyright message printing on run
* Move away from pybind11 in favour of raw numpy C API

0.1.4 (2021-02-09)
------------------

* Actually fix the build, maybe

0.1.3 (2021-02-09)
------------------

* Fix the build

0.1.2 (2021-02-09)
------------------

* Attempt C++ build

0.1.1 (2021-02-09)
------------------

* Trivial change to test updating package on PyPI.

0.1.0 (2021-02-09)
------------------

* First release on PyPI.
