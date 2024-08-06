from setuptools import Extension, setup

import numpy

setup(
    ext_modules=[
        Extension(
            "_mt2",
            ["src/main.cpp"],
            define_macros=[
                # For reasons explained in lester_mt2_bisect_v7.h, we need to manually
                # enable some inlining optimisations.
                ("ENABLE_INLINING", "1"),
                # Copyright printing is disabled here, since we include the necessary
                # citation information elsewhere.
                ("DISABLE_COPYRIGHT_PRINTING", "1"),
            ],
            include_dirs=[numpy.get_include()],
            language="c++",
            extra_compile_args=["-std=c++11"],
        ),
    ],
)
