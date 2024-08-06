from setuptools import Extension, setup

import numpy

__version__ = "1.2.1"

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

ext_modules = [
    Extension(
        "_mt2",
        ["src/main.cpp"],
        define_macros=[
            # Pass in the version info so we can expose it in the extension.
            ("VERSION_INFO", __version__),
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
]


setup(
    author="Thomas Gillam",
    author_email="tpgillam@googlemail.com",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    description="Stransverse mass computation as a numpy ufunc.",
    install_requires=["numpy>=1.19.3"],
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="mt2",
    name="mt2",
    packages=["mt2"],
    test_suite="tests",
    tests_require=["pytest>=3"],
    url="https://github.com/tpgillam/mt2",
    version=__version__,
    ext_modules=ext_modules,
    zip_safe=False,
)
