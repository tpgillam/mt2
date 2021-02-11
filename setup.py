from setuptools import setup, Extension

import numpy

__version__ = '0.1.5'

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

ext_modules = [
    Extension(
        "_mt2",
        ["src/main.cpp"],
        # Pass in the version info so we can expose it in the extension.
        # Also, for reasons explained in lester_mt2_bisect_v7.h, we need to manually enable some inlining optimisations.
        # Copyright printing is disabled here, since we include the necessary citation information elsewhere.
        define_macros=[
            ('VERSION_INFO', __version__),
            ('ENABLE_INLINING', 1),
            ('DISABLE_COPYRIGHT_PRINTING', 1),
        ],
        include_dirs=[numpy.get_include()],
        language='c++'
    ),
]

setup_requirements = ['pytest-runner', ]
test_requirements = ['pytest>=3', ]

setup(
    author="Thomas Gillam",
    author_email='tpgillam@googlemail.com',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Stransverse mass computation as a numpy ufunc.",
    install_requires=["numpy"],
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='mt2',
    name='mt2',
    packages=['mt2'],
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/tpgillam/mt2',
    version=__version__,
    ext_modules=ext_modules,
    zip_safe=False,
)
