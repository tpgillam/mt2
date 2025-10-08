from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

import numpy

__version__ = "1.3.1"


class _CustomBuildExt(build_ext):
    def build_extensions(self) -> None:
        for extension in self.extensions:
            assert len(extension.extra_compile_args) == 0
            if self.compiler.compiler_type == "msvc":
                extension.extra_compile_args = [
                    "/std:c++11",
                    "/W4",
                    "/WX",
                ]
            else:
                # Assume a unix-like compiler otherwise.
                extension.extra_compile_args = [
                    "-std=c++11",
                    "-Wall",
                    "-pedantic",
                    "-Werror",
                ]
        super().build_extensions()


setup(
    ext_modules=[
        Extension(
            "mt2._mt2",
            ["src/_mt2/main.cpp"],
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
        ),
    ],
    cmdclass={"build_ext": _CustomBuildExt},
)
