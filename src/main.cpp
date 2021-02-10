#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "lester_mt2_bisect_v7.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

double add(double i, double j) {
    return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(mt2, m) {
    m.doc() = R"pbdoc(
        mt2 calculator module
        -----------------------
        .. currentmodule:: mt2
        .. autosummary::
           :toctree: _generate
           disableCopyrightMessage
           get_mT2
    )pbdoc";

    m.def(
        "disableCopyrightMessage",
        &asymm_mt2_lester_bisect::disableCopyrightMessage,
        py::arg("printIfFirst") = false,
        R"pbdoc(
            Wrapper around asymm_mt2_lester_bisect::disableCopyrightMessage.
        )pbdoc"
    );

    m.def(
        "get_mT2",
        py::vectorize(asymm_mt2_lester_bisect::get_mT2),
        py::arg("mVis1"),
        py::arg("pxVis1"),
        py::arg("pyVis1"),
        py::arg("mVis2"),
        py::arg("pxVis2"),
        py::arg("pyVis2"),
        py::arg("pxMiss"),
        py::arg("pyMiss"),
        py::arg("mInvis1"),
        py::arg("mInvis2"),
        py::arg("desiredPrecisionOnMT2") = 0.0,
        py::arg("useDeciSectionsInitially") = true,
        R"pbdoc(
            Wrapper around asymm_mt2_lester_bisect::get_mT2.
        )pbdoc"
    );


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
