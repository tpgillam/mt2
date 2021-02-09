#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "lester_mt2_bisect_v4.h"

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

    m.def("get_mT2", py::vectorize(asymm_mt2_lester_bisect::get_mT2), R"pbdoc(
        Wrapper around asymm_mt2_lester_bisect::get_mT2.
    )pbdoc");


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
