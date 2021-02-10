#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>

#include "lester_mt2_bisect_v7.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

static void mt2_ufunc(
    char** args,
    npy_intp const* dimensions,
    npy_intp const* steps,
    void* data)
{
    npy_intp n = dimensions[0];

    char* mVis1 = args[0];
    char* pxVis1 = args[1];
    char* pyVis1 = args[2];
    char* mVis2 = args[3];
    char* pxVis2 = args[4];
    char* pyVis2 = args[5];
    char* pxMiss = args[6];
    char* pyMiss = args[7];
    char* mInvis1 = args[8];
    char* mInvis2 = args[9];
    char* desiredPrecisionOnMT2 = args[10];
    char* useDeciSectionsInitially = args[11];
    char* out = args[12];

    npy_intp mVis1_step = steps[0];
    npy_intp pxVis1_step = steps[1];
    npy_intp pyVis1_step = steps[2];
    npy_intp mVis2_step = steps[3];
    npy_intp pxVis2_step = steps[4];
    npy_intp pyVis2_step = steps[5];
    npy_intp pxMiss_step = steps[6];
    npy_intp pyMiss_step = steps[7];
    npy_intp mInvis1_step = steps[8];
    npy_intp mInvis2_step = steps[9];
    npy_intp desiredPrecisionOnMT2_step = steps[10];
    npy_intp useDeciSectionsInitially_step = steps[10];
    npy_intp out_step = steps[12];

    for (npy_intp i = 0; i < n; ++i) {
        *((double *)out) = asymm_mt2_lester_bisect::get_mT2(
            *(double *)mVis1,
            *(double *)pxVis1,
            *(double *)pyVis1,
            *(double *)mVis2,
            *(double *)pxVis2,
            *(double *)pyVis2,
            *(double *)pxMiss,
            *(double *)pyMiss,
            *(double *)mInvis1,
            *(double *)mInvis2,
            *(double *)desiredPrecisionOnMT2,
            *(npy_bool *)useDeciSectionsInitially
        );

        mVis1 += mVis1_step;
        pxVis1 += pxVis1_step;
        pyVis1 += pyVis1_step;
        mVis2 += mVis2_step;
        pxVis2 += pxVis2_step;
        pyVis2 += pyVis2_step;
        pxMiss += pxMiss_step;
        pyMiss += pyMiss_step;
        mInvis1 += mInvis1_step;
        mInvis2 += mInvis2_step;
        desiredPrecisionOnMT2 += desiredPrecisionOnMT2_step;
        useDeciSectionsInitially += useDeciSectionsInitially_step;
        out += out_step;
    }
}

/* This a pointer to the above function */
PyUFuncGenericFunction ufuncs[1] = {&mt2_ufunc};

/* These are the input and return dtypes of mt2_ufunc.*/
static char types[13] = {
    NPY_DOUBLE, // double mVis1,
    NPY_DOUBLE, // double pxVis1,
    NPY_DOUBLE, // double pyVis1,
    NPY_DOUBLE, // double mVis2,
    NPY_DOUBLE, // double pxVis2,
    NPY_DOUBLE, // double pyVis2,
    NPY_DOUBLE, // double pxMiss,
    NPY_DOUBLE, // double pyMiss,
    NPY_DOUBLE, // double mInvis1,
    NPY_DOUBLE, // double mInvis2,
    NPY_DOUBLE, // double desiredPrecisionOnMT2 = 0
    NPY_BOOL,  // bool useDeciSectionsInitially=true
    NPY_DOUBLE  // <result>
};


PyDoc_STRVAR(mt2_module_doc, "Provides the mt2 stransverse mass ufunc.");

static PyMethodDef methods[] = {
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_mt2",
    mt2_module_doc,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};

static void *data[1] = {NULL};

PyMODINIT_FUNC PyInit__mt2(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    if (!module) {
        return NULL;
    }

    import_array();
    import_ufunc();
    import_umath();

    PyObject* mt2_ufunc = PyUFunc_FromFuncAndData(
        ufuncs,  // func
        data,  // data. The documentation claims we can pass NULL here, but then it segfaults!
        types,  // types
        1,  // ntypes
        12,  // nin
        1,  // nout
        PyUFunc_None,  // identity
        "mt2_ufunc",  // name
        "Numpy ufunc to compute mt2",  // doc
        0  // unused
    );

    PyObject* module_dict = PyModule_GetDict(module);
    PyDict_SetItemString(module_dict, "mt2_ufunc", mt2_ufunc);
    PyDict_SetItemString(module_dict, "__version__", PyUnicode_FromString(MACRO_STRINGIFY(VERSION_INFO)));
    Py_DECREF(mt2_ufunc);

    return module;
}
