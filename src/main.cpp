#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>

#include "lester_mt2_bisect_v7.h"
#include "mt2_Lallyver2.h"
#include "mt2_bisect.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

static void mt2_lester_ufunc(
    char **args,
// const-correctness was introduced in numpy 1.19, but retain backward compatibility.
#ifdef NPY_1_19_API_VERSION
    npy_intp const *dimensions,
    npy_intp const *steps,
#else
    npy_intp *dimensions,
    npy_intp *steps,
#endif
    void *data)
{
    const npy_intp n = dimensions[0];

    char *mVis1 = args[0];
    char *pxVis1 = args[1];
    char *pyVis1 = args[2];
    char *mVis2 = args[3];
    char *pxVis2 = args[4];
    char *pyVis2 = args[5];
    char *pxMiss = args[6];
    char *pyMiss = args[7];
    char *mInvis1 = args[8];
    char *mInvis2 = args[9];
    char *desiredPrecisionOnMT2 = args[10];
    char *useDeciSectionsInitially = args[11];
    char *out = args[12];

    const npy_intp mVis1_step = steps[0];
    const npy_intp pxVis1_step = steps[1];
    const npy_intp pyVis1_step = steps[2];
    const npy_intp mVis2_step = steps[3];
    const npy_intp pxVis2_step = steps[4];
    const npy_intp pyVis2_step = steps[5];
    const npy_intp pxMiss_step = steps[6];
    const npy_intp pyMiss_step = steps[7];
    const npy_intp mInvis1_step = steps[8];
    const npy_intp mInvis2_step = steps[9];
    const npy_intp desiredPrecisionOnMT2_step = steps[10];
    const npy_intp useDeciSectionsInitially_step = steps[10];
    const npy_intp out_step = steps[12];

    for (npy_intp i = 0; i < n; ++i)
    {
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
            *(npy_bool *)useDeciSectionsInitially);

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

static void mt2_lally_ufunc(
    char **args,
// const-correctness was introduced in numpy 1.19, but retain backward compatibility.
#ifdef NPY_1_19_API_VERSION
    npy_intp const *dimensions,
    npy_intp const *steps,
#else
    npy_intp *dimensions,
    npy_intp *steps,
#endif
    void *data)
{
    const npy_intp n = dimensions[0];

    char *mVis1 = args[0];
    char *pxVis1 = args[1];
    char *pyVis1 = args[2];
    char *mVis2 = args[3];
    char *pxVis2 = args[4];
    char *pyVis2 = args[5];
    char *pxMiss = args[6];
    char *pyMiss = args[7];
    char *mInvis1 = args[8];
    char *mInvis2 = args[9];
    char *desiredPrecisionOnMT2 = args[10];
    char *out = args[11];

    const npy_intp mVis1_step = steps[0];
    const npy_intp pxVis1_step = steps[1];
    const npy_intp pyVis1_step = steps[2];
    const npy_intp mVis2_step = steps[3];
    const npy_intp pxVis2_step = steps[4];
    const npy_intp pyVis2_step = steps[5];
    const npy_intp pxMiss_step = steps[6];
    const npy_intp pyMiss_step = steps[7];
    const npy_intp mInvis1_step = steps[8];
    const npy_intp mInvis2_step = steps[9];
    const npy_intp desiredPrecisionOnMT2_step = steps[10];
    const npy_intp out_step = steps[11];

    for (npy_intp i = 0; i < n; ++i)
    {
        *((double *)out) = mt2_lally(
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
            *(double *)desiredPrecisionOnMT2);

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
        out += out_step;
    }
}

static void mt2_tombs_ufunc(
    char **args,
// const-correctness was introduced in numpy 1.19, but retain backward compatibility.
#ifdef NPY_1_19_API_VERSION
    npy_intp const *dimensions,
    npy_intp const *steps,
#else
    npy_intp *dimensions,
    npy_intp *steps,
#endif
    void *data)
{
    const npy_intp n = dimensions[0];

    char *mVis1 = args[0];
    char *pxVis1 = args[1];
    char *pyVis1 = args[2];
    char *mVis2 = args[3];
    char *pxVis2 = args[4];
    char *pyVis2 = args[5];
    char *pxMiss = args[6];
    char *pyMiss = args[7];
    char *mInvis1 = args[8];
    char *mInvis2 = args[9];
    char *desiredPrecisionOnMT2 = args[10];
    char *out = args[11];

    const npy_intp mVis1_step = steps[0];
    const npy_intp pxVis1_step = steps[1];
    const npy_intp pyVis1_step = steps[2];
    const npy_intp mVis2_step = steps[3];
    const npy_intp pxVis2_step = steps[4];
    const npy_intp pyVis2_step = steps[5];
    const npy_intp pxMiss_step = steps[6];
    const npy_intp pyMiss_step = steps[7];
    const npy_intp mInvis1_step = steps[8];
    const npy_intp mInvis2_step = steps[9];
    const npy_intp desiredPrecisionOnMT2_step = steps[10];
    const npy_intp out_step = steps[11];

    for (npy_intp i = 0; i < n; ++i)
    {
        *((double *)out) = mt2_bisect_impl(
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
            *(double *)desiredPrecisionOnMT2);

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
        out += out_step;
    }
}

/* This a pointer to mt2_lester_ufunc */
PyUFuncGenericFunction mt2_lester_ufuncs[1] = {&mt2_lester_ufunc};

/* These are the input and return dtypes of mt2_lester_ufunc.*/
static char mt2_lester_types[13] = {
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
    NPY_BOOL,   // bool useDeciSectionsInitially=true
    NPY_DOUBLE  // <result>
};

/* This a pointer to mt2_lally_ufunc */
PyUFuncGenericFunction mt2_lally_ufuncs[1] = {&mt2_lally_ufunc};

/* These are the input and return dtypes of mt2_lally_ufunc.*/
static char mt2_lally_types[12] = {
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
    NPY_DOUBLE  // <result>
};

/* This a pointer to mt2_tombs_ufunc */
PyUFuncGenericFunction mt2_tombs_ufuncs[1] = {&mt2_tombs_ufunc};

/* These are the input and return dtypes of mt2_tombs_ufunc.*/
static char mt2_tombs_types[12] = {
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
    NPY_DOUBLE  // <result>
};

PyDoc_STRVAR(mt2_module_doc, "Provides the mt2 stransverse mass ufunc.");

static PyMethodDef methods[] = {
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_mt2",
    mt2_module_doc,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL};

static void *data[1] = {NULL};

PyMODINIT_FUNC PyInit__mt2(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    if (!module)
    {
        return NULL;
    }

    import_array();
    import_ufunc();
    import_umath();

    PyObject *mt2_lester_ufunc = PyUFunc_FromFuncAndData(
        mt2_lester_ufuncs,                 // func
        data,                              // data. The documentation claims we can pass NULL here, but then it segfaults!
        mt2_lester_types,                  // types
        1,                                 // ntypes
        12,                                // nin
        1,                                 // nout
        PyUFunc_None,                      // identity
        "mt2_lester_ufunc",                // name
        "Numpy ufunc to compute mt2 (LN)", // doc
        0                                  // unused
    );

    PyObject *mt2_lally_ufunc = PyUFunc_FromFuncAndData(
        mt2_lally_ufuncs,                              // func
        data,                                          // data. The documentation claims we can pass NULL here, but then it segfaults!
        mt2_lally_types,                               // types
        1,                                             // ntypes
        11,                                            // nin
        1,                                             // nout
        PyUFunc_None,                                  // identity
        "mt2_lally_ufunc",                             // name
        "Numpy ufunc to compute mt2 (by Colin Lally)", // doc
        0                                              // unused
    );

    PyObject *mt2_tombs_ufunc = PyUFunc_FromFuncAndData(
        mt2_tombs_ufuncs,                                                    // func
        data,                                                                // data. The documentation claims we can pass NULL here, but then it segfaults!
        mt2_tombs_types,                                                     // types
        1,                                                                   // ntypes
        11,                                                                  // nin
        1,                                                                   // nout
        PyUFunc_None,                                                        // identity
        "mt2_tombs_ufunc",                                                   // name
        "Numpy ufunc to compute mt2 (LN algo, implemented by Rupert Tombs)", // doc
        0                                                                    // unused
    );

    PyObject *module_dict = PyModule_GetDict(module);
    PyDict_SetItemString(module_dict, "mt2_lester_ufunc", mt2_lester_ufunc);
    PyDict_SetItemString(module_dict, "mt2_lally_ufunc", mt2_lally_ufunc);
    PyDict_SetItemString(module_dict, "mt2_tombs_ufunc", mt2_tombs_ufunc);
    PyDict_SetItemString(module_dict, "__version__", PyUnicode_FromString(MACRO_STRINGIFY(VERSION_INFO)));
    Py_DECREF(mt2_lester_ufunc);
    Py_DECREF(mt2_lally_ufunc);
    Py_DECREF(mt2_tombs_ufunc);

    return module;
}
