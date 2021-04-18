from typing import Optional, Union

import numpy

# noinspection PyUnresolvedReferences
from _mt2 import mt2_tombs_ufunc

__author__ = "Thomas Gillam"
__email__ = "tpgillam@googlemail.com"
__version__ = "1.1.0"

__all__ = ["mt2", "mt2_ufunc"]


def mt2(
    m_vis_1: Union[float, numpy.ndarray],
    px_vis_1: Union[float, numpy.ndarray],
    py_vis_1: Union[float, numpy.ndarray],
    m_vis_2: Union[float, numpy.ndarray],
    px_vis_2: Union[float, numpy.ndarray],
    py_vis_2: Union[float, numpy.ndarray],
    px_miss: Union[float, numpy.ndarray],
    py_miss: Union[float, numpy.ndarray],
    m_invis_1: Union[float, numpy.ndarray],
    m_invis_2: Union[float, numpy.ndarray],
    desired_precision_on_mt2: Union[float, numpy.ndarray] = 0.0,
    use_deci_sections_initially: Union[bool, numpy.ndarray] = True,
    *,
    out: Optional[numpy.ndarray] = None
) -> Union[float, numpy.ndarray]:
    """
    Returns asymmetric mT2 (which is >=0), or a negative value if no solution exists.

    We broadcast over any arguments that are provided, following standard numpy
    conventions.

    If more flexibility is required, note that the underlying mt2_ufunc can be used
    directly, which specifies additional arguments (like `where`) in keeping with other
    numpy ufuncs.

    Args:
        m_vis_1: Mass of visible particle 1
        px_vis_1: x-momentum of visible particle 1
        py_vis_1: y-momentum of visible particle 1
        m_vis_2: Mass of visible particle 2
        px_vis_2: x-momentum of visible particle 2
        py_vis_2: y-momentum of visible particle 2
        px_miss: x component of missing momentum
        py_miss: y component of missing momentum
        m_invis_1: Assumed mass of invisible particle 1
        m_invis_2: Assumed mass of invisible particle 2
        desired_precision_on_mt2: This must be non-negative.  If set to zero (default)
            MT2 will be calculated to the highest precision available on the machine (or
            as close to that as the algorithm permits). If set to a positive value,
            MT2 (note that is MT2, not its square) will be calculated to
            within ±desiredPrecisionOnMT2.
            Note that by requesting precision of ±0.01 GeV on an MT2 value of 100 GeV
            can result in speedups of a factor of two to three.
        out: If specified, an array into which the output will be placed.
            Must have dtype numpy.float64.

    Returns:
        MT2 calculated for all inputs. If an array, will have shape that is the result
        of broadcasting all inputs.

        A negative value is returned for any element for which MT2 cannot be computed;
        for example, this will occur if the arguments specify an infeasible optimisation
        problem.
    """
    return mt2_tombs_ufunc(
        m_vis_1,
        px_vis_1,
        py_vis_1,
        m_vis_2,
        px_vis_2,
        py_vis_2,
        px_miss,
        py_miss,
        m_invis_1,
        m_invis_2,
        desired_precision_on_mt2,
        out,
    )


mt2_ufunc = mt2_tombs_ufunc
