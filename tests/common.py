from _mt2 import mt2_lally_ufunc, mt2_lester_ufunc, mt2_tombs_ufunc


def mt2_lally(*args, desired_precision_on_mt2=0.0, out=None):
    """
    Compute mT2 using the method and code by Colin Lally.

    https://arxiv.org/abs/1509.01831
    """
    return mt2_lally_ufunc(*args, desired_precision_on_mt2, out)


def mt2_lester(
    *args, desired_precision_on_mt2=0.0, use_deci_sections_initially=True, out=None
):
    return mt2_lester_ufunc(
        *args, desired_precision_on_mt2, use_deci_sections_initially, out
    )


def mt2_tombs(*args, desired_precision_on_mt2=0.0, out=None):
    return mt2_tombs_ufunc(*args, desired_precision_on_mt2, out)
