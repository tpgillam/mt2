import pathlib
import time

import numpy

from matplotlib import pyplot

from mt2._mt2 import mt2_lally_ufunc, mt2_lester_ufunc, mt2_tombs_ufunc


def mt2_lester(
    *args, desired_precision_on_mt2=0.0, use_deci_sections_initially=True, out=None
):
    return mt2_lester_ufunc(
        *args, desired_precision_on_mt2, use_deci_sections_initially, out
    )


def mt2_lally(*args, desired_precision_on_mt2=0.0, out=None):
    return mt2_lally_ufunc(*args, desired_precision_on_mt2, out)


def mt2_tombs(*args, desired_precision_on_mt2=0.0, out=None):
    return mt2_tombs_ufunc(*args, desired_precision_on_mt2, out)


def _run_profile(args, shape):
    # Pre-allocate output so that we are just timing MT2 itself, not the allocation of
    # the output buffers.
    out_lester = numpy.zeros(shape)
    out_lester_no_ds = numpy.zeros(shape)
    # out_lally = numpy.zeros((n1, n2))
    out_tombs = numpy.zeros(shape)

    # `val` has shape (n1, n2), since `mass_1` and `mass_2` broadcast.
    t_lester_start = time.time()
    mt2_lester(*args, out=out_lester)
    t_lester_end = time.time()

    t_lester_no_ds_start = time.time()
    mt2_lester(*args, out=out_lester_no_ds, use_deci_sections_initially=False)
    t_lester_no_ds_end = time.time()

    # t_lally_start = time.time()
    # mt2_lally(*args, out=out_lally)
    # t_lally_end = time.time()

    t_tombs_start = time.time()
    mt2_tombs(*args, out=out_tombs)
    t_tombs_end = time.time()

    # Check that we get the same thing.
    numpy.testing.assert_array_almost_equal(out_lester, out_lester_no_ds)
    # numpy.testing.assert_array_almost_equal(out_lester, out_lally)
    numpy.testing.assert_array_almost_equal(out_lester, out_tombs)

    t_lester = t_lester_end - t_lester_start
    t_lester_no_ds = t_lester_no_ds_end - t_lester_no_ds_start
    # t_lally = t_lally_end - t_lally_start
    t_tombs = t_tombs_end - t_tombs_start

    return t_lester, t_lester_no_ds, t_tombs


def _print_profile_results(results):
    t_lester, t_lester_no_ds, t_tombs = results
    print("Elapsed time Lester        : {} seconds".format(t_lester))
    print("Elapsed time Lester (no DS): {} seconds".format(t_lester_no_ds))
    # print("Elapsed time Lally         : {} seconds".format())
    print("Elapsed time Tombs         : {} seconds".format(t_tombs))


def _run_and_plot(args, shape):
    n = 1
    for dim in shape:
        n *= dim

    t_lester = []
    t_lester_no_ds = []
    t_tombs = []
    for _ in range(10000):
        t1, t2, t3 = _run_profile(args, shape)
        t_lester.append(t1 / n)
        t_lester_no_ds.append(t2 / n)
        t_tombs.append(t3 / n)

    pyplot.hist(t_lester, bins=100, label="lester", alpha=0.5, log=True)
    pyplot.hist(t_lester_no_ds, bins=100, label="lester (no DS)", alpha=0.5, log=True)
    pyplot.hist(t_tombs, bins=100, label="tombs", alpha=0.5, log=True)
    pyplot.xlabel("Time / evaluation / s")
    pyplot.legend()


def make_plot(output_dir: pathlib.Path) -> None:
    shape = (100,)

    args = (numpy.full(shape, 100), 410, 20, 150, -210, -300, -200, 280, 200, 200)
    pyplot.figure()
    _run_and_plot(args, shape)
    pyplot.title("Example 1")
    pyplot.tight_layout()
    pyplot.savefig(output_dir / "example_1.png")

    args = (numpy.full(shape, 10), 20, 30, 10, -20, -30, -5, -5, 4, 7)
    pyplot.figure()
    _run_and_plot(args, shape)
    pyplot.title("Example 2")
    pyplot.tight_layout()
    pyplot.savefig(output_dir / "example_2.png")


def print_timing_example():
    n1 = 400
    n2 = 400

    # Make mass_1 vary over the first axis, and mass_2 vary over the second axis.
    print("EXAMPLE FROM README:")
    mass_1 = numpy.linspace(1, 200, n1).reshape((-1, 1))
    mass_2 = numpy.linspace(1, 200, n2).reshape((1, -1))
    args = (
        100,
        410,
        20,  # Visible 1: mass, px, py
        150,
        -210,
        -300,  # Visible 2: mass, px, py
        -200,
        280,  # Missing transverse momentum: x, y
        mass_1,
        mass_2,  # Invisible 1 mass, invisible 2 mass
    )
    results = _run_profile(args, (n1, n2))
    _print_profile_results(results)
    print()
    print()

    # print("EXAMPLE FROM LESTER HEADER FILE:")
    mass_1 = numpy.full((n1, 1), 4)
    mass_2 = numpy.full((1, n2), 7)
    args = (10, 20, 30, 10, -20, -30, -5, -5, mass_1, mass_2)
    results = _run_profile(args, (n1, n2))
    _print_profile_results(results)


def main():
    # print_timing_example()
    output_dir = pathlib.Path.cwd() / "output" / "benchmark_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)
    make_plot(output_dir)


if __name__ == "__main__":
    main()
