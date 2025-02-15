from matplotlib import pyplot

from mt2 import mt2
from mt2._mt2 import mt2_lally_ufunc
from mt2.diagnostics import make_ellipses


def plot_lester_v_lally(args: tuple[float, ...]) -> None:
    """Make a plot comparing Lester & Lally algorithms to MT2.

    `args` should be a tuple of arguments that would be given to `mt2` or `mt2_lally`.
    """
    m_lester = mt2(*args)
    desired_precision_on_mt2 = 0.0
    m_lally = mt2_lally_ufunc(*args, desired_precision_on_mt2)

    params_1_lester, params_2_lester = make_ellipses(m_lester, *args)
    params_1_lally, params_2_lally = make_ellipses(m_lally, *args)

    perturbation = 0.001
    params_1_lesterpp, params_2_lesterpp = make_ellipses(
        m_lester * (1 + perturbation), *args
    )
    params_1_lestermm, params_2_lestermm = make_ellipses(
        m_lester * (1 - perturbation), *args
    )

    params_1_lallypp, params_2_lallypp = make_ellipses(
        m_lally * (1 + perturbation), *args
    )
    params_1_lallymm, params_2_lallymm = make_ellipses(
        m_lally * (1 - perturbation), *args
    )

    pyplot.figure(figsize=(8, 6))
    pyplot.title(f"MT2 Lester={m_lester:.3f};   Lally={m_lally:.3f}")

    pyplot.plot(*params_1_lester.get_points(100).T, c="C0", label="Lester 1")
    pyplot.plot(
        *params_1_lesterpp.get_points(100).T,
        c="C0",
        alpha=0.5,
        label=f"Lester 1 ({perturbation})",
    )
    pyplot.plot(
        *params_1_lestermm.get_points(100).T,
        c="C0",
        alpha=0.5,
    )
    pyplot.plot(*params_2_lester.get_points(100).T, c="C1", label="Lester 2")
    pyplot.plot(
        *params_2_lesterpp.get_points(100).T,
        c="C1",
        alpha=0.5,
        label=f"Lester 2 ({perturbation})",
    )
    pyplot.plot(
        *params_2_lestermm.get_points(100).T,
        c="C1",
        alpha=0.5,
    )

    pyplot.plot(*params_1_lally.get_points(100).T, c="C0", ls="dashed", label="Lally 1")
    pyplot.plot(*params_2_lally.get_points(100).T, c="C1", ls="dashed", label="Lally 2")
    pyplot.plot(
        *params_1_lallypp.get_points(100).T,
        c="C0",
        alpha=0.5,
        ls="dashed",
        label=f"Lally 1 ({perturbation})",
    )
    pyplot.plot(
        *params_1_lallymm.get_points(100).T,
        c="C0",
        alpha=0.5,
        ls="dashed",
    )
    pyplot.plot(*params_2_lally.get_points(100).T, c="C1", label="Lally 2")
    pyplot.plot(
        *params_2_lallypp.get_points(100).T,
        c="C1",
        alpha=0.5,
        ls="dashed",
        label=f"Lally 2 ({perturbation})",
    )
    pyplot.plot(
        *params_2_lallymm.get_points(100).T,
        c="C1",
        alpha=0.5,
        ls="dashed",
    )

    pyplot.legend()
    pyplot.xlabel("p1x")
    pyplot.ylabel("p1y")
    pyplot.tight_layout()


def main():
    bad_args = (
        29.359184426449335,
        -41.66748425813935,
        2.4727555091581337,
        68.86316293078755,
        88.10079057588052,
        -39.0193751229873,
        27.644883008122292,
        -29.00097683721164,
        69.60314069891137,
        28.917825385644313,
    )

    # m1 = 0.00473856926 GeV, p1x = 110.43799970914 GeV,
    # p1y = −213.46687262192 GeV, m2 = 0.0 GeV, p2x = −20.28455035002 GeV, p2y =
    # 235.76522546534 GeV, p/x = 111.30684472357 GeV, p/yT = 37.47049084405 GeV, μN 1 =
    # μN 2 = 0.0 GeV
    args_event_1 = (
        0.00473856926,
        110.43799970914,
        -213.46687262192,
        0.0,
        -20.28455035002,
        235.76522546534,
        111.30684472357,
        37.47049084405,
        0.0,
        0.0,
    )

    nice_args = (
        100,
        410,
        20,  # Visible 1: mass, px, py
        150,
        -210,
        -300,  # Visible 2: mass, px, py
        -200,
        280,  # Missing transverse momentum: x, y
        100,
        100,
    )  # Invisible 1 mass, invisible 2 mass

    plot_lester_v_lally(nice_args)
    plot_lester_v_lally(bad_args)
    # TODO This doesn't work yet, as we don't support degenerate ellipses.
    if False:
        plot_lester_v_lally(args_event_1)
    pyplot.show()


if __name__ == "__main__":
    main()
