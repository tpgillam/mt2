import time

import numpy

from mt2 import mt2, mt2_lally


def main():
    n1 = 400
    n2 = 400

    # Make mass_1 vary over the first axis, and mass_2 vary over the second axis
    mass_1 = numpy.linspace(1, 200, n1).reshape((-1, 1))
    mass_2 = numpy.linspace(1, 200, n2).reshape((1, -1))

    # Pre-allocate output so that we are just timing MT2 itself, not the allocation of
    # the output buffer.
    out = numpy.zeros((n1, n2))

    # `val` has shape (n1, n2), since `mass_1` and `mass_2` broadcast.
    t_lester_start = time.time()
    mt2(
        100, 410, 20,  # Visible 1: mass, px, py
        150, -210, -300,  # Visible 2: mass, px, py
        -200, 280,  # Missing transverse momentum: x, y
        mass_1, mass_2,  # Invisible 1 mass, invisible 2 mass
        out=out)
    t_lester_end = time.time()
    
    t_lally_start = time.time()
    mt2(
        100, 410, 20,  # Visible 1: mass, px, py
        150, -210, -300,  # Visible 2: mass, px, py
        -200, 280,  # Missing transverse momentum: x, y
        mass_1, mass_2,  # Invisible 1 mass, invisible 2 mass
        out=out)
    t_lally_end = time.time()
    
    print("Elapsed time Lester: {} seconds".format(t_lester_end - t_lester_start))
    print("Elapsed time Lally : {} seconds".format(t_lally_end - t_lally_start))


if __name__ == "__main__":
    main()
