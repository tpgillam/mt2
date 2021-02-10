=====
Usage
=====

To use mt2 in a project::

    import mt2
    mt2.get_mT2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB)

A short example to test correct evaluation::

    import mt2

    pxA =   410
    pyA =    20
    mVisA = 100
    chiA = 100

    pxB =  -210
    pyB =  -300
    mVisB = 150
    chiB = 100

    pxMiss = -200
    pyMiss =  280

    computed_mT2_value = mt2.get_mT2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB)
    expected_mT2_value = 412.628

    print("Expected mT2 = "+str(expected_mT2_value))
    print("Computed mT2 = "+str(computed_mT2_value))
