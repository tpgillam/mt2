def test_3():
    # Test lots of collinear endpoint cases.
    import random, math, pytest
    from mt2 import mt2
    random.seed(0) # If the test fails, we want to be able to repeat the test!
    for i in range(10000):
        m_visA = random.uniform(0,10)
        m_visB = random.uniform(0,10)
        m_invisA = random.uniform(0,10)
        m_invisB = random.uniform(0,10)
        m_parent = max(m_visA+m_invisA, m_visB+m_invisB) + random.uniform(0.1,10) # 0.1 prevents possibility of division by zero at MOOO
        p_parentA = random.uniform(0,10)
        p_parentB = random.uniform(0,10)
        E_parentA = math.sqrt(p_parentA**2 + m_parent**2)
        E_parentB = math.sqrt(p_parentB**2 + m_parent**2)
        beta_A = p_parentA/E_parentA  # MOO
        beta_B = p_parentB/E_parentB  # MOO
        gamma_A = 1.0/math.sqrt(1-beta_A**2)
        gamma_B = 1.0/math.sqrt(1-beta_B**2)
        pA = math.sqrt((m_parent-m_visA-m_invisA)*
                       (m_parent+m_visA-m_invisA)*
                       (m_parent-m_visA+m_invisA)*
                       (m_parent+m_visA+m_invisA))/(2*m_parent)
        pB = math.sqrt((m_parent-m_visB-m_invisB)*
                       (m_parent+m_visB-m_invisB)*
                       (m_parent-m_visB+m_invisB)*
                       (m_parent+m_visB+m_invisB))/(2*m_parent)
        p_visA_boosted = gamma_A * (beta_A*math.sqrt(m_visA**2 + pA**2) + pA)
        p_visB_boosted = gamma_B * (beta_B*math.sqrt(m_visB**2 + pB**2) + pB)
        p_invisA_boosted = gamma_A * (beta_A*math.sqrt(m_invisA**2 + pA**2) - pA)
        p_invisB_boosted = gamma_B * (beta_B*math.sqrt(m_invisB**2 + pB**2) - pB)
        p_miss = p_invisA_boosted + p_invisB_boosted
        theta = random.uniform(0, math.tau)
        c = math.cos(theta)
        s = math.sin(theta)
        pxmiss, pymiss = p_miss*c,  p_miss*s
        ax, ay = p_visA_boosted*c, p_visA_boosted*s
        bx, by = p_visB_boosted*c, p_visB_boosted*s
        val = mt2(m_visA, ax, ay,
                  m_visB, bx, by,
                  pxmiss, pymiss,
                  m_invisA, m_invisB)
        is_ok = (val == pytest.approx(m_parent, rel=1e-13))    # passes with rel=1e-12 but fails with rel=1e-13
        if (not is_ok):
            print(f"WARNING! Expected {m_parent} from collinear event but instead got {m_parent} for mt2({m_visA},{ax},{ay}, {m_visB},{bx},{by}, {pxmiss},{pymiss}, {m_invisA},{m_invisB}) in test case {i}.")
        assert is_ok
