/*
 * Asymmetric MT2 with the Lester-Nachman bisection algorithm.
 *
 * Please cite arxiv.org/abs/1411.4312 and arxiv.org/abs/hep-ph/9906349 .
 *
 * C++-subset version.
 */

/*
 * Includes
 *
 * cmath
 *     std::sqrt, std::fabs
 * limits
 *     std::numeric_limits
 */
#include <cmath>
#include <limits>


/* Macros */
#ifdef __GNUC__
#define mt2_rare(expr) __builtin_expect((expr), 0)
#else
#define mt2_rare(expr) (expr)
#endif


/* Types */
/*
 * Parametrize a conic section
 *
 *                [[cxx cxy cx ]    [x
 *      [x y 1] ·  [cxy cyy cy ]  ·  y  == 0
 *                 [cx  cy  c  ]]    1]
 *
 * where cx, cy, and c depend on two variables.
 */
template <typename T>
struct mt2_conic {
    T cxx;
    T cyy;
    T cxy;
    T cx[2];
    T cy[2];
    T c[2];
};

/* Three numbers; used to specify a quadratic. */
template <typename T>
struct mt2_trio {
    T c0;
    T c1;
    T c2;
};


/* Template declarations */
template <typename T>
static struct mt2_conic<T> mt2_ellipse(T m, T px, T py, T ssm, T sspx, T sspy);

template <typename T>
static struct mt2_conic<T> mt2_ellipse_rest(T m, T px, T py, T ssm);

template <typename T>
static struct mt2_trio<T> mt2_det(const struct mt2_conic<T> *a);

template <typename T>
static struct mt2_trio<T> mt2_lester(const struct mt2_conic<T> *a,
                                     const struct mt2_conic<T> *b);

template <typename T>
static bool mt2_disjoint(const struct mt2_trio<T> qs[4], T m, bool *error);

template <typename T>
static inline T mt2_eval_quadratic(const struct mt2_trio<T> *p, T x);

template <typename T>
static inline void mt2_swap(T *x, T *y);


/* Global variables */
static const float mt2_error = -1.38570487f;


/* Template definitions */
/*
 * Return asymmetric MT2, approximated by a bisection method.
 *
 * If a pair of particles with equal mass M which decayed semi-invisibly to
 * `a' + something invisible and `b' + something else invisible, then MT2 is a
 * greatest lower bound on M.
 *
 * Arguments:
 *     am, apx, apy:
 *         mass and transverse momentum components of one visible child
 *     bm, bpx, bpy:
 *         mass and transverse momentum components of the other visible child
 *     sspx, sspy:
 *         missing transverse momentum components
 *     ssam, ssbm
 *         masses of the invisible particles associated with `a' and `b'
 *
 *     Masses are assumed to be non-negative.
 *
 * Returns:
 *     An estimate of MT2, or a negative number if something goes wrong.
 */
template <typename T>
T
mt2_bisect_impl(T am, T apx, T apy,
                T bm, T bpx, T bpy,
                T sspx, T sspy,
                T ssam, T ssbm,
                T precision=0)
{
    /* This physical scale is used for initial bounding and input testing. */
    const auto scale = std::sqrt(0.125f*(
        sspx*sspx + sspy*sspy + (ssam*ssam + ssbm*ssbm)
        + ((apx*apx + apy*apy + am*am) + (bpx*bpx + bpy*bpy + bm*bm))
    ));

    const auto squeeze = 1 / scale;

    /* If scale is 0 or NAN, then mt2 is also. */
    if (mt2_rare(!(scale > 0)))
        return scale;

    /* Sort legs by lower bounds on the parent mass. */
    if (am + ssam > bm + ssbm) {
        mt2_swap(&am, &bm);
        mt2_swap(&apx, &bpx);
        mt2_swap(&apy, &bpy);
        mt2_swap(&ssam, &ssbm);
    }

    /* Squeeze towards 1 to reduce over/underflow risk. */
    am *= squeeze;
    apx *= squeeze;
    apy *= squeeze;
    bm *= squeeze;
    bpx *= squeeze;
    bpy *= squeeze;
    sspx *= squeeze;
    sspy *= squeeze;
    ssam *= squeeze;
    ssbm *= squeeze;

    /* At `lo', the ellipses will be disjoint. */
    auto lo = bm + ssbm;
    auto hi = lo + 1;

    /* Negative masses can cause negative bounds; avoid that case. */
    if (mt2_rare(!(lo > 0)))
        return mt2_error;

    /* Construct the ellipses and their properties as quadratics. */
    auto a_ellipse = mt2_ellipse_rest(am, -apx, -apy, ssam);
    auto b_ellipse = mt2_ellipse(bm, bpx, bpy, ssbm, sspx, sspy);

    const struct mt2_trio<T> quadratics[4] = {
        mt2_det(&a_ellipse),
        mt2_det(&b_ellipse),
        mt2_lester(&a_ellipse, &b_ellipse),
        mt2_lester(&b_ellipse, &a_ellipse),
    };

    /* Expand to find an upper bound. */
    for (;;) {
        bool error;
        const auto disjoint = mt2_disjoint(quadratics, hi, &error);

        if (mt2_rare(error | (hi >= std::numeric_limits<T>::max())))
            return mt2_error;

        if (!disjoint)
            break;

        lo = hi;
        hi *= 2;
    }

    /* Set termination tolerances. If precision is NAN, rel_tol is epsilon. */
    const auto epsilon = std::numeric_limits<T>::epsilon();
    const auto rel_tol = epsilon < precision ? precision : epsilon;
    const auto abs_tol = epsilon;

    /* Bisect; this loop is our fiery pit of hell. */
    for (;;) {
        const auto m = 0.5f*(lo + hi);

        if (mt2_rare(hi <= lo*(1 + 2*rel_tol) + 2*abs_tol))
            return m * scale;

        bool error;
        const auto disjoint = mt2_disjoint(quadratics, m, &error);

        if (disjoint)
            lo = m;
        else
            hi = m;

        if (mt2_rare(error))
            return lo * scale;
    }
}

/*
 * Return a parametrized ellipse for given kinematics.
 */
template <typename T>
static struct mt2_conic<T>
mt2_ellipse(T m, T px, T py, T ssm, T sspx, T sspy)
{
    struct mt2_conic<T> out;
    const auto tx = 2 * px;
    const auto ty = 2 * py;
    const auto m2sum = m*m + ssm*ssm;
    const auto m2dif = m*m - ssm*ssm;
    const auto gx = (m*m*4 + ty*ty)*sspx - tx*ty*sspy;
    const auto gy = (m*m*4 + tx*tx)*sspy - tx*ty*sspx;

    out.cxx = m*m*4 + ty*ty;
    out.cyy = m*m*4 + tx*tx;
    out.cxy = -tx*ty;
    out.cx[0] = -m2sum*tx - gx;
    out.cx[1] = tx;
    out.cy[0] = -m2sum*ty - gy;
    out.cy[1] = ty;
    out.c[0] = (
        sspx*(2*m2sum*tx + gx) + sspy*(2*m2sum*ty + gy)
        + (ssm*ssm*(tx*tx + ty*ty) - m2dif*m2dif)
    );
    out.c[1] = 2*(m2sum - (sspx*tx + sspy*ty));
    return out;
}

/*
 * Special case of `mt2_ellipse' with zero missing momenta.
 *
 * Yes this does make a measurable performance improvement; algebra with zeros
 * is not optimised away (see, e.g., the C17 standard, Appendix F.9.2).
 */
template <typename T>
static struct mt2_conic<T>
mt2_ellipse_rest(T m, T px, T py, T ssm)
{
    struct mt2_conic<T> out;
    const auto tx = 2 * px;
    const auto ty = 2 * py;
    const auto m2sum = m*m + ssm*ssm;
    const auto m2dif = m*m - ssm*ssm;

    out.cxx = m*m*4 + ty*ty;
    out.cyy = m*m*4 + tx*tx;
    out.cxy = -tx*ty;
    out.cx[0] = -m2sum*tx;
    out.cx[1] = tx;
    out.cy[0] = -m2sum*ty;
    out.cy[1] = ty;
    out.c[0] = ssm*ssm*(tx*tx + ty*ty) - m2dif*m2dif;
    out.c[1] = 2*m2sum;
    return out;
}

/*
 * Return the quadratic for the determinant of a parametrized conic.
 *
 * The quadratic part of the c parameter is -1.
 */
template <typename T>
static struct mt2_trio<T>
mt2_det(const struct mt2_conic<T> *a)
{
    struct mt2_trio<T> out;
    const auto xx = a->cxx;
    const auto yy = a->cyy;
    const auto xy = a->cxy;
    const auto x = a->cx;
    const auto y = a->cy;
    const auto c = a->c;

    out.c0 = (
        2*xy*x[0]*y[0]
        - (yy*x[0]*x[0] + xx*y[0]*y[0])
        + c[0]*(xx*yy - xy*xy)
    );

    out.c1 = (
        2*xy*(x[1]*y[0] + x[0]*y[1])
        - 2*(yy*x[0]*x[1] + xx*y[0]*y[1])
        + c[1]*(xx*yy - xy*xy)
    );

    out.c2 = (
        2*xy*x[1]*y[1]
        - (yy*x[1]*x[1] + xx*y[1]*y[1])
        - (xx*yy - xy*xy)
    );
    return out;
}

/*
 * Return the quadratic for the `Lester factor' of two parametrized conics.
 *
 * The quadratic part of the c parameter is -1.
 */
template <typename T>
static struct mt2_trio<T>
mt2_lester(const struct mt2_conic<T> *a, const struct mt2_conic<T> *b)
{
    struct mt2_trio<T> out;
    const auto axx = a->cxx;
    const auto ayy = a->cyy;
    const auto axy = a->cxy;
    const auto ax = a->cx;
    const auto ay = a->cy;
    const auto ac = a->c;
    const auto bxx = b->cxx;
    const auto byy = b->cyy;
    const auto bxy = b->cxy;
    const auto bx = b->cx;
    const auto by = b->cy;
    const auto bc = b->c;

    out.c0 = (
        bxx*(ayy*ac[0] - ay[0]*ay[0])
        + byy*(axx*ac[0] - ax[0]*ax[0])
        + bc[0]*(axx*ayy - axy*axy)
    ) + 2*(
        bx[0]*(axy*ay[0] - ayy*ax[0])
        + by[0]*(axy*ax[0] - axx*ay[0])
        + bxy*(ax[0]*ay[0] - axy*ac[0])
    );

    out.c1 = (
        bxx*(ayy*ac[1] - 2*ay[0]*ay[1])
        + byy*(axx*ac[1] - 2*ax[0]*ax[1])
        + bc[1]*(axx*ayy - axy*axy)
    ) + 2*(
        (bx[0]*(axy*ay[1] - ayy*ax[1]) + bx[1]*(axy*ay[0] - ayy*ax[0]))
        + (by[0]*(axy*ax[1] - axx*ay[1]) + by[1]*(axy*ax[0] - axx*ay[0]))
        + bxy*(ax[0]*ay[1] + ax[1]*ay[0] - axy*ac[1])
    );

    out.c2 = (
        -bxx*(ayy + ay[1]*ay[1])
        - byy*(axx + ax[1]*ax[1])
        - (axx*ayy - axy*axy)
    ) + 2*(
        bx[1]*(axy*ay[1] - ayy*ax[1])
        + by[1]*(axy*ax[1] - axx*ay[1])
        + bxy*(ax[1]*ay[1] + axy)
    );
    return out;
}

/*
 * Are our ellipses disjoint?
 *
 * Ellipse properties are specified as quadratics in mass `m' squared.
 */
template <typename T>
static bool
mt2_disjoint(const struct mt2_trio<T> quadratics[4], T m, bool *error)
{
    auto a_det = mt2_eval_quadratic(quadratics + 0, m*m);
    auto b_det = mt2_eval_quadratic(quadratics + 1, m*m);
    auto a_lester = mt2_eval_quadratic(quadratics + 2, m*m);
    auto b_lester = mt2_eval_quadratic(quadratics + 3, m*m);

    /* Sort sides. */
    if (std::fabs(a_det) < std::fabs(b_det)) {
        mt2_swap(&a_det, &b_det);
        mt2_swap(&a_lester, &b_lester);
    }

    /* Scale to 'monomial form'. */
    const auto a = a_lester / a_det;
    const auto b = b_lester / a_det;
    const auto c = b_det / a_det;

    *error = a_det == 0;

    /* Using some branching logic to aid early escapes. */
    return (
        (a*a > b*3) && 
        ((a < 0) || (b*b*4 > a*a*b + a*c*3)) &&
        (a*c*(b*18 - a*a*4) > c*c*27 + b*b*(b*4 - a*a))
    );
}

/*
 * Evaluate a quadratic with trio coefficients.
 *
 * Using Horner's style for fused-multiply-add opportunities.
 */
template <typename T>
static inline T
mt2_eval_quadratic(const struct mt2_trio<T> *p, T x)
{
    return p->c0 + x*(p->c1 + x*p->c2);
}

/* Using direct style for better results without fused-multiply-add. */
template <>
inline long double
mt2_eval_quadratic(const struct mt2_trio<long double> *p, long double x)
{
    return p->c0 + x*p->c1 + x*x*p->c2;
}

/* Swap values with pointer syntax. */
template <typename T>
void
mt2_swap(T *x, T *y)
{
    const auto tmp = *x;
    *x = *y;
    *y = tmp;
}

/* Clean-up */
#undef mt2_rare
