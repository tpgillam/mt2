/****************************************************************************************/
/*                                                                                      */
/*                  Finding Asymmetric MT2 by Characteristic Eq Cubic Discriminant    	*/
/*        		                                                          				*/
/*                  Copyright: Colin Lally                                              */
/*                  December 28, 2015, ver2                                          	*/
/*                                                                                    	*/
/****************************************************************************************

 *  version 2: arXiv:1509.01831v2
 *    * minor modifications to lines 395 to 415 to allow for a more accurate pass-through
 *      of current root upper bound value from Newton-Raphson to Regula Falsi methods if
 *      NR gets stuck in a non-converging cycle.
 *
 *  version 1: arXiv:1509.01831.v1
 *    * initial public release


*****************************************************************************************
  Notes:

    1. This is a simple, standalone header file containing the function mt2_lally.

    2. Code assumes that when run it is called by a program producing or containing the mass and momentum values for any number of events in the following order:

       a mass, pax momentum, pay momentum, b mass, pbx momentum, pby momentum, pmissx momentum, pmissy momentum, side a invisible mass, side b invisible mass, desired precision (default value of 0)

    3. Here is an example of it's use (note input variables have been deliberately kept the same as for the algorithm due to C.G.Lester contained in header file lester_mt2_bisect.h for ease of comparison - see http://arxiv.org/abs/1411.4312):

    double mVisA = 10; // mass of visible object on side A.  Should be >=0.
    double pxA = 20; // x momentum of visible object on side A.
    double pyA = 30; // y momentum of visible object on side A.

    double mVisB = 10; // mass of visible object on side B.  Should be >=0.
    double pxB = -20; // x momentum of visible object on side B.
    double pyB = -30; // y momentum of visible object on side B.

    double pxMiss = -5; // x component of missing transverse momentum.
    double pyMiss = -5; // y component of missing transverse momentum.

    double chiA = 4; // hypothesised mass of invisible on side A.  Must be >=0.
    double chiB = 7; // hypothesised mass of invisible on side B.  Must be >=0.

    double desiredPrecisionOnMt2 = 0; // Must be >=0.  If 0 alg aims for precision of 1E-15.  If >0, MT2 computed to supplied absolute precision (but note there is unlikely to be any speed benefit in doing so!)

    double MT2 =  mt2_lally(
           mVisA, pxA, pyA,
           mVisB, pxB, pyB,
           pxMiss, pyMiss,
           chiA, chiB,
           desiredPrecisionOnMt2);

*****************************************************************************************/

#include <iostream>
#include <cmath>

// First we set up a couple of structures to hold the cubic lambda function coefficients (which are each quadratic functions of delta) and the 8th-order discriminant polynomial coefficients
struct DiscriminantCoeffs
{
    double Coeffs8;
    double Coeffs7;
    double Coeffs6;
    double Coeffs5;
    double Coeffs4;
    double Coeffs3;
    double Coeffs2;
    double Coeffs1;
    double Coeffs0;
};

struct CubicCoeffs
{
    double Coeffa2;
    double Coeffa1;
    double Coeffa0;
    double Coeffb2;
    double Coeffb1;
    double Coeffb0;
    double Coeffc2;
    double Coeffc1;
    double Coeffc0;
    double Coeffd2;
    double Coeffd1;
    double Coeffd0;
};

// Define the key functions used by the algorithm - first one is self-evident, second one contains the Regula Falsi root-finding method, the third function is the check to see we have the right form of the lambda function
double NewtonRootFinder(double, double, const DiscriminantCoeffs &, const CubicCoeffs &, double);
double NewDeltaFinder(double, double, int, int, const DiscriminantCoeffs &, const CubicCoeffs &, double);
int lambdaSgnchanges(double, const CubicCoeffs &);

double mt2_lally(

    double ma, // side "a" visible variables
    double pax,
    double pay,

    double mb, // side "b" visible variables
    double pbx,
    double pby,

    double pmissx, // missing transverse momenta
    double pmissy,

    double mna, // side a invisible mass input/guess
    double mnb, // side b invisible mass input/guess

    const double desiredPrecisionOnMt2)
{
    ma = fabs(ma); // just in case input masses are negative - assume negative signs are erroneous
    mb = fabs(mb);
    double MT2 = 0.0;
    bool massless_flag = false;
    const double precision = std::max(desiredPrecisionOnMt2, 0.00000000000001); // if desiredPrecision set to zero, then use this precision of 1E-14

    double masq = ma * ma;
    double Easq = masq + pax * pax + pay * pay;

    double mbsq = mb * mb;
    double Ebsq = mbsq + pbx * pbx + pby * pby;

    // set a-side masses >= b-side masses or if Ma=Mb then set Ea >= Eb
    if (((ma + mna) < (mb + mnb)) || (((ma + mna) == (mb + mnb)) && (Easq < Ebsq)))
    {
        double temp;
        temp = pax;
        pax = pbx;
        pbx = temp;
        temp = pay;
        pay = pby;
        pby = temp;
        temp = Easq;
        Easq = Ebsq;
        Ebsq = temp;
        temp = masq;
        masq = mbsq;
        mbsq = temp;
        temp = ma;
        ma = mb;
        mb = temp;
        temp = mna;
        mna = mnb;
        mnb = temp;
    }

    const double Ea = sqrt(Easq);
    const double Eb = sqrt(Ebsq);

    const double mnasq = mna * mna;
    const double mnbsq = mnb * mnb;

    if ((ma == 0) && (mb == 0) && (mna == 0) && (mnb == 0))
    {
        massless_flag = true;
    }

    // find the coefficients for the two ellipse equations P(p1x,p1y) and Q(p1x, p1y)
    // coefficients for linear and constant terms are polynomials of delta=(Deltasq-masq)/(2Easq)
    // note some of these coefficients are slightly different from those that appear in Cheng & Han's code by factors of 2.
    // This is simply due to the definition of the coefficients of the conic i.e. coefficients B, D and E are defined as per this
    // implementation, but as 2*B, 2*D and 2*E as per Cheng & Han - makes no difference, purely down to convention used.

    const double Ap = 1 - pax * pax / Easq;
    const double Bp = -2. * pax * pay / Easq;
    const double Cp = 1 - pay * pay / Easq;
    const double Dp = -2. * pax;
    const double Ep = -2. * pay;
    const double Fp = -Easq; //quadratic term ie *delta^2
    const double Aq = 1 - pbx * pbx / Ebsq;
    const double Bq = -2. * pbx * pby / Ebsq;
    const double Cq = 1 - pby * pby / Ebsq;
    const double Dqii = 2. * Easq * pbx / Ebsq; //linear term ie *delta^1
    const double Dqi = (2. * (mnasq + masq - mnbsq - mbsq) * pbx) / (2. * Ebsq) - 2. * pmissx +
                       (2. * pbx * (pbx * pmissx + pby * pmissy)) / Ebsq; //constant term ie *delta^0
    const double Eqii = 2. * (Easq * pby) / Ebsq;                         //linear term ie *delta^1
    const double Eqi = (2. * (mnasq + masq - mnbsq - mbsq) * pby) / (2. * Ebsq) - 2. * pmissy +
                       (2. * pby * (pbx * pmissx + pby * pmissy)) / Ebsq;                                                     //constant term ie *delta^0
    const double Fqiii = -Easq * Easq / Ebsq;                                                                                 //quadratic term ie *delta^2
    const double Fqii = (-2. * Easq * ((mnasq + masq - mnbsq - mbsq) / (2. * Eb) + (pbx * pmissx + pby * pmissy) / Eb)) / Eb; //linear term ie *delta^1
    const double Fqi = mnbsq + pmissx * pmissx + pmissy * pmissy -
                       ((mnasq + masq - mnbsq - mbsq) / (2. * Eb) + (pbx * pmissx + pby * pmissy) / Eb) * ((mnasq + masq - mnbsq - mbsq) / (2. * Eb) + (pbx * pmissx + pby * pmissy) / Eb); //constant term ie *delta^0

    // find coefficients for cubic characteristic equation f(lambda) = det(lambda*P - Q)
    // (using the Chris Lester copyrighted observation (patent pending) that multiplying by 4 is faster than by 0.25!)

    const double a2 = (4. * Ap * Cp * Fp + Bp * Dp * Ep - Bp * Bp * Fp - Ap * Ep * Ep - Cp * Dp * Dp); //quadratic term of delta in lambda cubed
    const double a1 = 0.0;                                                                             //linear term of delta in lambda cubed
    const double a0 = (4. * Ap * Cp - Bp * Bp) * mnasq;                                                //constant term of delta in lambda cubed
    const double b2 = -4. * Ap * Cp * Fqiii - 4 * (Ap * Cq + Aq * Cp) * Fp - Bq * Dp * Ep + Bp * (-Dp * Eqii - Ep * Dqii) +
                      2. * Ap * Ep * Eqii + Aq * Ep * Ep + Bp * Bp * Fqiii + 2. * Bp * Bq * Fp + 2. * Cp * Dp * Dqii +
                      Cq * Dp * Dp;                                                                                                  //quadratic term of delta in lambda squared
    const double b1 = -4. * Ap * Cp * Fqii + Bp * (-Dp * Eqi - Ep * Dqi) + 2. * Ap * Ep * Eqi + Bp * Bp * Fqii + 2. * Cp * Dp * Dqi; //linear term of delta in lambda squared
    const double b0 = -4. * Ap * Cp * Fqi - 4. * (Ap * Cq + Aq * Cp) * mnasq + Bp * Bp * Fqi + 2. * Bp * Bq * mnasq;                 //constant term of delta in lambda squared
    const double c2 = 4. * Aq * Cq * Fp + 4. * (Ap * Cq + Aq * Cp) * Fqiii + Bp * Dqii * Eqii + Bq * (Ep * Dqii + Dp * Eqii) -
                      Ap * Eqii * Eqii - 2. * Aq * Ep * Eqii - 2. * Bp * Bq * Fqiii - Bq * Bq * Fp - Cp * Dqii * Dqii - 2. * Cq * Dp * Dqii; //quadratic term of delta in lambda
    const double c1 = 4. * (Ap * Cq + Aq * Cp) * Fqii + Bp * (Dqii * Eqi + Dqi * Eqii) + Bq * (Ep * Dqi + Dp * Eqi) -
                      2. * Ap * Eqi * Eqii - 2. * Aq * Ep * Eqi - 2. * Bp * Bq * Fqii - 2. * Cp * Dqi * Dqii - 2. * Cq * Dp * Dqi;                                     //linear term of delta in lambda
    const double c0 = 4. * Aq * Cq * mnasq + 4. * (Ap * Cq + Aq * Cp) * Fqi + Bp * Dqi * Eqi - Ap * Eqi * Eqi - 2. * Bp * Bq * Fqi - Bq * Bq * mnasq - Cp * Dqi * Dqi; //constant term of delta in lambda
    const double d2 = -4. * Aq * Cq * Fqiii - Bq * Dqii * Eqii + Aq * Eqii * Eqii + Cq * Dqii * Dqii + Bq * Bq * Fqiii;                                                //quadratic term of delta in constant term in f(lambda)
    const double d1 = -4. * Aq * Cq * Fqii - Bq * (Dqi * Eqii + Dqii * Eqi) + 2. * Aq * Eqi * Eqii + 2. * Cq * Dqi * Dqii + Bq * Bq * Fqii;                            //linear term of delta in constant term in f(lambda)
    const double d0 = -4. * Aq * Cq * Fqi - Bq * Dqi * Eqi + Aq * Eqi * Eqi + Cq * Dqi * Dqi + Bq * Bq * Fqi;                                                          //constant term of delta in constant term in f(lambda)

    // Now we find the discriminant of f(lambda)
    // The discriminant will be an eight-order polynomial in delta, and the lowest positive root of this function will be the delta value that gives MT2
    // The coefficients for the eight-order polynomial are as follows (skip if you are of a sensitive nature!):

    double disc8 = 18 * a2 * b2 * c2 * d2 - 4 * b2 * b2 * b2 * d2 + b2 * b2 * c2 * c2 - 4 * a2 * c2 * c2 * c2 - 27 * a2 * a2 * d2 * d2; //eight-order term of delta
    double disc7 = 18 * (a2 * b1 + b2 * a1) * c2 * d2 + 18 * (c2 * d1 + d2 * c1) * a2 * b2 - 4 * b2 * b2 * b2 * d1 - 12 * b2 * b2 * b1 * d2 + 2 * b2 * b2 * c1 * c2 +
                   2 * b2 * b1 * c2 * c2 - 4 * a1 * c2 * c2 * c2 - 12 * a2 * c2 * c2 * c1 - 54 * a2 * a2 * d2 * d1 - 54 * d2 * d2 * a2 * a1; //seventh-order term of delta
    double disc6 = 18 * (a2 * b2 * (c2 * d0 + c1 * d1 + c0 * d2) + c2 * d2 * (a2 * b0 + a1 * b1 + a0 * b2) + (a2 * b1 + b2 * a1) * (c2 * d1 + d2 * c1)) -
                   4 * (b2 * b2 * b2 * d0 + 3 * b2 * b2 * b1 * d1 + d2 * (3 * b2 * b2 * b0 + 3 * b2 * b1 * b1)) + 2 * b2 * b2 * c2 * c0 + 2 * b2 * b0 * c2 * c2 +
                   4 * b2 * b1 * c2 * c1 + b2 * b2 * c1 * c1 + c2 * c2 * b1 * b1 - 4 * (a2 * (3 * c2 * c2 * c0 + 3 * c2 * c1 * c1) + 3 * a1 * c2 * c2 * c1 + a0 * c2 * c2 * c2) - 54 * (a2 * a2 * d2 * d0 + d2 * d2 * a2 * a0 + 2 * a2 * a1 * d2 * d1) - 27 * (a2 * a2 * d1 * d1 + a1 * a1 * d2 * d2); //sixth-order term of delta
    double disc5 = 18 * (a2 * b2 * (c1 * d0 + c0 * d1) + c2 * d2 * (a1 * b0 + a0 * b1) + (a2 * b1 + b2 * a1) * (c2 * d0 + c1 * d1 + c0 * d2) +
                         (c2 * d1 + d2 * c1) * (a2 * b0 + a1 * b1 + a0 * b2)) -
                   4 * (3 * b2 * b2 * b1 * d0 + d1 * (3 * b2 * b2 * b0 + 3 * b2 * b1 * b1) +
                        d2 * (6 * b2 * b1 * b0 + b1 * b1 * b1)) +
                   2 * b2 * b2 * c1 * c0 + 2 * b1 * b0 * c2 * c2 + 4 * b2 * b1 * c2 * c0 + 4 * c2 * c1 * b2 * b0 +
                   2 * b2 * b1 * c1 * c1 + 2 * c2 * c1 * b1 * b1 - 4 * (3 * a0 * c2 * c2 * c1 + a1 * (3 * c2 * c2 * c0 + 3 * c2 * c1 * c1) + a2 * (6 * c2 * c1 * c0 + c1 * c1 * c1)) -
                   54 * (a2 * a2 * d1 * d0 + d2 * d2 * a1 * a0 + 2 * a2 * a1 * d2 * d0 + 2 * a2 * a0 * d2 * d1 + a2 * a1 * d1 * d1 + a1 * a1 * d2 * d1); //fifth-order term of delta
    double disc4 = 18 * (a2 * b2 * c0 * d0 + a0 * b0 * c2 * d2 + (a2 * b0 + a1 * b1 + a0 * b2) * (c2 * d0 + c1 * d1 + c0 * d2) + (a2 * b1 + a1 * b2) * (c1 * d0 + c0 * d1) +
                         (a1 * b0 + a0 * b1) * (c1 * d2 + c2 * d1)) -
                   4 * (d0 * (3 * b2 * b2 * b0 + 3 * b2 * b1 * b1) + d1 * (6 * b2 * b1 * b0 + b1 * b1 * b1) +
                        d2 * (3 * b2 * b0 * b0 + 3 * b1 * b1 * b0)) +
                   b2 * b2 * c0 * c0 + c2 * c2 * b0 * b0 + 4 * b2 * b1 * c1 * c0 + 4 * b1 * b0 * c2 * c1 + 4 * b2 * b0 * c2 * c0 +
                   2 * b2 * b0 * c1 * c1 + 2 * c2 * c0 * b1 * b1 + b1 * b1 * c1 * c1 - 4 * (a0 * (3 * c2 * c2 * c0 + 3 * c2 * c1 * c1) + a1 * (6 * c2 * c1 * c0 + c1 * c1 * c1) + a2 * (3 * c2 * c0 * c0 + 3 * c0 * c1 * c1)) - 27 * (a2 * a2 * d0 * d0 + d2 * d2 * a0 * a0 + 4 * a2 * a1 * d1 * d0 + 4 * a1 * a0 * d2 * d1 + 4 * a2 * a0 * d2 * d0 + a1 * a1 * d1 * d1 + 2 * a2 * a0 * d1 * d1 + 2 * a1 * a1 * d2 * d0); //fourth-order term of delta
    double disc3 = 18 * (c0 * d0 * (a2 * b1 + a1 * b2) + a0 * b0 * (c2 * d1 + c1 * d2) + (a2 * b0 + a1 * b1 + a0 * b2) * (c1 * d0 + c0 * d1) +
                         (a1 * b0 + a0 * b1) * (c2 * d0 + c1 * d1 + c0 * d2)) -
                   4 * (d0 * (6 * b2 * b1 * b0 + b1 * b1 * b1) + d1 * (3 * b2 * b0 * b0 + 3 * b1 * b1 * b0) +
                        3 * b1 * b0 * b0 * d2) +
                   2 * b2 * b1 * c0 * c0 + 2 * b0 * b0 * c2 * c1 + 4 * b2 * b0 * c1 * c0 + 4 * b1 * b0 * c2 * c0 + 2 * c1 * c0 * b1 * b1 + 2 * b1 * b0 * c1 * c1 -
                   4 * (a1 * (3 * c2 * c0 * c0 + 3 * c1 * c1 * c0) + a0 * (6 * c2 * c1 * c0 + c1 * c1 * c1) + 3 * a2 * c1 * c0 * c0) -
                   54 * (a2 * a1 * d0 * d0 + a0 * a0 * d2 * d1 + 2 * a2 * a0 * d1 * d0 + 2 * a1 * a0 * d2 * d0 + a1 * a0 * d1 * d1 + a1 * a1 * d1 * d0); //third-order term of delta
    double disc2 = 18 * (c0 * d0 * (a2 * b0 + a1 * b1 + a0 * b2) + a0 * b0 * (c2 * d0 + c1 * d1 + c0 * d2) + (a1 * b0 + a0 * b1) * (c1 * d0 + c0 * d1)) -
                   4 * (d0 * (3 * b2 * b0 * b0 + 3 * b1 * b1 * b0) + 3 * b1 * b0 * b0 * d1 + b0 * b0 * b0 * d2) + 2 * b2 * b0 * c0 * c0 + 2 * b0 * b0 * c2 * c0 +
                   4 * b1 * b0 * c1 * c0 + c0 * c0 * b1 * b1 + b0 * b0 * c1 * c1 - 4 * (a2 * c0 * c0 * c0 + 3 * a1 * c1 * c0 * c0 + a0 * (3 * c2 * c0 * c0 + 3 * c1 * c1 * c0)) -
                   54 * (a2 * a0 * d0 * d0 + a0 * a0 * d2 * d0 + 2 * a1 * a0 * d1 * d0) - 27 * (a0 * a0 * d1 * d1 + a1 * a1 * d0 * d0); //second-order term of delta
    double disc1 = 18 * (c0 * d0 * (a1 * b0 + a0 * b1) + a0 * b0 * (c1 * d0 + c0 * d1)) - 4 * (3 * b1 * b0 * b0 * d0 + b0 * b0 * b0 * d1) +
                   2 * b1 * b0 * c0 * c0 + 2 * b0 * b0 * c1 * c0 - 4 * (a1 * c0 * c0 * c0 + 3 * a0 * c1 * c0 * c0) - 54 * (a1 * a0 * d0 * d0 + a0 * a0 * d1 * d0); //first-order term of delta
    double disc0 = 18 * a0 * b0 * c0 * d0 - 4 * b0 * b0 * b0 * d0 + b0 * b0 * c0 * c0 - 4 * a0 * c0 * c0 * c0 - 27 * a0 * a0 * d0 * d0;                            //constant term of delta

    // Next we need to check for unbalanced events, we will use different methods depending on elliptic or parabolic (ie massless) situations
    const double tinyValue = std::min(precision, 0.00000000000001);
    bool unbalanced_flag = false;
    double mt2_unbalanced_value = 0.0;
    double delta0 = ma * mna / Easq; //The minimum mass square delta value to have two ellipses (i.e. the "turn-on" mass delta of the larger mass ellipse) - will give a lower bound for MT2
    if (delta0 != delta0)
    { // in case delta0 is undefined/infinite - can happen if input data is suspect i.e. side a momenta and mass (and thus Ea) all = 0
        delta0 = 0.0;
    }
    delta0 = std::max(delta0, tinyValue); //want to avoid lower bound of exactly zero if not unbalanced event

    double delta = 0.0;

    // find the intersection mu_Y mass of the smaller mass ellipse with the "turn-on" mass of the larger ellipse, this will give both a possible upper bound for MT2 and
    // allow a check for the unbalanced configuration in the elliptical case i.e. if intersection mu_Y mass is smaller than delta0, then we have an unbalanced situation.
    // (note intersection of ellipses methodology used is that due to Walker, J.W. see arXiv 1311:6219 for more details)
    double deltaIntersect = 0;
    if (!massless_flag)
    {
        double p1x_a = 1.0E+20 * pax / fabs(pax); //x-coordinate of "turn-on" mass of larger mass ellipse set initially to some large but not infinite quantity
        double p1y_a = 1.0E+20 * pay / fabs(pay); //y-coordinate of "turn-on" mass of larger mass ellipse
        if (ma != 0)
        {
            p1x_a = (mna / ma) * pax; //actual x-coordinate of "turn-on" mass of larger mass ellipse
            p1y_a = (mna / ma) * pay; //actual y-coordinate of "turn-on" mass of larger mass ellipse
        }
        double alpha = 2. * Fqiii;
        double beta = (Dqii * p1x_a + Eqii * p1y_a + Fqii);
        double gamma = (Aq * p1x_a * p1x_a + Bq * p1x_a * p1y_a + Cq * p1y_a * p1y_a + Dqi * p1x_a + Eqi * p1y_a + Fqi);
        deltaIntersect = -beta / alpha + sqrt(beta * beta / (alpha * alpha) - 2. * gamma / alpha); // delta value when smaller ellipse intersects heavier ellipse "turn on" value
        if (deltaIntersect != deltaIntersect)
        { // in case deltaIntersect is undefined/infinite - can happen if input data is suspect i.e. zero transverse momenta
            deltaIntersect = 0.0;
        }
        if (deltaIntersect <= delta0)
        {
            mt2_unbalanced_value = (ma + mna); // mu_Y "turn on" value of heavier ellipse
            unbalanced_flag = true;
        }
    }
    bool quasiUnbalanced = false;
    if ((ma + mb + mna + mnb) < 0.01)
    { // to capture "quasi-unbalanced" cases as well - these are likely to be very small MT2 values
        //  Do "usual" check for unbalanced massless events (see Lester, C.G. MT2 Special Cases paper, arXiv:1103.5682, for more detail)
        const double eap = (-pax * pmissy + pay * pmissx);
        const double ebp = (-pbx * pmissy + pby * pmissx);
        const double eahbh = sin(atan2(pay, pax) - atan2(pby, pbx)); // note to self: in C++, atan2 is defined as atan2(y-coord, x-coord) i.e. y-coord first!
        if ((eap / eahbh >= 0) && (ebp / eahbh <= 0))
        {
            if (massless_flag)
            {
                unbalanced_flag = true;
            }
            else
            {
                quasiUnbalanced = true;
            }
        }
    }
    // Assuming we do not have an unbalanced situation...
    if (!unbalanced_flag)
    {
        // Next we find intersection mu_Y mass of heavier ellipse with lighter ellipse vertex point thus giving another reasonable guess at an uppermost bound on MT2
        double p1x_b = pmissx;
        double p1y_b = pmissy;
        if (mb > 0)
        {
            p1x_b = pmissx - (mnb / mb) * pbx * Eb / Ea; //x-coordinate of "turn-on" mass of smaller mass ellipse
            p1y_b = pmissy - (mnb / mb) * pby * Eb / Ea; //y-coordinate of "turn-on" mass of smaller mass ellipse
        }
        double deltaIntersectTwo = sqrt((mnasq / Easq) + (p1x_b * p1x_b / Easq) + (p1y_b * p1y_b / Easq)) - (pax * p1x_b / Easq) - (pay * p1y_b / Easq);
        double deltaMax = 0;
        if (massless_flag)
        {
            deltaMax = deltaIntersectTwo; // the kinematic uppermost bound for MT2 in the parabolic case
            disc0 = disc4;
            disc1 = disc5;
            disc2 = disc6;
            disc3 = disc7;
            disc4 = disc8; // discriminant is a quartic so reassign coefficients to their correct place i.e. x^8 coeff is really a x^4 coeff
            disc8 = 0;
            disc7 = 0;
            disc6 = 0;
            disc5 = 0;
        }
        else
        {
            deltaMax = std::min(deltaIntersect, deltaIntersectTwo); // the kinematic uppermost bound for MT2 in elliptical case
        }
        DiscriminantCoeffs discPolynomial = {disc8, disc7, disc6, disc5, disc4, disc3, disc2, disc1, disc0}; //finally can unambiguously set the DiscriminantCoeffs structure

        if (deltaMax != deltaMax)
        {
            deltaMax = 5.0; // highly unlikely but just in case we have an improbable parabola situation where both intersection points are at infinity and deltaMax = nan (5.0 is just some arbitrary positive number)
        }
        int g_bisectDivisor = 2;                                                        // this will be used in Regula Falsi function if the Newton-Raphson method does not find a valid root
        int g_bisectMaxLoops = 50;                                                      // this will be used in Regula Falsi function if the Newton-Raphson method does not find a valid root
        CubicCoeffs cubicPolynomial = {a2, a1, a0, b2, b1, b0, c2, c1, c0, d2, d1, d0}; //set lambda cubic coeffs
        if (!quasiUnbalanced)
        {
            //  Use Newton-Raphson first (for speed) - if it doesn't work we'll use modified Regula Falsi (for safety) method later to obtain root of discriminant
            delta = NewtonRootFinder(delta0, deltaMax, discPolynomial, cubicPolynomial, precision);
        }
        else
        {
            delta = deltaMax;
            g_bisectDivisor = 10;  // probably a near-zero MT2 value, so let's find that zero value quickly!
            g_bisectMaxLoops = 15; // if dividing by 10 each time, only need 15 iterations for machine precision
        }

        /*
    If the ellipses can intersect on their "far side" as well as at their initial periphery then the discriminant can have three (possibly more) positive roots. We want the lowest positive root for MT2.
    There is a chance therefore that the UB may be higher than the second highest or highest root and thus we have not found the lowest positive root.
    Note this is much more likely in parabolic (ie massless or near massless) situations (it is unlikely, but not impossible, for it to happen in elliptical situations)

	Easiest way is to check is if the lambda function has only one positive root (means an ellipse intersection on the periphery) so check this first by counting no of sign changes in function
    (note three positive roots in lambda means an intersection of one conic "inside" i.e on the "far side of" another conic)
*/

        int nooflambdasgnchanges = lambdaSgnchanges(delta, cubicPolynomial);

        if ((nooflambdasgnchanges > 1) || (delta == deltaMax))
        { // last condition ensures a new root is searched for if Newton iteration did not find a valid root

            // If so then we are at a higher root and need to reduce down to the lowest positive root

            delta = NewDeltaFinder(delta0, delta, g_bisectDivisor, g_bisectMaxLoops, discPolynomial, cubicPolynomial, precision);
        }

        MT2 = sqrt(2. * delta * Easq + masq + mnasq);
    }
    if (unbalanced_flag)
    {
        MT2 = mt2_unbalanced_value;
    }
    if (MT2 != MT2)
    { // in case MT2 is undefined/infinite - would only happen if input data is suspect i.e. all zero values etc.
        MT2 = 0.0;
    }
    return MT2;
}

double NewtonRootFinder(double LB, double UB, const DiscriminantCoeffs &discCoeffs, const CubicCoeffs &cubeCoeffs, double accuracy)
{
    int maxIterations = 45; // for vast majority of events should find the root well before this, typically reach required accuracy after 10 iterations  - 45 has been found to be a reasonable number of iterations before N-R method should be abandoned.
    bool solutionFound = false;
    bool outsideLB = false;
    bool outsideUB = false;
    const double UBsq = UB * UB;
    const double UBsqsq = UBsq * UBsq;
    double y_UB = discCoeffs.Coeffs8 * UBsqsq * UBsqsq + discCoeffs.Coeffs7 * UBsqsq * UBsq * UB + discCoeffs.Coeffs6 * UBsqsq * UBsq + discCoeffs.Coeffs5 * UBsqsq * UB +
                  discCoeffs.Coeffs4 * UBsqsq + discCoeffs.Coeffs3 * UBsq * UB + discCoeffs.Coeffs2 * UBsq + discCoeffs.Coeffs1 * UB + discCoeffs.Coeffs0; // value of discriminant at lower bound

    double xNR = UB; // starting point for root value in Newton-Raphson method
    double x1 = UB;
    double originalUB = UB;
    double y_Newt = y_UB;
    double yprime = 0.0;
    double stuckInLoopArray[5] = {-99.0, -98.0, -97.0, -96.0, -95.0}; //just some random negative numbers that delta would never be.
    bool stuckInLoop = false;
    int k = 1;
    for (; k < maxIterations; k++)
    {
        double xNRsq = xNR * xNR;
        double xNRsqsq = xNRsq * xNRsq;
        yprime = 8 * discCoeffs.Coeffs8 * xNRsqsq * xNRsq * xNR + 7 * discCoeffs.Coeffs7 * xNRsqsq * xNRsq + 6 * discCoeffs.Coeffs6 * xNRsqsq * xNR +
                 5 * discCoeffs.Coeffs5 * xNRsqsq + 4 * discCoeffs.Coeffs4 * xNRsq * xNR + 3 * discCoeffs.Coeffs3 * xNRsq + 2 * discCoeffs.Coeffs2 * xNR + discCoeffs.Coeffs1; // value of derivative of discriminant at xNR

        if (!(fabs(yprime) > 0))
        {
            //                // Divide by zero error shouldn't happen - will have to resort to Regula Falsi method
            solutionFound = true;
            x1 = originalUB;
            break;
        }
        x1 = xNR - y_Newt / yprime; // Newton's computation
        if (!(outsideLB && (x1 < LB)) && !(k > maxIterations) && !(outsideUB && (x1 > UB)))
        { // Check if Newton's method still OK i.e. is staying within the original bounds and not heading towards a different root
            if (x1 < LB)
            {
                outsideLB = true; //flagging possible failure of Newton's method
                outsideUB = false;
            }
            else if (x1 > UB)
            {
                outsideUB = true; //flagging possible failure of Newton's method
                outsideLB = false;
            }
            stuckInLoopArray[4] = stuckInLoopArray[3]; //keeping track of last five values for x1 in case we get into a cycle - not common but still a well known issue with the NR method
            stuckInLoopArray[3] = stuckInLoopArray[2];
            stuckInLoopArray[2] = stuckInLoopArray[1];
            stuckInLoopArray[1] = stuckInLoopArray[0];
            stuckInLoopArray[0] = y_Newt / yprime;

            if ((stuckInLoopArray[0] == stuckInLoopArray[4]) || (stuckInLoopArray[0] == stuckInLoopArray[3]) || (stuckInLoopArray[0] == stuckInLoopArray[2]))
            {
                if ((stuckInLoopArray[2] == stuckInLoopArray[0]) && (stuckInLoopArray[4] != stuckInLoopArray[0]))
                {
                }
                else
                {
                    stuckInLoop = true;
                }
            }
            double x1sq = x1 * x1;
            double x1sqsq = x1sq * x1sq;
            y_Newt = discCoeffs.Coeffs8 * x1sqsq * x1sqsq + discCoeffs.Coeffs7 * x1sqsq * x1sq * x1 + discCoeffs.Coeffs6 * x1sqsq * x1sq + discCoeffs.Coeffs5 * x1sqsq * x1 +
                     discCoeffs.Coeffs4 * x1sqsq + discCoeffs.Coeffs3 * x1sq * x1 + discCoeffs.Coeffs2 * x1sq + discCoeffs.Coeffs1 * x1 + discCoeffs.Coeffs0; // value of discriminant at new guess, x1
            if ((((fabs(x1 - xNR) / fabs(x1) < accuracy) || (fabs(y_Newt / yprime) < accuracy)) && (k > 2)) || (stuckInLoop && (stuckInLoopArray[4] != 99.0)))
            { // if the result has met the desired precision && just doublecheck is a root!
                solutionFound = true;
                if (stuckInLoop)
                { //just take the largest of last five values as a new upper bound before returning to main
                    x1 = std::max(x1, std::max(xNR, std::max(xNR + stuckInLoopArray[2], std::max(xNR + stuckInLoopArray[3] + stuckInLoopArray[2], xNR + stuckInLoopArray[4] + stuckInLoopArray[3] + stuckInLoopArray[2]))));
                }
                break;
            }
            xNR = x1;
        }
        else
        {
            break;
        }
    }
    if ((solutionFound == false) || (x1 < 0))
    {
        if (xNR == x1)
        {
            xNR = xNR + stuckInLoopArray[0];
        }
        const double xUB = std::max(xNR + stuckInLoopArray[4], std::max(xNR, x1));
        int nooflambdasgnchanges = lambdaSgnchanges(xUB, cubeCoeffs);
        if ((nooflambdasgnchanges > 1) && (xUB >= 0))
        {
            x1 = std::min(xUB, originalUB); // assume xUB is a reasonable new UB and return to main
        }
        else
        {
            x1 = originalUB;
        }
    }
    return x1;
}

double NewDeltaFinder(double l_delta0, double l_delta, int bisectDivisor, int bisectMaxLoops, const DiscriminantCoeffs &discPolynomial, const CubicCoeffs &cubicPolynomial, double accuracy)
{
    double FunctionVal(double, const DiscriminantCoeffs &);                  // simple function to evaluate f(LB)*f(UB) if negative we know there is a root inbetween
    double RFRootFinder(double, double, const DiscriminantCoeffs &, double); // Regula Falsi function

    // We will search for a root that will leave the lambda function with only one positive root
    double delta2 = l_delta;
    int sgnLoop = 1;
    double deltaMaxOld = l_delta;
    double deltaMaxNew = (deltaMaxOld + l_delta0) / bisectDivisor;
    double l_nooflambdasgnchanges = 2;

    while (sgnLoop <= bisectMaxLoops)
    { //we will reduce the interval until no. of lambda function sgn changes is 1 (or zero) and then we have both a new LB and the previous value is a new UB

        l_nooflambdasgnchanges = lambdaSgnchanges(deltaMaxNew, cubicPolynomial);
        if ((l_nooflambdasgnchanges <= 1))
        {
            sgnLoop++;
            break;
        }
        deltaMaxOld = deltaMaxNew;
        deltaMaxNew = (deltaMaxOld + l_delta0) / bisectDivisor;
        sgnLoop++;
    }
    if ((sgnLoop > 1) && (sgnLoop < (bisectMaxLoops + 1)))
    { // if we have bisected this number of times we are at machine precision, we are at the correct delta value

        // Now we need to check that deltaMaxNew is the new LB ie f(LB)*f(deltaMaxNew) should be > 0 (else it's actually a new UB)
        double rootBounds = 0;
        double newdeltaUB = 0;
        double newdeltaLB = 0;

        rootBounds = FunctionVal(deltaMaxNew, discPolynomial) * FunctionVal(l_delta0, discPolynomial);
        if (rootBounds > 0)
        { // we have the new LB, now need to find the new UB
            newdeltaLB = deltaMaxNew;
            newdeltaUB = deltaMaxOld;
            double checkdelta = 0;
            int counter = 1;
            bool foundnewUB = false;
            // Bisect between LB and UB until we have only the lowest root between them
            while (!foundnewUB)
            {
                checkdelta = (newdeltaUB + newdeltaLB) / 2.0;
                rootBounds = FunctionVal(checkdelta, discPolynomial) * FunctionVal(newdeltaLB, discPolynomial);
                if (rootBounds < 0)
                {
                    //Need to transpose function candidate UB to check there is only one root below the candidate UB value and above the candidate LB (awkward check but couldn't think of more efficient way)
                    const double checkdeltasq = checkdelta * checkdelta;
                    const double checkdeltasqsq = checkdeltasq * checkdeltasq;
                    const double transBdelta = (discPolynomial.Coeffs7 + 8 * checkdelta * discPolynomial.Coeffs8);
                    const double transCdelta = (discPolynomial.Coeffs6 + 28 * checkdeltasq * discPolynomial.Coeffs8 + 7 * checkdelta * discPolynomial.Coeffs7);
                    const double transDdelta = (discPolynomial.Coeffs5 + 56 * checkdeltasq * checkdelta * discPolynomial.Coeffs8 + 21 * checkdeltasq * discPolynomial.Coeffs7 + 6 * checkdelta * discPolynomial.Coeffs6);
                    const double transEdelta = (discPolynomial.Coeffs4 + 70 * checkdeltasq * checkdeltasq * discPolynomial.Coeffs8 + 35 * checkdeltasq * checkdelta * discPolynomial.Coeffs7 + 15 * checkdeltasq * discPolynomial.Coeffs6 + 5 * checkdelta * discPolynomial.Coeffs5);
                    const double transFdelta = (discPolynomial.Coeffs3 + 56 * checkdeltasq * checkdeltasq * checkdelta * discPolynomial.Coeffs8 + 35 * checkdeltasq * checkdeltasq * discPolynomial.Coeffs7 + 20 * checkdeltasq * checkdelta * discPolynomial.Coeffs6 + 10 * checkdeltasq * discPolynomial.Coeffs5 + 4 * checkdelta * discPolynomial.Coeffs4);
                    const double transGdelta = (discPolynomial.Coeffs2 + 28 * checkdeltasq * checkdeltasq * checkdeltasq * discPolynomial.Coeffs8 + 21 * checkdeltasq * checkdeltasq * checkdelta * discPolynomial.Coeffs7 + 15 * checkdeltasq * checkdeltasq * discPolynomial.Coeffs6 + 10 * checkdeltasq * checkdelta * discPolynomial.Coeffs5 + 6 * checkdeltasq * discPolynomial.Coeffs4 + 3 * checkdelta * discPolynomial.Coeffs3);
                    const double transHdelta = (discPolynomial.Coeffs1 + 8 * checkdeltasq * checkdeltasq * checkdeltasq * checkdelta * discPolynomial.Coeffs8 + 7 * checkdeltasq * checkdeltasq * checkdeltasq * discPolynomial.Coeffs7 + 6 * checkdeltasq * checkdeltasq * checkdelta * discPolynomial.Coeffs6 + 5 * checkdeltasq * checkdeltasq * discPolynomial.Coeffs5 + 4 * checkdeltasq * checkdelta * discPolynomial.Coeffs4 + 3 * checkdeltasq * discPolynomial.Coeffs3 + 2 * checkdelta * discPolynomial.Coeffs2);
                    const double transIdelta = (discPolynomial.Coeffs0 + checkdeltasq * checkdeltasq * checkdeltasq * checkdeltasq * discPolynomial.Coeffs8 + checkdeltasq * checkdeltasq * checkdeltasq * checkdelta * discPolynomial.Coeffs7 + checkdeltasq * checkdeltasq * checkdeltasq * discPolynomial.Coeffs6 + checkdeltasq * checkdeltasq * checkdelta * discPolynomial.Coeffs5 + checkdeltasq * checkdeltasq * discPolynomial.Coeffs4 + checkdeltasq * checkdelta * discPolynomial.Coeffs3 + checkdeltasq * discPolynomial.Coeffs2 + checkdelta * discPolynomial.Coeffs1);
                    int nscdelta = 0;
                    if (discPolynomial.Coeffs8 * transBdelta < 0)
                        nscdelta++;
                    if (transBdelta * transCdelta < 0)
                        nscdelta++;
                    if (transCdelta * transDdelta < 0)
                        nscdelta++;
                    if (transDdelta * transEdelta < 0)
                        nscdelta++;
                    if (transEdelta * transFdelta < 0)
                        nscdelta++;
                    if (transFdelta * transGdelta < 0)
                        nscdelta++;
                    if (transGdelta * transHdelta < 0)
                        nscdelta++;
                    if (transHdelta * transIdelta < 0)
                        nscdelta++;

                    const double newdeltaLBsq = newdeltaLB * newdeltaLB;
                    const double newdeltaLBsqsq = newdeltaLBsq * newdeltaLBsq;
                    const double transBLB = (discPolynomial.Coeffs7 + 8 * newdeltaLB * discPolynomial.Coeffs8);
                    const double transCLB = (discPolynomial.Coeffs6 + 28 * newdeltaLBsq * discPolynomial.Coeffs8 + 7 * newdeltaLB * discPolynomial.Coeffs7);
                    const double transDLB = (discPolynomial.Coeffs5 + 56 * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs8 + 21 * newdeltaLBsq * discPolynomial.Coeffs7 + 6 * newdeltaLB * discPolynomial.Coeffs6);
                    const double transELB = (discPolynomial.Coeffs4 + 70 * newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs8 + 35 * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs7 + 15 * newdeltaLBsq * discPolynomial.Coeffs6 + 5 * newdeltaLB * discPolynomial.Coeffs5);
                    const double transFLB = (discPolynomial.Coeffs3 + 56 * newdeltaLBsq * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs8 + 35 * newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs7 + 20 * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs6 + 10 * newdeltaLBsq * discPolynomial.Coeffs5 + 4 * newdeltaLB * discPolynomial.Coeffs4);
                    const double transGLB = (discPolynomial.Coeffs2 + 28 * newdeltaLBsq * newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs8 + 21 * newdeltaLBsq * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs7 + 15 * newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs6 + 10 * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs5 + 6 * newdeltaLBsq * discPolynomial.Coeffs4 + 3 * newdeltaLB * discPolynomial.Coeffs3);
                    const double transHLB = (discPolynomial.Coeffs1 + 8 * newdeltaLBsq * newdeltaLBsq * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs8 + 7 * newdeltaLBsq * newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs7 + 6 * newdeltaLBsq * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs6 + 5 * newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs5 + 4 * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs4 + 3 * newdeltaLBsq * discPolynomial.Coeffs3 + 2 * newdeltaLB * discPolynomial.Coeffs2);
                    const double transILB = (discPolynomial.Coeffs0 + newdeltaLBsq * newdeltaLBsq * newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs8 + newdeltaLBsq * newdeltaLBsq * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs7 + newdeltaLBsq * newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs6 + newdeltaLBsq * newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs5 + newdeltaLBsq * newdeltaLBsq * discPolynomial.Coeffs4 + newdeltaLBsq * newdeltaLB * discPolynomial.Coeffs3 + newdeltaLBsq * discPolynomial.Coeffs2 + newdeltaLB * discPolynomial.Coeffs1);
                    int nscLB = 0;
                    if (discPolynomial.Coeffs8 * transBLB < 0)
                        nscLB++;
                    if (transBLB * transCLB < 0)
                        nscLB++;
                    if (transCLB * transDLB < 0)
                        nscLB++;
                    if (transDLB * transELB < 0)
                        nscLB++;
                    if (transELB * transFLB < 0)
                        nscLB++;
                    if (transFLB * transGLB < 0)
                        nscLB++;
                    if (transGLB * transHLB < 0)
                        nscLB++;
                    if (transHLB * transILB < 0)
                        nscLB++;
                    if (nscLB - nscdelta == 1)
                    { //ie only one root between LB and UB
                        foundnewUB = true;
                    }
                    newdeltaUB = checkdelta;
                }
                else if (lambdaSgnchanges(checkdelta, cubicPolynomial) > 1)
                {
                    newdeltaUB = checkdelta;
                }
                else
                {
                    newdeltaLB = checkdelta;
                }
                counter++;
                // want to avoid possibility of infinite loop
                if (counter > 50)
                {
                    foundnewUB = true; //  this highly unlikely - means two roots within 1E-15 of each other
                    newdeltaUB = checkdelta;
                }
            }
        }
        else
        {
            // If massless or near massless, and the initial LB is zero, then the RF method may not work efficiently (as the LB will generally be stuck near zero for many iterations)
            // We want to therefore find a better LB by decreasing the distance between it and the new UB whilst still keeping the root between them
            newdeltaUB = deltaMaxNew;
            bool foundNewLB = false;
            int LBcounter = 1;
            double checkLB = l_delta0;
            double checkLBprevious = 0;
            while (!foundNewLB)
            {
                checkLBprevious = checkLB;
                checkLB = checkLB + (newdeltaUB - l_delta0) / 5;

                rootBounds = FunctionVal(checkLB, discPolynomial) * FunctionVal(newdeltaUB, discPolynomial);

                if (rootBounds > 0)
                {
                    foundNewLB = true;
                }
                LBcounter++;
                if (LBcounter == 5)
                {
                    foundNewLB = true;
                }
            }
            newdeltaLB = checkLBprevious;
        }
        delta2 = RFRootFinder(newdeltaLB, newdeltaUB, discPolynomial, accuracy);
    }
    else if (sgnLoop == (bisectMaxLoops + 1))
    {
        delta2 = deltaMaxNew; //will be of order zero as of order massless and unbalanced
    }
    return delta2;
}

// Regula Falsi Method for root finding - used if original (NR) interation did not find correct root
double RFRootFinder(double LB, double UB, const DiscriminantCoeffs &discCoeffs, double accuracy)
{
    double x0 = UB;             // starting guess for root value in RF method (using the kinematic lower bound)
    int maxIterationsRF = 1000; // second time around we definitely have bounds solvable by RF, so give it whatever time it needs (within reason!)
    double adjustvalue = 0.0;
    const double LBsq = LB * LB;
    const double LBsqsq = LBsq * LBsq;
    const double UBsq = UB * UB;
    const double UBsqsq = UBsq * UBsq;
    double y_LB = discCoeffs.Coeffs8 * LBsqsq * LBsqsq + discCoeffs.Coeffs7 * LBsqsq * LBsq * LB + discCoeffs.Coeffs6 * LBsqsq * LBsq + discCoeffs.Coeffs5 * LBsqsq * LB +
                  discCoeffs.Coeffs4 * LBsqsq + discCoeffs.Coeffs3 * LBsq * LB + discCoeffs.Coeffs2 * LBsq + discCoeffs.Coeffs1 * LB + discCoeffs.Coeffs0; // value of discriminant at upper bound - multiplication faster than using pow function?
    double y_UB = discCoeffs.Coeffs8 * UBsqsq * UBsqsq + discCoeffs.Coeffs7 * UBsqsq * UBsq * UB + discCoeffs.Coeffs6 * UBsqsq * UBsq + discCoeffs.Coeffs5 * UBsqsq * UB +
                  discCoeffs.Coeffs4 * UBsqsq + discCoeffs.Coeffs3 * UBsq * UB + discCoeffs.Coeffs2 * UBsq + discCoeffs.Coeffs1 * UB + discCoeffs.Coeffs0; // value of discriminant at lower bound
    double y_x0 = y_UB;
    int l = 1;
    for (; l < maxIterationsRF; l++)
    {
        if (fabs(y_UB - y_LB) > 0)
        {
            adjustvalue = -y_x0 * (UB - LB) / (y_UB - y_LB); // RF computation
            x0 = x0 + adjustvalue;
            double x0sq = x0 * x0;
            double x0sqsq = x0sq * x0sq;
            y_x0 = discCoeffs.Coeffs8 * x0sqsq * x0sqsq + discCoeffs.Coeffs7 * x0sqsq * x0sq * x0 + discCoeffs.Coeffs6 * x0sqsq * x0sq + discCoeffs.Coeffs5 * x0sqsq * x0 +
                   discCoeffs.Coeffs4 * x0sqsq + discCoeffs.Coeffs3 * x0sq * x0 + discCoeffs.Coeffs2 * x0sq + discCoeffs.Coeffs1 * x0 + discCoeffs.Coeffs0;
        }
        else
        {
            adjustvalue = 0.0; // assume have found root
            LB = UB;
        }
        if (((fabs(adjustvalue) < accuracy) && (fabs((UB / LB) - 1) < 0.01)) || (l == maxIterationsRF))
        { // if the result has met the desired precision and check it looks like a root! (or just way out of time!)
            break;
        }

        //      Tighten running bounds on root
        if ((y_x0 * y_UB) < 0)
        {
            LB = UB;
            y_LB = y_UB;
        }
        else if (y_x0 != 0)
        {                                       // should always be true unless found exact root
            y_LB = y_LB * y_UB / (y_UB + y_x0); // RF iteration modified using pegasus method

            /*
    This is a more complicated alternative to the Pegasus method, but was found to not be any more efficient
            if((1 - (y_x0/y_UB)) >= 0) { // RF iteration modified using A&B method
				y_LB = (1 - (y_x0/y_UB));
			} else {
                y_LB = 0.5;
            }
*/
        }
        else
        {
            y_LB = 0; // assume have found root
        }
        UB = x0;
        y_UB = y_x0;
    }
    return x0;
}

int lambdaSgnchanges(double deltaFunc, const CubicCoeffs &cubeCoeffs)
{
    const double l3 = cubeCoeffs.Coeffa2 * deltaFunc * deltaFunc + cubeCoeffs.Coeffa1 * deltaFunc + cubeCoeffs.Coeffa0;
    const double l2 = cubeCoeffs.Coeffb2 * deltaFunc * deltaFunc + cubeCoeffs.Coeffb1 * deltaFunc + cubeCoeffs.Coeffb0;
    const double l1 = cubeCoeffs.Coeffc2 * deltaFunc * deltaFunc + cubeCoeffs.Coeffc1 * deltaFunc + cubeCoeffs.Coeffc0;
    const double l0 = cubeCoeffs.Coeffd2 * deltaFunc * deltaFunc + cubeCoeffs.Coeffd1 * deltaFunc + cubeCoeffs.Coeffd0;

    int nsc = 0;
    if (l3 * l2 < 0)
        nsc++;
    if (l2 * l1 < 0)
        nsc++;
    if (l1 * l0 < 0)
        nsc++;
    return nsc;
}

double FunctionVal(double LB, const DiscriminantCoeffs &discCoeffs)
{
    const double LBsq = LB * LB;
    const double LBsqsq = LBsq * LBsq;
    const double funcValue = discCoeffs.Coeffs8 * LBsqsq * LBsqsq + discCoeffs.Coeffs7 * LBsqsq * LBsq * LB + discCoeffs.Coeffs6 * LBsqsq * LBsq + discCoeffs.Coeffs5 * LBsqsq * LB +
                             discCoeffs.Coeffs4 * LBsqsq + discCoeffs.Coeffs3 * LBsq * LB + discCoeffs.Coeffs2 * LBsq + discCoeffs.Coeffs1 * LB + discCoeffs.Coeffs0;
    return funcValue;
}
