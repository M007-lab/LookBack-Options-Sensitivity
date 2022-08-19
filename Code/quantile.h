//
// Created by toor on 4/19/21.
//

#ifndef LOOKBACKOPTION_QUANTILE_H
#define LOOKBACKOPTION_QUANTILE_H

#include <cmath>
#include <iostream>
#include <random>

static double a[] = {
        2.50662823884, -18.61500062529,
        41.39119773534, -25.44106049637 };
static double b[] = {
        -8.47351093090, 23.08336743743,
        -21.06224101826, 3.13082909833 };
static double c[] = {
        0.3374754822726147, 0.9761690190917186,
        0.1607979714918209, 0.0276438810333863,
        0.0038405729373609, 0.0003951896511919,
        0.0000321767881768, 0.0000002888167364,
        0.0000003960315187 };


double normal_cdf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

double quantile_normal(double u) {
    double y = u - 0.5;
    if (fabs(y) < 0.42) {
        double r = y*y, nume = a[3], denom = b[3];
        for (int i = 2; i >= 0; i--) {
            nume *= r; nume += a[i];
            denom *= r; denom += b[i];
        }
        return y * nume / (denom * r  + 1.);
    }
    else {
        double r = (u > 0.5) ? log(-log(1-u)) : log(-log(u));
        double x = c[8];
        for (int i = 7; i >= 0; i--) {
            x *= r; x += c[i];
        }
        return (y < 0) ? -x : x;
    }
};


struct gauss_cond {
    gauss_cond(double phi_a = 0, double phi_b = 1)
            : phi_a(phi_a), phi_b(phi_b), U(0,1) {}; //, s(0,1) {};
    template <typename TGen>
    double operator()(TGen & gen) {
        return quantile_normal(phi_a + U(gen)*(phi_b - phi_a));
        //return quantile(s, phi_a + U()*(phi_b - phi_a));
    };
private:
    double phi_a, phi_b;
    std::uniform_real_distribution<> U;
    //      boost::math::normal::normal_distribution s;
};

#endif //LOOKBACKOPTION_QUANTILE_H
