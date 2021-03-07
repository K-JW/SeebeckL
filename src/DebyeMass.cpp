/*******************************************************************************
 *
 *  Author: KANG, Jin-Wen
 *  E-Mail: kangjinwen@vip.qq.com
 *  Date:   2020-05-07 14:42 +0800
 *  Desc:   Debye 屏蔽质量及相关.
 * 
 ******************************************************************************/

#include "DebyeMass.hpp"
#include "Distribution.hpp"
#include "Constants.hpp"

#include <cmath>
using std::log;
using std::abs;
using std::pow;

using namespace SeebeckL;

double SeebeckL::alpha_sB(double eB)
{
    return 1 / (1 / alpha_s0 + 11 * Nc / (12 * M_PI) * 
        log((pow(LambdaQCD, 2) + pow(MB, 2)) / pow(mu0, 2)) + 1 / (3 * M_PI * sigma)
        * 2 * (abs(qu * eB) + abs(qd * eB) + abs(qs * eB)));
}

static double sum_dist_for_M_DebyeSqure(double T, double eB, double mu, double q, int l)
{
    double value = 0.0;
    for(int n = 0; n <=l; n++)
    {
        value += (2 - delta0l(n)) * (
            1 / (exp((sqrt(2 * n * abs(q * eB)) + mu) / T) + 1) +
            1 / (exp((sqrt(2 * n * abs(q * eB)) - mu) / T) + 1)
        );
    }
    return value;
}

double SeebeckL::M_DebyeSqure(double T, double eB, double mu, int l)
{
    double part1 = 4 * M_PI * alpha_sB(eB);
    double part2 = pow(T, 2) * Nc / 3;
    double part3 = abs(qu * eB) * sum_dist_for_M_DebyeSqure(T, eB, mu, qu, l) +
        abs(qd * eB) * sum_dist_for_M_DebyeSqure(T, eB, mu, qd, l) + 
        abs(qs * eB) * sum_dist_for_M_DebyeSqure(T, eB, mu, qs, l);
    return part1 * part2 + part1 * 2 * part3 / (4 * M_PI * M_PI);
}

static double A1(double pz, double q, int l, double T, double eB, double m, double mu)
{
    return pz / energy(pz, m, q, eB, l) * (dist0_quark(T, m, mu, pz, q, eB, l) * 
        (1 - dist0_quark(T, m, mu, pz, q, eB, l)) + dist0_antiquark(T, m, mu, pz, q, eB, l) *
        (1 - dist0_antiquark(T, m, mu, pz, q, eB, l)));
}

static double A2(double pz, double q, int l, double T, double eB, double m, double mu, double xi)
{
    return xi * pow(pz, 3) / (2 * pow(energy(pz, m, q, eB, l), 3)) * (
        dist0_quark(T, m, mu, pz, q, eB, l) * (1 - dist0_quark(T, m, mu, pz, q, eB, l)) +
        dist0_antiquark(T, m, mu, pz, q, eB, l) * (1 - dist0_antiquark(T, m, mu, pz, q, eB, l))
    );
}

static double A3(double pz, double q, int l, double T, double eB, double m, double mu, double xi)
{
    return xi * pow(pz, 3) / (2 * T * pow(energy(pz, m, q, eB, l), 2)) * (
        dist0_quark(T, m, mu, pz, q, eB, l) * (1 - dist0_quark(T, m, mu, pz, q, eB, l)) *
        (1 - 2 * dist0_quark(T, m, mu, pz, q, eB, l)) + dist0_antiquark(T, m, mu, pz, q, eB, l) * (
            1 - dist0_antiquark(T, m, mu, pz, q, eB, l)
        ) * (1 - 2 * dist0_antiquark(T, m, mu, pz, q, eB, l))
    );
}

// l>=0
double SeebeckL::M_DebyeBSqure(double T, double eB, double m, double mu, double xi, int l_prime)
{
    //
}