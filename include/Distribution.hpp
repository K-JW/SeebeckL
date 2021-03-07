/*******************************************************************************
 *
 *  Author: KANG, Jin-Wen
 *  E-Mail: kangjinwen@vip.qq.com
 *  Date:   2020-05-07 14:14 +0800
 *  Desc:   定义分布函数.
 * 
 ******************************************************************************/

#include <cmath>

using std::sqrt;
using std::pow;
using std::abs;
using std::exp;

namespace SeebeckL
{   
    // 定义 delta_{0l} 函数
    int delta0l(int l) { if(l == 0) return 1; else return 0; }
    
    // 定义能量计算函数
    double energy(double p, double m, double q, double eB, int l)
    {
        return sqrt(pow(p, 2) + pow(m, 2) + 2 * l * abs(q * eB));
    }

    // 定义 f_g^0 
    double dist0_gluon(double T, double e) { return 1 / (exp(e / T) - 1); }

    // f_q^{0}
    double dist0_quark(double T, double m, double mu, double p, double q, double eB, int l)
    {
        return 1 / (exp((energy(p, m, q, eB, l) - mu) / T) + 1);
    }

    // f_qbar^0
    double dist0_antiquark(double T, double m, double mu, double p, double q, double eB, int l)
    {
        return 1 / (exp((energy(p, m, q, eB, l) + mu) / T) + 1);
    }

    // f_q^{\xi}
    double dist_quark(double T, double m, double mu, double p, double xi, 
        double q, double eB, int l)
    {
        return dist0_quark(T, m, mu, p, q, eB, l) - xi * pow(p, 2) / (2 * 
            energy(p, m, q, eB, l) * T) * dist0_quark(T, m, mu, p, q, eB, l) * (
                1 - dist0_quark(T, m, mu, p, q, eB, l)
            );
    }

    // f_qbar^{\xi}
    double dist_antiquark(double T, double m, double mu, double p, double xi, 
        double q, double eB, int l)
    {
        return dist0_antiquark(T, m, mu, p, q, eB, l) - xi * pow(p, 2) / (2 *
            energy(p, m, q, eB, l) * T) * dist0_antiquark(T, m, mu, p, q, eB, l) * (
                1 - dist0_antiquark(T, m, mu, p, q, eB, l)
            );
    }

}