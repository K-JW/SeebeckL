/*******************************************************************************
 *
 *  Author: KANG, Jin-Wen
 *  E-Mail: kangjinwen@vip.qq.com
 *  Date:   2020-05-07 14:00 +0800
 *  Desc:   此文件内定义运行过程中用到的常数.
 * 
 ******************************************************************************/
#include <cmath>
// #include <gsl/gsl_math.h>

using std::log;

namespace SeebeckL
{
    const int Nc = 3;
    const double M0u = 0.003, M0d = 0.005, M0s = 0.08;
    const double CR = (Nc * Nc - 1) / (2 * Nc);
    const double MB = 1, sigma = 0.18, LambdaV = 0.385, mu0 = 1.1;
    const double LambdaQCD = 0.2;
    const double alpha_s0 = 12 * M_PI / (11 * Nc * log((mu0 * mu0 + MB * MB) / 
        (LambdaV * LambdaV)));
    const double qu = 2.0/3.0, qd = -1.0/3.0, qs = -1.0/3.0;
    const double qub = -2.0/3.0, qdb = +1.0/3.0, qsb = +1.0/3.0;
}