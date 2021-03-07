/*******************************************************************************
 *
 *  Author: KANG, Jin-Wen
 *  E-Mail: kangjinwen@vip.qq.com
 *  Date:   2020-05-07 14:35 +0800
 *  Desc:   Debye 屏蔽质量及相关.
 * 
 ******************************************************************************/

namespace SeebeckL
{
    double alpha_sB(double eB);

    double M_DebyeSqure(double T, double eB, double mu, int l);

    double M_DebyeBSqure(double T, double eB, double m, double mu, double xi, int l_prime);

    double alpha_effB(double T, double eB, double m, double xi, int l);

}