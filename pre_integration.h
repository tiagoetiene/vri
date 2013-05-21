#ifndef PRE_INTEGRATION_H
#define PRE_INTEGRATION_H

template<typename Real>
inline
Real pre_attenuation(const Solution<Real> &solve, const Real d, const unsigned n)
{
    Real integral = 1.0;
    for(unsigned j = 0; j < n; ++j)
    {
        integral *= solve.attenuation_1st(solve.X(j*d));
//        integral *= solve.attenuation_2nd(solve.X(j*d));
    }
    return integral;
}

template<typename Real>
inline
Real pre_outer(const Solution<Real> &solve,
           const Real d,
           const unsigned n)
{
    Real I = 0.0;

    for(unsigned i = 0; i < n-1; ++i)
    {
        I += solve.emission_1st(solve.X(i*d)) * pre_attenuation(solve, d, i);
//        I += solve.emission_2nd(solve.X(i*d)) * pre_attenuation(solve, d, i);
//        I += solve.emission_2nd_approx(solve.X(i*d)) * pre_attenuation(solve, d, i);
    }

    return I;
}

#endif // PRE_INTEGRATION_H
