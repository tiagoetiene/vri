#include <iostream>
#include <vector>
#include <numeric>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <tr1/random>

#include "solutions.h"

#define _USE_MATH_DEFINES

enum Method
{
    MONTE_CARLO,
    RIEMANN,
    TRAPEZOID,
    LINEAR,
    QUADRATIC,
    CUBIC,
    QUARTIC,
    QUINTIC,
    EXACT,
    GAUSS_QUADRATURE,
    GAUSS_QUADRATURE_5,
    SIMPSON,
    BOOLE
};


template<typename Real>
inline
const std::vector<Real>& inner(const Solution<Real> &solve,
                               Real d, unsigned i, Method method,
                               const std::vector<Real>& pos_array,
                               std::vector<Real>& integrands)
{
    if(method == MONTE_CARLO)
    {
        size_t j = integrands.size();
        while(j < pos_array.size() and pos_array[j] <= i * d)
            integrands.push_back(d * solve.T(pos_array[j++]));

    }
    else if(method == RIEMANN)
    {
        if(i >= 2)
        {
            size_t j = integrands.size();
            if(j < i-1)
                integrands.push_back(solve.T(j * d) * d);
        }
    }
    else if(method == TRAPEZOID)
    {
        size_t j = integrands.size()+1;
        if(j < i+1)
            integrands.push_back((solve.T(j*d) + solve.T((j-1)*d)) * d * 0.5);
    }
    else if(method == SIMPSON)
    {
        size_t j = integrands.size()+1;
        if(i >= 1)
        {
            integrands.push_back((solve.T((j-1)*d) + 4.0 * solve.T((j - 0.5)*d) + solve.T((j+0)*d)) * d /6.0);
        }
    }
    else if(method == GAUSS_QUADRATURE)
    {
        static Real P[] = {-1.0 / sqrtl(3.0), +1.0 / sqrtl(3.0)};
        static Real W[] = {1.0, 1.0};

        if(i >= 1)
        {
            Real a = (i-1)*d;
            Real b = i*d;
            Real int_ab = 0.0;
            for(unsigned j = 0; j < 2; ++j) int_ab += solve.T( 0.5 * (b-a) * P[j] + 0.5 * (a + b)) * W[j];
            integrands.push_back( 0.5*(b-a) * int_ab );
        }

    }
    else if(method == GAUSS_QUADRATURE_5)
    {
        static Real P[] = {-sqrtl(15.0) / 5.0, 0.0, +sqrtl(15.0) / 5.0};
        static Real W[] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

        if(i >= 1)
        {
            Real a = (i-1)*d;
            Real b = i*d;
            Real int_ab = 0.0;
            for(unsigned j = 0; j < 3; ++j) int_ab += solve.T( 0.5 * (b-a) * P[j] + 0.5 * (a + b)) * W[j];
            integrands.push_back( 0.5*(b-a) * int_ab );
        }
    }
    else
        assert(0);

    return integrands;
}

template<typename Real>
inline
Real exponential(const std::vector<Real>& integrands, Real& S, Method method)
{
    static size_t last_size = 0;
    if(integrands.size() == 0 or integrands.size() == last_size)
        return S;

    last_size = integrands.size();
    const Real s = integrands[last_size-1];

    if (method == LINEAR)
        S = S * (1);
    else if (method == QUADRATIC)
        S = S * (1 - s);
    else if (method == CUBIC)
        S = S * (1 - s + 0.5 * s * s);
    else if (method == QUARTIC)
        S = S * (1 - s + 0.5 * s * s - (1.0/6.0) * s * s * s);
    else if (method == QUINTIC)
        S = S * (1 - s + 0.5 * s * s - (1.0/6.0) * s * s * s + (1.0/24.0) * s * s * s * s);
    else if (method == EXACT)
        S = exp(-std::accumulate(integrands.begin(), integrands.end(), Real(0.0)));
    else
        assert(0);

    return S;
}

template<typename Real>
Real outer(const Solution<Real> &solve, Real d, unsigned n, const Method outer_method, const Method inner_method, const Method exp_method,
           const std::vector<Real>& samples)
{
    std::vector<Real> integrands;
    Real alpha = 1.0;
    Real I = 0.0;

    integrands.reserve(n);
    if (outer_method == RIEMANN)
    {
        for(unsigned i = 1; i < n; ++i)
            I += solve.C(i * d) * solve.T(i * d) * d * exponential(inner(solve, d, i, inner_method, samples, integrands), alpha, exp_method);
    }
    else if(outer_method == TRAPEZOID)
    {
        Real A, B;
        A = solve.C(0 * d) * solve.T(0 * d) * exponential(inner(solve, d, 0, inner_method, samples, integrands), alpha, exp_method);
        for(unsigned i = 1; i < n; ++i)
        {
            B = solve.C(i*d) * solve.T(i*d) * exponential(inner(solve, d, i, inner_method, samples, integrands), alpha, exp_method);
            I += (A+B) * d * 0.5;
            A = B;
        }
    }
    else if(outer_method == SIMPSON)
    {
        Real fa, fm, fb;
        unsigned a, b, c;

        a = 0;
        fa = solve.C(a * d) * solve.T(a * d) * exponential(inner(solve, d, a, inner_method, samples, integrands), alpha, exp_method);
        for(unsigned i = 2; i < n; i += 2)
        {
            b = i-1;
            c = i-0;
            fm = solve.C(b * d) * solve.T(b * d) * exponential(inner(solve, d, b, inner_method, samples, integrands), alpha, exp_method);
            fb = solve.C(c * d) * solve.T(c * d) * exponential(inner(solve, d, c, inner_method, samples, integrands), alpha, exp_method);
            I += fa + 4.0 * fm + fb;
            fa = fb;
        }
        I *= (2.0 * d) / 6.0;
    }
    else if(outer_method == BOOLE)
    {
        static const Real W[] = {7.0, 32.0, 12.0, 32.0, 7.0};
        Real f[5];

        f[0] = solve.C(0 * d) * solve.T(0 * d) * exponential(inner(solve, d, 0, inner_method, samples, integrands), alpha, exp_method);
        for(unsigned i = 4; i < n; i += 4)
        {
            for(int k = 3; k >= 0; --k)
            {
                unsigned l = i-k;
                f[4-k] = solve.C(l * d) * solve.T(l * d) * exponential(inner(solve, d, l, inner_method, samples, integrands), alpha, exp_method);
            }

            for(unsigned k = 0; k < 5; ++k)
                I += W[k] * f[k];

            f[0] = f[4];
        }
        I *= (2.0 * d) / 45.0;
    }
    else
        assert(0);

    return I;
}

typedef long double Real;
std::tr1::random_device rd;
std::tr1::subtract_with_carry_01<Real, 48, 10, 24> gen(rd());

int main()
{
    std::vector<Real> I;
    VRI_solution_00<Real> _solution;
    const Solution<Real> &solve(_solution);

    Method exp_method, inner_method, outer_method;

//    exp_method = LINEAR;
//    exp_method = QUADRATIC;
//    exp_method = CUBIC;
//    exp_method = QUARTIC;
//    exp_method = QUINTIC;
    exp_method = EXACT;

//    inner_method = MONTE_CARLO;
//    inner_method = RIEMANN;
//    inner_method = TRAPEZOID;
    inner_method = SIMPSON;
//    inner_method = GAUSS_QUADRATURE;
//    inner_method = GAUSS_QUADRATURE_5;

//    outer_method = RIEMANN;
//    outer_method = TRAPEZOID;
    outer_method = SIMPSON;
//    outer_method = BOOLE;

    //Number of tests to be made
    unsigned N = 8;

    std::cout << std::setprecision(20);

    //Domain size and step size
    Real D = 1.0;
    Real d = 0.125;
    std::tr1::uniform_real<> dis(0.0, D);
    for(unsigned test = 0; test < N; ++test)
    {
        unsigned n = unsigned((D / d) + 1);
        std::cout << 1.0 / (n-1) << " " << std::flush;
        std::cerr << "(" << n-1 << "," << std::flush;

        std::vector<Real> samples;
        if(inner_method == MONTE_CARLO)
        {
            // ... Create a new array everytime
            samples.resize(n);
            for(Real& v : samples)
                v = dis(gen);
            std::sort(samples.begin(), samples.end());
        }

        Real sol = 0.0, num = 0.0;
        sol = solve.sol(D);
        num = outer(solve, d, n, outer_method, inner_method, exp_method, samples);

        I.push_back(fabs(sol-num));
        d = d * 0.5;

        std::cerr << I[I.size()-1] << ") " << std::flush;
    }
    std::cout << std::endl;
    std::cerr << std::endl;

    for(const Real& err : I)
        std::cout << err << " ";
    std::cout << std::endl;


    return 0;
}

