#ifndef INTEGRATION_H
#define INTEGRATION_H

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

inline
Method getMethod(const std::string &m)
{
    if(m == "MONTE_CARLO")          return MONTE_CARLO;
    if(m == "RIEMANN")              return RIEMANN;
    if(m == "TRAPEZOID")            return TRAPEZOID;
    if(m == "LINEAR")               return LINEAR;
    if(m == "QUADRATIC")            return QUADRATIC;
    if(m == "CUBIC")                return CUBIC;
    if(m == "QUARTIC")              return QUARTIC;
    if(m == "QUINTIC")              return QUINTIC;
    if(m == "EXACT")                return EXACT;
    if(m == "GAUSS_QUADRATURE")     return GAUSS_QUADRATURE;
    if(m == "GAUSS_QUADRATURE_5")   return GAUSS_QUADRATURE_5;
    if(m == "SIMPSON")              return SIMPSON;
    if(m == "BOOLE")                return BOOLE;
    assert(0 and "Integration method not found");
    return Method(0);
}

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

#endif // INTEGRATION_H
