#ifndef SOLUTIONS_EXP_H
#define SOLUTIONS_EXP_H

#include <numeric>
#include <cmath>

const long double PI = std::atan(1.0)*4.0;

template<typename Real>
struct Point {
    Real x, y, z;
};

template<typename real>
inline
real magnitue(const real *v)
{
    return std::sqrt(std::inner_product(v+0, v+3, v+0, real(0.0)));
}

template<typename Real>
struct Solution
{
    Solution(const Real *start, const Real *end)
    {
        for(unsigned i = 0; i < 3; ++i)
        {
            m_start[i] = start[i];
            m_end[i]   = end[i];
            m_step[i] = end[i]-start[i];
        }
    }

    virtual inline Real sol(Real) const { assert(0); return Real(0.0); }
    virtual inline Real sol_vri(Real) const { assert(0); return Real(0.0); }
    virtual inline Real s(const Point<Real>&) const { assert(0); return Real(0.0); }
    virtual inline Real emission_1st(const Point<Real>&) const { assert(0); return Real(0.0); }
    virtual inline Real emission_2nd(const Point<Real>&) const { assert(0); return Real(0.0); }
    virtual inline Real emission_2nd_approx(const Point<Real>&) const { assert(0); return Real(0.0); }
    virtual inline Real attenuation_1st(const Point<Real>&) const { assert(0); return Real(0.0); }
    virtual inline Real attenuation_2nd(const Point<Real>&) const { assert(0); return Real(0.0); }
    virtual inline Real T(Real) const { assert(0); return Real(0.0); }
    virtual inline Real C(Real) const { assert(0); return Real(0.0); }
    virtual inline Point<Real> X(Real lambda) const
    {
        Point<Real> x;
        x.x = m_start[0] + lambda*m_step[0];
        x.y = m_start[1] + lambda*m_step[1];
        x.z = m_start[2] + lambda*m_step[2];
        return x;
    }

    Real m_start[3];
    Real m_end[3];
    Real m_step[3];
};

template<typename Real>
struct VRI_solution_00 : public Solution<Real>
{
    VRI_solution_00(const Real *start, const Real *end) :
        Solution<Real>::Solution(start, end)
    {
    }
    inline Real sol(Real l) const
    {
        return 2-exp(-0.5*sin(l*l))*(sin(l*l)+2);
    }
    inline Real s(const Point<Real>& x) const
    {
        return x.x;
    }
    inline Real T(Real l) const
    {
        return s(this->X(l))*cos(s(this->X(l))*s(this->X(l)));
    }
    inline Real C(Real l) const
    {
        return sin(s(this->X(l))*s(this->X(l)));
    }
};


template<typename Real>
struct VRI_solution_01 : public Solution<Real>
{
    inline Real sol(Real l) const
    {
        return -exp(-l)*( atan(sin(l)/cos(l))-exp(l)+1 );
    }
    inline Real s(const Point<Real>& x) const
    {
        return tan(x.x);
    }
    inline Real T(Real l) const
    {
        return 1.0;
    }
    inline Real C(Real l) const
    {
        return atan(s(X(l)));
    }
};

template<typename Real>
struct VRI_solution_02 : public Solution<Real>
{
    inline Real sol(Real l) const
    {
        return 1-exp(-sin(l))*(sin(l)+1);
    }

    inline Real s(const Point<Real>& x) const
    {
        return x.x;
    }

    inline Real T(Real l) const
    {
        return cos(s(X(l)));
    }

    inline Real C(Real l) const
    {
        return sin(s(X(l)));
    }
};

template<typename Real>
struct VRI_solution_03 : public Solution<Real>
{
    inline Real sol(Real l) const
    {
        return 1.0-exp(-l);
    }

    inline Real s(const Point<Real>& x) const
    {
        return x.x;
    }

    inline Real T(Real l) const
    {
        return 1.0;
    }

    inline Real C(Real l) const
    {
        return 1.0;
    }
};
template<typename Real>
struct Exp_solution_00 : public Solution<Real>
{
    mutable Real d = std::numeric_limits<Real>::infinity();
    inline Real sol(Real l) const
    {
        return exp(1-exp(l));
    }

    inline Real s(const Point<Real>& x) const
    {
        return exp(x.x);
    }
    inline Real attenuation_1st(const Point<Real>& x) const
    {
        Point<Real> x_n = {0.0, 0.0, 0.0};
        x_n.x = x.x + d;
        Real sf = s(x_n);
        Real sb = s(x);
        return exp(-(d*(sf+sb))/2);
    }
    inline Real attenuation_2nd(const Point<Real>& _x) const
    {
        Point<Real> x_n = {0.0, 0.0, 0.0};
        Point<Real> x_m = {0.0, 0.0, 0.0};
        x_n.x = _x.x + d;
        x_m.x = _x.x + 0.5*d;
        Real sf = s(x_n);
        Real sm = s(x_m);
        Real sb = s(_x);
        Real x = _x.x;

        return exp(((8*sm-4*sf-4*sb)*x*x*x+(12*d*sm-3*d*sf-9*d*sb)*x*x-6*d*d*sb*x-4*d*d*d*sm-d*d*d*sf-d*d*d*sb)/(6*d*d)-((8*sm-4*sf-4*sb)*x*x*x+(12*d*sm-3*d*sf-9*d*sb)*x*x-6*d*d*sb*x)/(6*d*d));
    }
    inline Real T(Real l) const
    {
        return s(X(l));
    }
};

template<typename Real>
struct Exp_solution_02 : public Solution<Real>
{
    Exp_solution_02(const Real *start, const Real *end) :
        Solution<Real>::Solution(start, end)
    {
    }

    mutable Real d = std::numeric_limits<Real>::infinity();
    inline Real sol(Real l) const
    {
        return 1.0 / (l + 1.0);
    }
    inline Real sol_vri(Real D) const
    {
        return log(D+1);
    }
    inline Real s(const Point<Real>& x) const
    {
        return 1.0/(x.x + 1.0);
    }
    inline Real emission_1st(const Point<Real>& x) const
    {
        Point<Real> x_n = {0.0, 0.0, 0.0};
        x_n.x = x.x + d;
        Real sf = s(x_n);
        Real sb = s(x);

        return (sqrt(PI)*sqrt(d*sb-d*sf)*exp(-(d*sf*sf)/(2*sf-2*sb))*(sqrt(2)*erf((d*sf)/(sqrt(2)*sqrt(d*sb-d*sf)))-sqrt(2)*erf((d*sb)/(sqrt(2)*sqrt(d*sb-d*sf)))))/(2*sf-2*sb);
    }
    inline Real f(Real sb, Real sm, Real sf, Real x, Real t) const
    {
        return exp(-(4*sm*x*x*x)/(3*d*d)+(2*sf*x*x*x)/(3*d*d)+(2*sb*x*x*x)/(3*d*d)+(4*sm*t*x*x)/(d*d)-(2*sf*t*x*x)/(d*d)-(2*sb*t*x*x)/(d*d)-(2*sm*x*x)/d+(sf*x*x)/(2*d)+(3*sb*x*x)/(2*d)-(4*sm*t*t*x)/(d*d)+(2*sf*t*t*x)/(d*d)+(2*sb*t*t*x)/(d*d)+(4*sm*t*x)/d-(sf*t*x)/d-(3*sb*t*x)/d+sb*x+(4*sm*t*t*t)/(3*d*d)-(2*sf*t*t*t)/(3*d*d)-(2*sb*t*t*t)/(3*d*d)-(2*sm*t*t)/d+(sf*t*t)/(2*d)+(3*sb*t*t)/(2*d)-sb*t);
    }
    inline Real integrate(Real sb, Real sm, Real sf, Real x) const
    {
        unsigned N = 1001;
        Real delta = d / (N-1);
        Real integral = 0.0;
        for(unsigned i = 0; i < N-1; ++i)
            integral += 0.5 * delta * (f(sb, sm, sf, x, x + i*delta) + f(sb, sm, sf, x, x + (i+1)*delta));
        return integral;
    }
    inline Real emission_2nd(const Point<Real>& x) const
    {
        Point<Real> x_n = {0.0, 0.0, 0.0};
        Point<Real> x_m = {0.0, 0.0, 0.0};
        x_n.x = x.x + d;
        x_m.x = x.x + 0.5*d;
        Real sf = s(x_n);
        Real sm = s(x_m);
        Real sb = s(x);

        return integrate(sb, sm, sf, x.x);
    }
    inline Real emission_2nd_approx(const Point<Real>& /*x*/) const
    {
        return d;
    }
    inline Real attenuation_1st(const Point<Real>& x) const
    {
        Point<Real> x_n = {0.0, 0.0, 0.0};
        x_n.x = x.x + d;
        Real sf = s(x_n);
        Real sb = s(x);
        return exp(-(d*(sf+sb))/2);
    }
    inline Real attenuation_2nd(const Point<Real>& _x) const
    {
        Point<Real> x_n = {0.0, 0.0, 0.0};
        Point<Real> x_m = {0.0, 0.0, 0.0};
        x_n.x = _x.x + d;
        x_m.x = _x.x + (2.0*d)/3.0;
//        x_m.x = _x.x + 0.5*d;
        Real sf = s(x_n);
        Real sm = s(x_m);
        Real sb = s(_x);
        Real x = _x.x;

        return exp(((6*sm-4*sf-2*sb)*x*x*x+(9*d*sm-4*d*sf-5*d*sb)*x*x-4*d*d*sb*x-3*d*d*d*sm-d*d*d*sb)/(4*d*d)-((6*sm-4*sf-2*sb)*x*x*x+(9*d*sm-4*d*sf-5*d*sb)*x*x-4*d*d*sb*x)/(4*d*d));
//        return exp(((8*sm-4*sf-4*sb)*x*x*x+(12*d*sm-3*d*sf-9*d*sb)*x*x-6*d*d*sb*x-4*d*d*d*sm-d*d*d*sf-d*d*d*sb)/(6*d*d)-((8*sm-4*sf-4*sb)*x*x*x+(12*d*sm-3*d*sf-9*d*sb)*x*x-6*d*d*sb*x)/(6*d*d));
    }
    inline Real T(Real l) const
    {
        return s(this->X(l));
    }
    inline Real C(Real l) const
    {
        return 1./s(this->X(l));
    }
};

#endif // SOLUTIONS_H
