#ifndef SOLUTIONS_EXP_H
#define SOLUTIONS_EXP_H

#include <numeric>

template<typename Real>
struct Point {
    Real x, y, z;
};

template<typename Real>
struct Solution
{
    virtual Real sol(Real l) const = 0;
    virtual Point<Real> X(Real l) const = 0;
    virtual Real s(const Point<Real>& x) const = 0;
    virtual Real T(Real l) const = 0;
    virtual Real C(Real l) const = 0;
};

template<typename Real>
struct VRI_solution_00 : public Solution<Real>
{
    inline Real sol(Real l) const
    {
    //    return -exp(-l)*( atan(sin(l)/cos(l))-exp(l)+1 );
        return 2-exp(-0.5*sin(l*l))*(sin(l*l)+2);
    //    return 1-exp(-sin(l))*(sin(l)+1);
    //    return 1.0-exp(-l);
    }

    inline Point<Real> X(Real l) const
    {
        return {l, 0.0, 0.0};
    }

    inline Real s(const Point<Real>& x) const
    {
    //    return tan(x.x);
        return x.x;
        //return x.x;
        //return x.x;
    }

    inline Real T(Real l) const
    {
    //    return 1.0;
        return s(X(l))*cos(s(X(l))*s(X(l)));
    //    return cos(s(X(l)));
    //    return 1.0;
    }

    inline Real C(Real l) const
    {
    //    return atan(s(X(l)));
        return sin(s(X(l))*s(X(l)));
        //return sin(s(X(l)));
        //return 1.0;
    }
};


#endif // SOLUTIONS_H
