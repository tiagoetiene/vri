#include <iostream>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <tr1/random>
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp>

#include "solutions.h"
#include "integration.h"
#include "pre_integration.h"

typedef long double Real;
std::tr1::random_device rd;
std::tr1::subtract_with_carry_01<Real, 48, 10, 24> gen(rd());

namespace po = boost::program_options;

int main(int argc, const char *argv[])
{

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("start", po::value< std::string >()->default_value("0 0 0"), "ray starting point")
        ("end", po::value< std::string >()->default_value("1 0 0"), "ray ending point")
        ("inner", po::value< std::string>()->default_value("RIEMANN"), "inner integral numerical integration method: RIEMANN, MONTE_CARLO, TRAPEZOID, GAUSS_QUADRATURE, GAUSS_QUADRATURE_5, SIMPSON, BOOLE")
        ("outer", po::value< std::string>()->default_value("RIEMANN"), "outer integral numerical integration method: RIEMANN, MONTE_CARLO, TRAPEZOID, GAUSS_QUADRATURE, GAUSS_QUADRATURE_5, SIMPSON, BOOLE")
        ("exp", po::value< std::string>()->default_value("QUADRATIC"), "exponential approximation method: LINEAR, QUADRATIC, CUBIC, QUARTIC, QUINTIC, EXACT")
        ("step-size", po::value< float >()->default_value(0.125E+0), "step size along the parameterized ray. The ray is parameterized by as X = start + delta * (end - start), where delta is the step size")
        ("input", po::value< std::string >(), "input nrrd scalar field")
        ("color", po::value< std::string >(), "input nrrd color transfer function")
        ("transparency", po::value< std::string >(), "input nrrd extinction coefficient")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    std::cerr << "VRI (Volume Rendering Integral Discretization Schemes)"
              << std::endl;

    std::string start_s = vm["start"].as<std::string>();
    std::string end_s = vm["end"].as<std::string>();

    float start_f[3];
    float end_f[3];

    sscanf(start_s.c_str(), "%f %f %f", &start_f[0], &start_f[1], &start_f[2]);
    sscanf(end_s.c_str(), "%f %f %f", &end_f[0], &end_f[1], &end_f[2]);

    Real start[3] = {start_f[0], start_f[1], start_f[2]};
    Real end[3] = {end_f[0], end_f[1], end_f[2]};

    std::cerr << "\t* Ray starting point: "
              << start[0] << " " << start[1] << " " << start[2] << std::endl;

    std::cerr << "\t* Ray end point     : "
              << end[0] << " " << end[1] << " " << end[2] << std::endl;

    std::cerr << "\t* Inner integral numerical discretization method: "
              << vm["inner"].as<std::string>() << std::endl;
    std::cerr << "\t* Outer integral numerical discretization method: "
              << vm["outer"].as<std::string>() << std::endl;
    std::cerr << "\t* Exponential discretization method             : "
              << vm["exp"].as<std::string>() << std::endl;

    bool pre_integrated_test = false;

    std::vector<Real> I;

    if(!pre_integrated_test)
    {
        VRI_solution_00<Real> solve(start, end);

        Method  exp_method = getMethod( vm["exp"].as<std::string>() ),
                inner_method = getMethod( vm["inner"].as<std::string>() ),
                outer_method = getMethod( vm["outer"].as<std::string>() );

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
    }
    else
    {
        Exp_solution_02<Real> solve(start, end);

        //Number of tests to be made
        unsigned N = 8;

        std::cout << std::setprecision(20);

        //Domain size and step size
        Real D = 1.0;
        Real d = 0.125;
        for(unsigned test = 0; test < N; ++test)
        {
            unsigned n = unsigned((D / d) + 1);
            std::cout << 1.0 / (n-1) << " " << std::flush;
            std::cerr << "(" << n-1 << "," << std::flush;

            solve.d = d;
    //        Real sol = solve.sol(D);
    //        Real num = attenuation(solve, d, n);
            Real sol = solve.sol_vri(D);
            Real num = pre_outer(solve, d, n);

            I.push_back(fabs(sol-num));
            d = d * 0.5;

            std::cerr << I[I.size()-1] << ") " << std::flush;
        }
    }

    std::cout << std::endl;
    std::cerr << std::endl;

    for(const Real& err : I)
        std::cout << err << " ";
    std::cout << std::endl;

    return 0;
}

