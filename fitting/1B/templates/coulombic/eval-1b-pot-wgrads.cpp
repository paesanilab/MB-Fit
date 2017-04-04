#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include "io-xyz.h"
#include "x1b-v1x.h"

int main(int argc, char** argv) {

    if (argc != 3) {
	std::cerr << "usage: eval-1b-pot params.nc monomer.xyz" << std::endl;
	return 0;
    }

    std::cout << std::scientific << std::setprecision(9);

    h3o::x1b_v1x pot_1b;

    std::vector<std::string> elements;
    std::vector<double> crd;

    try {
	++argv;
	--argc;
	pot_1b.load_poly_dat(*argv);

	++argv;
	--argc;
	std::ifstream ifs(*argv);

	if (!ifs)
	    throw std::runtime_error("could not load the XYZ file");

	std::string comment;
	kit::io::load_xyz(ifs, comment, elements, crd);
    } catch (const std::exception& e) {
	std::cerr << " ** ERROR ** : " << e.what() << std::endl;
	return 1;
    }

    double Enogrd = pot_1b(&(crd[0]));
    std::cout << "Enogrd = " << Enogrd  << std::endl;
    
    double  grd[12];
    double Egrd = pot_1b(&(crd[0]), &(grd[0]));
    std::cout << "Egrd = " << Egrd  << std::endl;
   
    const double eps = 1.0e-4;

    for (int n = 0; n < 12; ++n) {
        double tmp[12];
        const double x_orig = crd[n];

        crd[n] = x_orig + eps;
	const double Ep = pot_1b(&(crd[0]));

        crd[n] = x_orig - eps;
	const double Em = pot_1b(&(crd[0]));

        const double gfd = 0.5*(Ep - Em)/eps;
        crd[n] = x_orig;

	std::cout << "grd[" << n << "] = " << grd[n] << "\t";
        std::cout << "gfd = " << gfd << '\t';
        std::cout << " error = " << std::fabs(grd[n] - gfd) << '\n';
    }

    return 0;
}
