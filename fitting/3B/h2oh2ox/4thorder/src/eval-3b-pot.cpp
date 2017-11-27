#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include "io-xyz.h"
#include "x3b-h2o-ion-v1.h"

#ifdef GET_GRADS
#include "x3b-h2o-ion-v1x.h"
#endif

#include "poly-3b-h2o-ion-v1x.h"

int main(int argc, char** argv) {

    if (argc != 3) {
	std::cerr << "usage: eval-3b-pot v1.nc trimer.xyz" << std::endl;
	return 0;
    }

    std::cout << std::scientific << std::setprecision(9);

    h2o_ion::x3b_h2o_ion_v1 pot_3b;

    std::vector<std::string> elements;
    std::vector<double> crd;

    try {
	++argv;
	--argc;
	pot_3b.load_netcdf(*argv);

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

    double E = pot_3b.value(&(crd[0]));
    std::cout << "E = " << E  << std::endl;
    /*h2o_ion::x3b_h2o_ion_v1 pot_3b;

    double terms[924];

    double a[924];
    double ksi[21];
    double g[21];

    for (int i = 0; i < 924; i++) {
	a[i] = i+1;
    }

    for (int i = 0; i < 21; i++) {
	ksi[i] = std::exp(3*i+1);
    }	
    double xvalue = h2o_ion::poly_3b_h2o_ion_v1x::eval(a,ksi,g);

    pot_3b.doevpoly(ksi,terms);
    double fitvalue = 0.0;
    for (int i = 0; i < 924; i++) {
	fitvalue += a[i]*terms[i];
    }


    std::cout << "xvalue = " << xvalue << std::endl;
    std::cout << "fitvalue = " << fitvalue << std::endl;*/
    return 0;
}
