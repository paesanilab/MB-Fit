#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include "io-xyz.h"
#include "x1b-v1.h"

//#ifdef GET_GRADS
//#include "x1b-h2o-cl-v1x.h"
//#endif

//#include "poly-1b-h2o-cl-v1x.h"

int main(int argc, char** argv) {

    if (argc != 3) {
	std::cerr << "usage: eval-1b-pot params.nc monomer.xyz" << std::endl;
	return 0;
    }

    std::cout << std::scientific << std::setprecision(9);

    h3o::x1b_v1 pot_1b;

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

    double E = pot_1b.value(&(crd[0]));
    std::cout << "E = " << E  << std::endl;
    return 0;
}
