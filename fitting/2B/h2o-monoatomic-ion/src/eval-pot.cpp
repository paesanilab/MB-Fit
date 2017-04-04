#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "io-xyz.h"
#include "xyz-water-utils.h"
#include "x2b-h2o-ion-v1.h"

int main(int argc, char** argv)
{
    if (argc != 3) {
	std::cerr << "usage: eval v1.nc dimer.xyz" << std::endl;
	return 0;
    }

    std::cout << std::scientific << std::setprecision(9);

    h2o_ion::x2b_h2o_ion_v1<5> pot;

    std::vector<std::string> elements;
    std::vector<double> crd;

    try {
	++argv;
       	--argc;
	pot.load_netcdf(*argv);

	++argv;
       	--argc;
	std::ifstream ifs(*argv);

	if(!ifs)
	    throw std::runtime_error("could not open the XYZ file");

	std::string comment;
	kit::io::load_xyz(ifs, comment, elements, crd);
    } catch (const std::exception& e) {
	std::cerr << " ** Error ** : " << e.what() << std::endl;
	return 1;
    }

    double E = pot.value(&(crd[0]));

    std::cout << "E = " << E << std::endl;

    return 0;
}


