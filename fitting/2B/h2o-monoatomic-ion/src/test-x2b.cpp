#include <iostream>

#include "x2b-h2o-ion-v1.h"

#include "fit-utils.h"
#include "training-set-h2o-ion.h"

#define PR(x) std::cout << #x << "= " << (x) << std::endl;

//
// To test all the parts that the fitting code uses.
//

// Print arrays of doubles
void prarr(double* arr, int size, std::string name) {
	for (int i = 0; i < size; i++) {
		std::cout << name << "[" << i << "] = " << arr[i] << std::endl;
	}
}

int main() {

//
// TODO: Test I/O
//

//
// TODO: Test x2b-base functions.
//

	double crds[] = {
		0.000000000e+00,  0.000000000e+00,  1.353710000e-01,
		0.000000000e+00,  7.586280000e-01, -4.555280000e-01,
		0.000000000e+00, -7.586280000e-01, -4.555280000e-01,
		0.000000000e+00,  0.000000000e+00, -8.200000000e+00 };
	
	typedef h2o_ion::x2b_h2o_ion_v1<4> model_type;
	model_type model;

//
// TODO: Test polynomial
//


//
// Test x2b-h2o-cl-v1 functions.
//

	double nlpparams[7];
	double basis[105];

	PR(model.name())
	PR(model_type::poly_type::size)
	PR(model.nparams())
	
	model.get_nonlinear_parameters(nlpparams);
	prarr(nlpparams, 7, "nlpparams");

	//double basis_value = model.basis(crds, basis);
	//PR(basis_value)

	//for (int i = 0; i < 105; i++) {
	//	std::cout << "basis[" << i << "] = " << basis[i] << std::endl;
	//}

	double value = model.value(crds);
	PR(value)

}
