#include <iostream>
#include <vector>
#include "x3b-h2o-ion-v1.h"
#include "training-set-h2o-ion.h"

namespace {


} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: ./test-x3b ts1"
                  << std::endl;
        return 0;
    }

    std::vector<h2o_ion::trimer> training_set;
    typedef h2o_ion::x3b_h2o_ion_v1 model_type;
	static model_type model;

    size_t ntrimers = h2o_ion::load_trimers(argv[1], training_set);
    
    for(int i = 0; i < 21; i++) {
	std::cout << "training_set[0].xyz[" << i << "] = " << training_set[0].xyz[i] << "\n";
    }

	double value = model.value(training_set[0].xyz);
	std::cout << "value = " << value << std::endl;

    return 0;
}
