from . import utils_nb_fitting 
import itertools as it

def write_monomer_class_header(mon_index):
    """
    Writes the monomer header for the fitting code

    Args:
        mon_index - The index of the monomer you want to write the header for.
    """

    filename = "mon" + str(mon_index) + ".h"

    header_text = """#ifndef MON""" + str(mon_index) + """_H
#define MON""" + str(mon_index) + """_H
#include "molecule.h"
#include "constants.h"

#include <algorithm>
#include <cstddef>

namespace x   {

   class mon""" + str(mon_index) + """ : public molecule {
      public :
      mon""" + str(mon_index) + """();
      ~mon""" + str(mon_index) + """();

      mon""" + str(mon_index) + """(double* crd);
      double* set_sitecrds(double* xyz);
      double* set_charges(double* xyz);
      double* set_polfacs(double* atmpolar);
      double* set_pol();
      void allocate();

      int get_nsites();
      int get_realsites();
      double* get_sitecrds();
      double* get_charges();
      double* get_polfacs();
      double* get_pol();

      excluded_set_type::iterator get_begin_12();
      excluded_set_type::iterator get_begin_13();
      excluded_set_type::iterator get_begin_14();
      excluded_set_type::iterator get_end_12();
      excluded_set_type::iterator get_end_13();
      excluded_set_type::iterator get_end_14();

   };
}
#endif
"""

    with open(filename, 'w') as header_file:
        header_file.write(header_text)


def write_monomer_class_cpp(mon_index, nsites, natoms, excl12, excl13, excl14, chg, pol, polfac):
    """
    Writes the monomer cpp for the fitting code

    Args:
        mon_index - The index of the monomer you want to write the cpp for
        nsites    - Number of sites (electrostatic sites) of the monomer
        natoms    - Number of real atoms in the monomer
        excl12    - List with the pairs that are excluded at 1-2 distance
        excl13    - List with the pairs that are excluded at 1-3 distance
        excl14    - List with the pairs that are excluded at 1-4 distance
        chg       - Charges of each site
        pol       - Polarizability of each site
        polfac    - Polarizability factor of each site
    """

    filename = "mon" + str(mon_index) + ".cpp"

    cpp_text = """#include "mon""" + str(mon_index) + """.h"

namespace x  {

  mon""" + str(mon_index) + """::mon""" + str(mon_index) + """() { }

  mon""" + str(mon_index) + """::~mon""" + str(mon_index) + """() {
    delete[] memory;
  }

  mon""" + str(mon_index) + """::mon""" + str(mon_index) + """( double* crd) {

    nsites = """ + str(nsites) + """;
    realsites = """ + str(natoms) + """;
    is_w = 0;
    allocate();

    sitecrds = set_sitecrds(crd);
    atmpolar = set_pol();
    charge = set_charges(crd);
    polfac = set_polfacs(atmpolar);

    // Excluded pairs
"""

    for p in excl12:
        cpp_text += '    excluded12.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n'
    for p in excl13:
        cpp_text += '    excluded13.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n'
    for p in excl14:
        cpp_text += '    excluded14.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n'
    cpp_text += """

  }

  double* mon""" + str(mon_index) + """::set_charges(double* atmcrds) {
    charge = memory;
"""

    for i in range(len(chg)):
        cpp_text += '    charge[' + str(i) + '] = ' + str(chg[i]) + '*constants::CHARGECON;\n'

    cpp_text += """

    return charge;
  }

  double* mon""" + str(mon_index) + """::set_pol() {
    atmpolar = memory + nsites + nsites*3;
"""

    for i in range(len(pol)):
        cpp_text += '    atmpolar[' + str(i) + '] = ' + str(pol[i]) + ';\n'

    cpp_text += """
    return atmpolar;
  }

  double* mon""" + str(mon_index) + """::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;

"""
    for i in range(len(polfac)):
        cpp_text += '    polfac[' + str(i) + '] = ' + str(polfac[i]) + ';\n'

    cpp_text += """
    return polfac;
  }

  double* mon""" + str(mon_index) + """::set_sitecrds(double* atmcrds) {
    sitecrds = memory + nsites;

    return atmcrds;
  }

  void mon""" + str(mon_index) + """::allocate() {
    memory = new double [nsites  //charge
      + nsites*3  //site coordinates
      + nsites  //polarizabilities
      + nsites];  //polfac
    }

    int mon""" + str(mon_index) + """::get_nsites() { return nsites; }
    int mon""" + str(mon_index) + """::get_realsites() { return realsites; }
    double* mon""" + str(mon_index) + """::get_charges() { return charge; }
    double* mon""" + str(mon_index) + """::get_sitecrds() { return sitecrds; }
    double* mon""" + str(mon_index) + """::get_pol() { return atmpolar; }
    double* mon""" + str(mon_index) + """::get_polfacs() { return polfac; }

    excluded_set_type::iterator mon""" + str(mon_index) + """::get_begin_12() { return excluded12.begin(); }
    excluded_set_type::iterator mon""" + str(mon_index) + """::get_begin_13() { return excluded13.begin(); }
    excluded_set_type::iterator mon""" + str(mon_index) + """::get_begin_14() { return excluded14.begin(); }
    excluded_set_type::iterator mon""" + str(mon_index) + """::get_end_12() { return excluded12.end(); }
    excluded_set_type::iterator mon""" + str(mon_index) + """::get_end_13() { return excluded13.end(); }
    excluded_set_type::iterator mon""" + str(mon_index) + """::get_end_14() { return excluded14.end(); }

} // namespace x
"""

    with open(filename, 'w') as cpp_file:
        cpp_file.write(cpp_text)


def write_mbpol_monomer(mon_index):
    """
    Writes the mb-pol monomer cpp for the fitting code
    NOTE: If used, order of atoms MUST be OHH

    Args:
        mon_index - The index of the monomer you want to write the cpp for.
    """

    filename = "mon" + str(mon_index) + ".cpp"

    cpp_text = """
#include "mon""" + str(mon_index) + """.h"

namespace {

const double gammaM = 0.426706882;
const double gamma1 = 1.0 - gammaM;
const double gamma2 = gammaM/2;

inline void compute_M_site_crd
    (const double O[3], const double H1[3], const double H2[3], double M[3])
{
    for (size_t i = 0; i < 3; ++i)
        M[i] = gamma1*O[i] + gamma2*(H1[i] + H2[i]);
}

} // namespace

namespace x {

mon""" + str(mon_index) + """::mon""" + str(mon_index) + """() { }

mon""" + str(mon_index) + """::~mon""" + str(mon_index) + """() {
  delete[] memory;
}

mon""" + str(mon_index) + """::mon""" + str(mon_index) + """(double* crd) {
    nsites = 4;
    realsites = 3;
    is_w = 1;
    allocate();
    sitecrds = set_sitecrds(crd);
    atmpolar = set_pol();
    charge = set_charges(crd);
    polfac = set_polfacs(atmpolar);

    excluded12.clear();
    excluded13.clear();

    excluded12.insert(std::make_pair(0, 1)); // O - H1
    excluded12.insert(std::make_pair(0, 2)); // O - H2
    excluded12.insert(std::make_pair(0, 3)); // O - M
    excluded13.insert(std::make_pair(1, 2)); // H1 - H2
    excluded12.insert(std::make_pair(1, 3)); // H1 - M
    excluded12.insert(std::make_pair(2, 3)); // H2 - M
}

double* mon""" + str(mon_index) + """::set_sitecrds(double* atmcrds) {
    sitecrds = memory + nsites;
    // assumes O H H
    compute_M_site_crd(atmcrds, atmcrds + 3, atmcrds + 6, sitecrds + 9);
    std::copy(atmcrds, atmcrds + 9, sitecrds);
    return sitecrds;
}

double* mon""" + str(mon_index) + """::set_charges(double* atmcrds) {
    double chgtmp[3];
    h2o::ps::dms_nasa(0.0, 0.0, 0.0, atmcrds, chgtmp, 0, false);
    const double tmp = 0.5*gammaM/(1.0 - gammaM);
    charge = memory;
    charge[0] = 0.0;                        // O
    charge[1] = constants::CHARGECON*(chgtmp[1] + tmp*(chgtmp[1] + chgtmp[2])); // H1
    charge[2] = constants::CHARGECON*(chgtmp[2] + tmp*(chgtmp[1] + chgtmp[2])); // H2
    charge[3] = constants::CHARGECON*(chgtmp[0]/(1.0 - gammaM));       // M

    //std::cout << "charge[1] = " << charge[1] << std::endl;

    return charge;
}

double* mon""" + str(mon_index) + """::set_pol() {
    atmpolar = memory + nsites + nsites*3;
    atmpolar[0] = 1.310; // polarO
    atmpolar[1] = 0.294; // polarH
    atmpolar[2] = 0.294; // polarH
    atmpolar[3] = 0.0;   // polarM
    return atmpolar;
}

double* mon""" + str(mon_index) + """::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;
    polfac[0] = 1.310; // polarO
    polfac[1] = 0.294; // polarH
    polfac[2] = 0.294; // polarH
    polfac[3] = 1.31;   // polarM
    return polfac;
}

void mon""" + str(mon_index) + """::allocate() {
    memory = new double [nsites // charges
  + nsites*3              // sitecrds
  + nsites                // polarizabilities
  + nsites];              // polfacs
}

int mon""" + str(mon_index) + """::get_nsites() { return nsites; }
int mon""" + str(mon_index) + """::get_realsites() { return realsites; }
double* mon""" + str(mon_index) + """::get_sitecrds() { return sitecrds; }
double* mon""" + str(mon_index) + """::get_charges() { return charge; }
double* mon""" + str(mon_index) + """::get_polfacs() { return polfac; }
double* mon""" + str(mon_index) + """::get_pol() { return atmpolar; }

excluded_set_type::iterator mon""" + str(mon_index) + """::get_begin_12() { return excluded12.begin(); }
excluded_set_type::iterator mon""" + str(mon_index) + """::get_begin_13() { return excluded13.begin(); }
excluded_set_type::iterator mon""" + str(mon_index) + """::get_begin_14() { return excluded14.begin(); }
excluded_set_type::iterator mon""" + str(mon_index) + """::get_end_12() { return excluded12.end(); }
excluded_set_type::iterator mon""" + str(mon_index) + """::get_end_13() { return excluded13.end(); }
excluded_set_type::iterator mon""" + str(mon_index) + """::get_end_14() { return excluded14.end(); }


} // namespace x

"""
    with open(filename, 'w') as cpp_file:
        cpp_file.write(cpp_text)


def write_fit_polynomial_holder_header(system_name, number_of_monomers, non_linear_parameters, ri, ro):

    """
    Writes the polynomial holder for the fitting

    Args:
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        number_of_monomers     - Number of monomers in the system
        non_linear_parameters  - All the non linear parameters
        ri                     - Inner cutoff (not used if 1b)
        ro                     - Outer cutoff (not used in 2b)
    """

    defflag = "MBNRG_" + str(number_of_monomers) + "B_" + system_name + "_FIT"

    system_keyword_mbnrg = "mbnrg_" + str(number_of_monomers) + "b_" + system_name

    system_keyword_poly = "poly_" + str(number_of_monomers) + "b_" + system_name

    system_keyword_polyholder = "polyholder_" + str(number_of_monomers) + "b_" + system_name

    defflag = defflag.replace("(", "_o_").replace(")", "_c_")
    system_keyword_mbnrg = system_keyword_mbnrg.replace("(", "_o_").replace(")", "_c_")
    system_keyword_poly = system_keyword_poly.replace("(", "_o_").replace(")", "_c_")
    system_keyword_polyholder = system_keyword_polyholder.replace("(", "_o_").replace(")", "_c_")

    header_mbnrg_fit = system_keyword_mbnrg + "_fit.h"
    header_poly_fit = system_keyword_poly + "_fit.h"

    # Generate arguments for the eval
    args_eval = "const std::vector<double> mon1"
    for i in range(2, number_of_monomers + 1):
        args_eval += ", const std::vector<double> mon" + str(i)

    header_text = """#ifndef """ + defflag + """
#define """ + defflag + '''

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include <netcdf.h>


#include "constants.h"
#include "variable.h"
#include "water_monomer_lp.h"
#include "''' + header_poly_fit  + '''"

////////////////////////////////////////////////////////////////////////////////

namespace ''' + system_keyword_mbnrg + """_fit {

struct """ + system_keyword_polyholder + """_fit {

    typedef mb_system_fit::poly_model poly_type_fit;

    void cart_to_vars(const std::vector<double> xyz, std::vector<double> &vars, double& s, double& gs);

    double calculate(const std::vector<double> xyz);
    size_t get_nvars() {return poly_type_fit::n_vars;}
    size_t get_num_nonlinear_params() {return """ + str(len(non_linear_parameters)) + """;}
    void set_nonlinear_parameters(const double*);
    void set_linear_parameters(const double*);
    void get_nonlinear_parameters(double*);
    bool nonlinear_parameters_out_of_range();
    void load_netcdf(const char* fn);

    inline double calculate_nl_part_of_terms(const std::vector<double> xyz, std::vector<double> &nl_part_of_terms);
    inline size_t nparams() { return poly_type_fit::size;}

    inline void as_cdl(std::ostream&);
    void write_cdl(std::ostream&, unsigned, const double*);

  private:

"""
    for nl_param in non_linear_parameters:
        header_text += "    double m_" + nl_param + ";\n"

    header_text += """protected:
    double m_ri = """ + str(ri) + """;
    double m_ro = """ + str(ro) + """;

    double f_switch(const double, double&);

private:
    double m_poly[poly_type_fit::size];

};

// For fitting
inline void """ + system_keyword_polyholder + """_fit::as_cdl(std::ostream& os) {
    write_cdl(os, poly_type_fit::size, m_poly);
}


inline double """ + system_keyword_polyholder + """_fit::calculate_nl_part_of_terms(std::vector<double> xyz, std::vector<double> &nl_part_of_terms) {

    double s, gs;
    std::vector<double> v(poly_type_fit::n_vars,0.0);

    cart_to_vars(xyz, v, s, gs);

    poly_type_fit polyn;
    polyn.eval(v.data(), nl_part_of_terms.data());
    for (unsigned n = 0; n < poly_type_fit::size; ++n)
        nl_part_of_terms[n] *= s;

    return 0;
}

inline double """ + system_keyword_polyholder + """_fit::calculate(std::vector<double> xyz) {

    double s, gs;
    std::vector<double> v(poly_type_fit::n_vars,0.0);
    std::vector<double> nl_part_of_terms(poly_type_fit::size,0.0);

    cart_to_vars(xyz, v, s, gs);

    poly_type_fit polyn;
    polyn.eval(v.data(), nl_part_of_terms.data());
    double energy = 0.0;
    for (unsigned n = 0; n < poly_type_fit::size; ++n) {
        nl_part_of_terms[n] *= s*m_poly[n];
        energy += nl_part_of_terms[n];
    }

    return energy;
}

//----------------------------------------------------------------------------//

} // namespace """ + system_keyword_mbnrg + """_fit

////////////////////////////////////////////////////////////////////////////////

#endif

"""

    with open(header_mbnrg_fit, "w") as header_file:
        header_file.write(header_text)


def get_individual_atoms_with_type(monomer_atom_types, vsites):
    """
    From the monomer atom types ([A,1,B,1],[C,2],...) generates a list of lists that contain information about atom type, index, and monomer identity.

    Args:
        monomer_atom_types    - List with the atom symmetry and the number of atoms of that symmetry, e.g. [[A,1,B,1],[C,2]] for A1B2_C2
        vsites                - Virtual site labels
    """
    mon_id = 'a'
    atoms = []
    for monomer in monomer_atom_types:
        num_ats = int(len(monomer)/2 + 0.49999)
        atoms.append([])
        for i in range(num_ats):
            atom = monomer[2*i]
            if not atom in vsites:
                for j in range(monomer[2*i + 1]):
                    atoms[-1].append([atom,j+1,mon_id])
        mon_id = chr(ord(mon_id)+1)
    return atoms

def get_coords_var_name(symmetry_class, atom_index, fragment_index, extension=""):
    return "_".join(["coords", symmetry_class, str(atom_index), fragment_index, extension]
                    if extension != ""
                    else ["coords", symmetry_class, str(atom_index), fragment_index])


def get_C6_var_name(symmetry_class1, symmetry_class2, extension=""):
    return "_".join(["m", "C6", symmetry_class1, symmetry_class2, extension]
                    if extension != ""
                    else ["m", "C6", symmetry_class1, symmetry_class2])


def get_d6_var_name(symmetry_class1, symmetry_class2, extension=""):
    return "_".join(["m", "d6", symmetry_class1, symmetry_class2, extension]
                    if extension != ""
                    else ["m", "d6", symmetry_class1, symmetry_class2])


def get_A_var_name(symmetry_class1, symmetry_class2, extension=""):
    return "_".join(["m", "A", symmetry_class1, symmetry_class2, extension]
                    if extension != ""
                    else ["m", "A", symmetry_class1, symmetry_class2])


def get_b_var_name(symmetry_class1, symmetry_class2, extension=""):
    return "_".join(["m", "b", symmetry_class1, symmetry_class2, extension]
                    if extension != ""
                    else ["m", "b", symmetry_class1, symmetry_class2])

def get_pointer_setup_string(symmetry_parser, vsites, xyz_var_name, extension = "", nspaces = 4, is_const = True):
    """
    Returns a string with C++ code that initializes all the pointers to the coordinates or gradients.

    Args:
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        vsites                 - Virtual site labels
        xyz_var_name           - Name of the variable in C++ that contains the data for coordinates, gradients or any 3N array
        extension              - Extra characters that need to be added at the end of the variable name (A_1_a --> A_1_a_g)
                Default: ""
        nspaces                - Number of spaces in front of the variable declaration; 4 spaces = 1 indentation level
        is_const               - Specifies if the pointer is a constant or not
    """

    crd_shift = 0
    spaces = " "*nspaces
    mon_id = 'a'
    string_pointers = ""
    string_pointers_vs = ""
    prefix = ""

    if is_const:
        prefix = "const "

    for symmetry_class, atom_index, fragment_index in symmetry_parser.get_atoms():
        if not symmetry_class in vsites:
            string_pointers += "{}{}double* {} = {} + {};\n".format(spaces, prefix, get_coords_var_name(symmetry_class, atom_index, fragment_index, extension), xyz_var_name, crd_shift)
            string_pointers += "\n"
            crd_shift += 3
        else:
            string_pointers_vs += "{}{}double {}[3];\n".format(spaces, "", get_coords_var_name(symmetry_class, atom_index, fragment_index, extension))
            string_pointers_vs += "{}{}std::fill({},{} + 3, 0.0);\n".format(spaces,"",get_coords_var_name(symmetry_class, atom_index, fragment_index, extension),get_coords_var_name(symmetry_class, atom_index, fragment_index, extension))
            string_pointers_vs += "\n"
        mon_id = chr(ord(mon_id)+1)

    return string_pointers, string_pointers_vs


def get_variables_string(variables, name_vout, name_vstruct, nspaces = 4):
    """
    Returns a string with C++ code that initializes all the variables for the polynomials

    Args:
        variables              - List of lists with the information in the poly.in file
        name_vout              - Name of the variable that will be at the left side of the =
        name_vstruct           - Name of the variable that is the array of "variables"
        monomer_atom_types     - List with the atom symmetry and the number of atoms of that symmetry, e.g. [[A,1,B,1],[C,2]] for A1B2_C2
        nspaces                - Number of spaces in front of the variable declaration; 4 spaces = 1 indentation level
    """

    variables_string = ""
    spaces = " "*nspaces
    for variable_index, (atom1_sym, atom1_ind, atom1_frag, atom2_sym, atom2_ind, atom2_frag, category, type) in enumerate(variables):
        types_1 = utils_nb_fitting.get_atom_types(atom1_sym + str(atom1_ind))
        types_2 = utils_nb_fitting.get_atom_types(atom2_sym + str(atom2_ind))

        nl_param_k = "m_k_" + category.replace("-", "_").replace("+", "_")
        nl_param_d = "m_d_" + category.replace("-", "_").replace("+", "_")
        atom1_coords = get_coords_var_name(atom1_sym, atom1_ind, atom1_frag)
        atom2_coords = get_coords_var_name(atom2_sym, atom2_ind, atom2_frag)

        if "0" in type:
            variables_string += spaces + "{}[{}] = {}[{}].v_{}({}, {}, {}, {});\n".format(name_vout, variable_index, name_vstruct, variable_index, type, nl_param_d, nl_param_k, atom1_coords, atom2_coords)
        else:
            variables_string += spaces + "{}[{}] = {}[{}].v_{}({}, {}, {});\n".format(name_vout, variable_index, name_vstruct, variable_index, type, nl_param_k, atom1_coords, atom2_coords)

    return variables_string


def get_grad_var_string(variables, name_vin, name_g, nspaces = 4):
    """
    Returns a string with C++ code that initializes all the gradient redistribution for the variables of the polynomials

    Args:
        variables              - List of lists with the information in the poly.in file
        name_vin               - Name of the variable that contains the variable struct
        name_g                 - Name of the array that has the gradients of the variables
        nspaces                - Number of spaces in front of the variable declaration; 4 spaces = 1 indentation level
    """
    variables_string = ""
    spaces = " "*nspaces
    for i in range(len(variables)):
        var = variables[i]
        var1 = var[0] + str(var[1])
        var2 = var[3] + str(var[4])
        types_1 = utils_nb_fitting.get_atom_types(var1)
        types_2 = utils_nb_fitting.get_atom_types(var2)
        atom1_name = "{}_{}_{}".format(types_1[0], types_1[1], var[2])
        atom2_name = "{}_{}_{}".format(types_2[0], types_2[1], var[5])

        variables_string += spaces + "{0}[{1}].grads({2}[{1}], {3}, {4}, {5}, {6});\n".format(
                name_vin,
                i,
                name_g,
                get_coords_var_name(var[0], var[1], var[2], extension="g"),
                get_coords_var_name(var[3], var[4], var[5], extension="g"),
                get_coords_var_name(var[0], var[1], var[2]),
                get_coords_var_name(var[3], var[4], var[5]))

    return variables_string


def write_fit_polynomial_holder_cpp(system_name, symmetry_parser, number_of_monomers, number_of_atoms, vsites, use_lonepairs, non_linear_parameters, variables, number_of_variables, ri, ro, k_min_intra, k_max_intra, k_min, k_max, d_min_intra, d_max_intra, d_min, d_max):
    """
    Writes the polynomial holder C++ file

    Args:
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        number_of_monomers     - Number of monomers in the system
        number_of_atoms        - Number of real atoms in the monomer
        vsites                 - Virtual site labels
        use_lonepairs          - List with 0 and 1 of the same size as monomers. If 0, no lone pairs will be used or declared. If 1, lone pairs will be used. NOTE: TODO As for 09/25/2019, only water can have lone pairs.
        non_linear_parameters  - All the non linear parameters
        variables              - List of lists with the information in the poly.in file
        number_of_variables    - Number of variables
        ri                     - Inner cutoff (not used if 1b)
        ro                     - Outer cutoff (not used in 2b)
        k_min_intra            - Minimum value allowed for k_intra
        k_max_intra            - Maximum value allowed for k_intra
        k_min                  - Minimum value allowed for k
        k_max                  - Maximum value allowed for k
        d_min_intra            - Minimum value allowed for d_intra
        d_max_intra            - Maximum value allowed for d_intra
        d_min                  - Minimum value allowed for d
        d_max                  - Maximum value allowed for d
    """

    system_keyword_mbnrg = "mbnrg_" + str(number_of_monomers) + "b_" + system_name
    system_keyword_poly = "poly_" + str(number_of_monomers) + "b_" + system_name
    system_keyword_polyholder = "polyholder_" + str(number_of_monomers) + "b_" + system_name

    system_keyword_mbnrg = system_keyword_mbnrg.replace("(", "_o_").replace(")", "_c_")
    system_keyword_poly = system_keyword_poly.replace("(", "_o_").replace(")", "_c_")
    system_keyword_polyholder = system_keyword_polyholder.replace("(", "_o_").replace(")", "_c_")

    header_mbnrg_fit = system_keyword_mbnrg + "_fit.h"
    header_poly_fit = system_keyword_poly + "_fit.h"
    header_poly = system_keyword_poly + ".h"

    cpp_mbnrg_fit = system_keyword_mbnrg + "_fit.cpp"

    ff = open(cpp_mbnrg_fit,'w')

    a = '''#include "''' + header_mbnrg_fit + '''"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int kode) {

    std::cerr << " ** Fatal Error in ''' + system_keyword_polyholder + '''_fit::load_netcdf() ** :"
              << nc_strerror(kode) << std::endl;
    std::exit(EXIT_FAILURE);
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace ''' + system_keyword_mbnrg + """_fit {

void """ + system_keyword_polyholder + """_fit::set_linear_parameters(const double* xxx) {
  std::copy(xxx, xxx + nparams(), m_poly);
}


"""
    ff.write(a)

    ff.write('void ' + system_keyword_polyholder + '_fit::set_nonlinear_parameters(const double* xxx) { \n ')
    for nl in non_linear_parameters:
        ff.write("    m_" + nl + " = *xxx++; \n")
    a = """
}

//----------------------------------------------------------------------------//

"""
    ff.write(a)
    ff.write('void ' + system_keyword_polyholder + '_fit::get_nonlinear_parameters(double* xxx) \n { \n ')
    for nl in non_linear_parameters:
        ff.write("   *xxx++ = m_" + nl + "; \n")
    a = """
}

void """ + system_keyword_polyholder + """_fit::write_cdl(std::ostream& os, unsigned npoly, const double* poly) {

    using namespace std;

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific
       << "netcdf """ + system_keyword_mbnrg + """_fit {" << endl
       << "  // global attributes " << endl
       << "  :name = \\" """ + system_keyword_mbnrg + """_fit\\";" << endl;

    // """ + system_keyword_mbnrg + """_fit::as_cdl(os);
"""
    ff.write(a)
    ff.write('    os ')
    for nl in non_linear_parameters:
        ff.write('       << "  :' + nl + ' = " << setw(22) << m_' + nl +  ' << "; // A^(-1))" << endl \n')
    a = """
         << "  :ri = " << setw(22) << m_ri << "; // A" << endl
         << "  :ro = " << setw(22) << m_ro << "; // A" << endl
         << "  dimensions:" << endl
         << "  poly = " << npoly << ";" << endl
         << "  variables:" << endl
         << "    double poly(poly);" << endl
         << "data:" << endl;

    os << "poly =" << endl;

    {
        unsigned n(0);
        for (n = 0; n + 1 < npoly; ++n)
            os << std::setw(22) << poly[n] << ", // " << n << endl;
        os <<  std::setw(22) << poly[n] << "; // " << n << endl << endl;
    }

    os << "}" << endl;
    os.flags(saved_flags);
}

//----------------------------------------------------------------------------//
"""
    ff.write(a)
    ff.write('bool ' + system_keyword_polyholder + '_fit::nonlinear_parameters_out_of_range() { \n')
    ff.write('    const double k_min =  ' + k_min + ' ;\n')
    ff.write('    const double k_max =  ' + k_max + ' ;\n')
    ff.write('    const double k_min_intra =  ' + k_min_intra + ' ;\n')
    ff.write('    const double k_max_intra =  ' + k_max_intra + ' ;\n\n')

    ff.write('    const double d_min =  ' + d_min + ' ;\n')
    ff.write('    const double d_max =  ' + d_max + ' ;\n')
    ff.write('    const double d_min_intra =  ' + d_min_intra + ' ;\n')
    ff.write('    const double d_max_intra =  ' + d_max_intra + ' ;\n\n')


    ff.write('return false')
    for nl in non_linear_parameters:
        if nl.startswith('d'):
            if (nl.startswith('d_intra')):
                ff.write('\n       || m_' + nl + ' < d_min_intra ')
                ff.write('\n       || m_' + nl + ' > d_max_intra ')
            else:
                ff.write('\n       || m_' + nl + ' < d_min ')
                ff.write('\n       || m_' + nl + ' > d_max ')
        else:
            if (nl.startswith('k_intra')):
                ff.write('\n       || m_' + nl + ' < k_min_intra ')
                ff.write('\n       || m_' + nl + ' > k_max_intra ')
            else:
                ff.write('\n       || m_' + nl + ' < k_min ')
                ff.write('\n       || m_' + nl + ' > k_max ')

    ff.write('; \n')
    a = """
}

"""
    ff.write(a)

    all_distances = utils_nb_fitting.get_list_of_numeric_pairs("d",number_of_monomers)
    all_switches = utils_nb_fitting.get_list_of_numeric_pairs("sw",number_of_monomers)

    ff.write('void ' + system_keyword_polyholder + '_fit::cart_to_vars(const std::vector<double> xyz, std::vector<double> &vars, double& s, double& gs) { \n\n')

    # Write the distances
    for d in all_distances:
        shift1 = 0
        shift2 = 0
        for i in range(d[1] - 1):
            shift1 += number_of_atoms[i]
        for i in range(d[2] -1):
            shift2 += number_of_atoms[i]
        shift1 *= 3
        shift2 *= 3
        a = """
    const double """ + d[0] + """[3] =
                         {xyz[""" + str(shift1 + 0) + """] - xyz[""" + str(shift2 + 0) + """],
                          xyz[""" + str(shift1 + 1) + """] - xyz[""" + str(shift2 + 1) + """],
                          xyz[""" + str(shift1 + 2) + """] - xyz[""" + str(shift2 + 2) + """]};

    const double """ + d[0] + """rsq = """ + d[0] + """[0]*""" + d[0] + """[0] + """ + d[0] + """[1]*""" + d[0] + """[1] + """ + d[0] + """[2]*""" + d[0] + """[2];
    const double """ + d[0] + """r = std::sqrt(""" + d[0] + """rsq);
"""
        ff.write(a)

    # Write Check for distances but skip if distances is empty
    if len(all_distances) > 0:
        condition = "true "
        for d in all_distances:
            condition += " && " + d[0] + "r > m_ro "
    else:
        condition = "false "

    a = """
    vars = std::vector<double>(""" + str(number_of_variables) + """,0.0);

    if (""" + condition + """) {
         s=0.0;
         gs = 0.0;
         return;
    }

    double sw = 0.0;


"""

    ff.write(a)

    # Get the pointers to the atoms
    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(symmetry_parser, vsites, "xyz.data()")

    ff.write(pointer_to_coordinates)
    ff.write(pointer_to_vsites)

    # Write water monomer setup if needed
    # w12, w13 qne wcross will be unused if none of the monomers is water. Is OK.

    a = """
    double w12 =     -9.721486914088159e-02;  //from MBpol
    double w13 =     -9.721486914088159e-02;
    double wcross =   9.859272078406150e-02;
"""
    ff.write(a)

    fragments = symmetry_parser.get_sub_parsers()

    # FIXME Only monomer that accepts lone pairs, for now, is MBpol water.
    char_code = 'a'
    for i in range(len(use_lonepairs)):
        if use_lonepairs[i] != 0:
            atoms = list(fragments[i].get_atoms())
            a = "    monomer m{};\n".format(str(i + 1))
            a += "    m{}.setup({}, w12, wcross, {}, {});\n".format(str(i + 1), get_coords_var_name(atoms[0][0], 1, char_code),
                                                                    get_coords_var_name(atoms[3][0], 1, char_code),
                                                                    get_coords_var_name(atoms[4][0], 2, char_code))
            ff.write(a)
        char_code = chr(ord(char_code)+1)

    ff.write("\n    variable vs[" + str(number_of_variables) + "];\n\n")
    string_vars = get_variables_string(variables, "vars", "vs")
    ff.write(string_vars)

    # write the switches
    # for now, will be set as the sum of the products of each n-1 combination
    # for a 3b, it will be the sum of the product of each pair
    sw_string = "\n"
    for sw, d in zip(all_switches, all_distances):
        sw_string += "    double g{} = 0.0;\n".format(sw[0])
        sw_string += "    double {} = f_switch({}, g{});\n".format(sw[0],d[0] + "r" ,sw[0])

    ff.write(sw_string)

    all_switches_only = []
    for sw in all_switches:
        all_switches_only.append(sw[0])
    sw_comb = list(it.combinations(all_switches_only,number_of_monomers - 1))
    sw_string = " + "
    terms = []
    if len(sw_comb[0]) == 0:
        sw_string = "1.0"
    else:
        for swc in sw_comb:
            term = "*"
            term = term.join(swc)
            terms.append(term)

        sw_string = sw_string.join(terms)

    a = """
    sw = """ + sw_string + """;

    s = sw;
    gs = 0.0;
}

void """ + system_keyword_polyholder + """_fit::load_netcdf(const char* fn)
{
    assert(fn);

    int rc, ncid;
    if ((rc = nc_open(fn, NC_NOWRITE, &ncid)))
        error(rc);

#   define RETRIEVE(name) \\
    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \\
        error(rc);
"""
    ff.write(a)
    for nl in non_linear_parameters:
        ff.write("    RETRIEVE(" + nl + ")\n")
    a = """

    RETRIEVE(ri)
    RETRIEVE(ro)

#   undef RETRIEVE

    int varid;

    if ((rc = nc_inq_varid(ncid, "poly", &varid)))
        error(rc);

    for (size_t n = 0; n < poly_type_fit::size; ++n) {
        if ((rc = nc_get_var1_double(ncid, varid, &n, m_poly + n)))
            error(rc);
    }

    if ((rc = nc_close(ncid)))
        error(rc);
}

//----------------------------------------------------------------------------//

double """ + system_keyword_polyholder + """_fit::f_switch(const double r, double& g)
{
    if (r > m_ro) {
        g = 0.0;
        return 0.0;
    } else if (r > m_ri) {
        const double t1 = M_PI/(m_ro - m_ri);
        const double x = (r - m_ri)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}

} // namespace """ + system_keyword_mbnrg + """_fit

"""
    ff.write(a)

    ff.close()


def write_buckingham_header(symmetry_parser, virtual_sites_poly, A_buck, b_buck):
    """
    Writes the buckingham header

    Args:
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        virtual_sites_poly     - Virtual site labels
        A_buck                 - List of all the A coefficients for the system. If 1b, they should correspond to the A of the homodimer. If 2b, the A of the dimer itself.
        b_buck                 - List of all the b (d6) coefficients for the system. If 1b, they should correspond to the b (d6) of the homodimer. If 2b, the b (d6) of the dimer itself.
    """
    hname = "buckingham.h"
    ff = open(hname,'w')
    a = """
#ifndef BUCKINGHAM_H
#define BUCKINGHAM_H

#include <cmath>
#include <algorithm>
#include <vector>

struct mbnrg_buck {
  mbnrg_buck();
  mbnrg_buck(std::vector<double> c);
  ~mbnrg_buck();

"""
    ff.write(a)

    # Need to select if we have 1b or 2b
    if symmetry_parser.get_num_fragments() == 1:
        pairs = symmetry_parser.get_pairs(vsites=virtual_sites_poly)
    else:
        pairs = symmetry_parser.get_intermolecular_pairs(vsites=virtual_sites_poly)

    # Write A and b declaration
    a = ""
    b = ""
    for pair_index, pair in enumerate(pairs):
        a += "  double " + get_A_var_name(pair[0], pair[1]) + " = " + str(A_buck[pair_index]) + ";\n"
        b += "  double " + get_b_var_name(pair[0], pair[1]) + " = " + str(b_buck[pair_index]) + ";\n"

    ff.write(a + "\n")
    ff.write(b)

    a = """

  std::vector<double> xyz;

  double get_buckingham();

  std::vector<double> get_nonlinear_terms();

  inline void set_nonlinear_parameters(std::vector<double> b) {
"""
    ff.write(a)

    a = ""
    for pair_index, pair in enumerate(pairs):
        a += "    " + get_b_var_name(pair[0], pair[1]) + " = b[" + str(pair_index) + "];\n"

    ff.write(a)

    a = """
  }

  inline void set_linear_parameters(std::vector<double> a) {
"""
    ff.write(a)

    a = ""
    for pair_index, pair in enumerate(pairs):
        a += "    " + get_A_var_name(pair[0], pair[1]) + " = a[" + str(pair_index) + "];\n"

    ff.write(a)

    a = """
  }

  inline double buck_energy(const double a, const double b,
                     const double* p1, const double* p2,
                     double* g1, double* g2)
  {

    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double fac = a*exp(-b*r);
    const double grd = b/r*fac;

    g1[0] -= dx*grd;
    g2[0] += dx*grd;

    g1[1] -= dy*grd;
    g2[1] += dy*grd;

    g1[2] -= dz*grd;
    g2[2] += dz*grd;

    return fac;
  }

  inline double buck_energy(const double a, const double b,
                     const double* p1, const double* p2)
  {

    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double fac = a*exp(-b*r);

    return fac;
  }
};

#endif

"""
    ff.write(a)
    ff.close()


def write_buckingham_cpp(symmetry_parser, virtual_sites_poly, excl12 = None, excl13 = None, excl14 = None):
    """
    Writes the buckingham C++ file for the fitting

    Args:
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        virtual_sites_poly     - Virtual site labels
        excl12                 - List with the pairs that are excluded at 1-2 distance
        excl13                 - List with the pairs that are excluded at 1-3 distance
        excl14                 - List with the pairs that are excluded at 1-4 distance
    """
    # We will need the pairs:
    # Need to select if we have 1b or 2b
    if symmetry_parser.get_num_fragments() == 1:
        pairs = symmetry_parser.get_pairs(vsites=virtual_sites_poly)
    else:
        pairs = symmetry_parser.get_intermolecular_pairs(vsites=virtual_sites_poly)

    cppname = "buckingham.cpp"
    ff = open(cppname,'w')
    a = """#include "buckingham.h"

mbnrg_buck::mbnrg_buck() {}
mbnrg_buck::~mbnrg_buck() {}

mbnrg_buck::mbnrg_buck(std::vector<double> c) {
  xyz = c;
}

std::vector<double> mbnrg_buck::get_nonlinear_terms() {
    std::vector<double> nl_terms(""" + str(len(pairs)) + """,0.0);
    """
    ff.write(a)

    # Get pointer to coordinates only for real sites
    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(symmetry_parser, virtual_sites_poly, "xyz.data()")

    ff.write(pointer_to_coordinates)

    a = ""

    excl = []
    for excl_level in [excl12, excl13, excl14]:
        for excl_pair in excl_level:
            excl.append(excl_pair)

    # case where we have a monomer
    if symmetry_parser.get_num_fragments() == 1:
        mon = list(symmetry_parser.get_atoms())
        for atom1_index, (atom1_symmetry, atom1_sym_index, atom1_frag_index) in enumerate(mon):
            for atom2_index, (atom2_symmetry, atom2_sym_index, atom2_frag_index) in enumerate(mon[atom1_index+1:]):
                index_pair = [atom1_index, atom2_index + 1 + atom1_index]
                if index_pair not in excl:
                    pair = tuple(sorted([atom1_symmetry, atom2_symmetry]))
                    vector_index = str(pairs.index(pair))

                    a += "    nl_terms[" + vector_index + "] += buck_energy(1.0, {}, {}, {});\n".format(get_b_var_name(pair[0], pair[1]), get_coords_var_name(atom1_symmetry, atom1_sym_index, atom1_frag_index), get_coords_var_name(atom2_symmetry, atom2_sym_index, atom2_frag_index))

    # case in which we have a dimer
    elif symmetry_parser.get_num_fragments() == 2:
        mons = list(symmetry_parser.get_atoms())
        for atom1_index, (atom1_symmetry, atom1_sym_index, atom1_frag_index) in enumerate(mons):
            for atom2_index, (atom2_symmetry, atom2_sym_index, atom2_frag_index) in enumerate(mons[atom1_index+1:]):
                if atom1_frag_index[0] == atom2_frag_index[0]:
                    continue
                if atom1_symmetry in virtual_sites_poly or atom2_symmetry in virtual_sites_poly:
                    continue
                pair = tuple(sorted([atom1_symmetry, atom2_symmetry]))
                vector_index = str(pairs.index(pair))
                a += "    nl_terms[" + vector_index + "] += buck_energy(1.0, {}, {}, {});\n".format(get_b_var_name(pair[0], pair[1]), get_coords_var_name(atom1_symmetry, atom1_sym_index, atom1_frag_index), get_coords_var_name(atom2_symmetry, atom2_sym_index, atom2_frag_index))

    ff.write(a)

    a = """

    return nl_terms;
}

double mbnrg_buck::get_buckingham() {

    double buck = 0.0;

"""
    ff.write(a)

    # Get pointer to coordinates only for real sites
    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(symmetry_parser, virtual_sites_poly, "xyz.data()")

    ff.write(pointer_to_coordinates)

    a = ""

    excl = []
    for excl_level in [excl12, excl13, excl14]:
        for excl_pair in excl_level:
            excl.append(excl_pair)

    # case where we have a monomer
    if symmetry_parser.get_num_fragments() == 1:
        mon = list(symmetry_parser.get_atoms())
        for atom1_index, (atom1_symmetry, atom1_sym_index, atom1_frag_index) in enumerate(mon):
            for atom2_index, (atom2_symmetry, atom2_sym_index, atom2_frag_index) in enumerate(mon[atom1_index+1:]):
                index_pair = [atom1_index, atom2_index + 1 + atom1_index]
                if index_pair not in excl:
                    pair = tuple(sorted([atom1_symmetry, atom2_symmetry]))
                    vector_index = str(pairs.index(pair))

                    a += "    buck += buck_energy({}, {}, {}, {});\n".format(get_A_var_name(pair[0], pair[1]), get_b_var_name(pair[0], pair[1]), get_coords_var_name(atom1_symmetry, atom1_sym_index, atom1_frag_index), get_coords_var_name(atom2_symmetry, atom2_sym_index, atom2_frag_index))

    # case in which we have a dimer
    elif symmetry_parser.get_num_fragments() == 2:
        mons = list(symmetry_parser.get_atoms())
        for atom1_index, (atom1_symmetry, atom1_sym_index, atom1_frag_index) in enumerate(mons):
            for atom2_index, (atom2_symmetry, atom2_sym_index, atom2_frag_index) in enumerate(mons[atom1_index+1:]):
                if atom1_frag_index[0] == atom2_frag_index[0]:
                    continue
                if atom1_symmetry in virtual_sites_poly or atom2_symmetry in virtual_sites_poly:
                    continue
                pair = tuple(sorted([atom1_symmetry, atom2_symmetry]))
                a += "    buck += buck_energy({}, {}, {}, {});\n".format(get_A_var_name(pair[0], pair[1]), get_b_var_name(pair[0], pair[1]), get_coords_var_name(atom1_symmetry, atom1_sym_index, atom1_frag_index), get_coords_var_name(atom2_symmetry, atom2_sym_index, atom2_frag_index))

    ff.write(a)

    a = """
    return buck;
}

"""

    ff.write(a)
    ff.close()


def write_dispersion_header(symmetry_parser, virtual_sites_poly, c6, d6):
    """
    Writes the header file for the dispersion energy calculation in the fitting

    Args:
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        virtual_sites_poly    - Virtual site labels
        c6                    - List of all the C6 coefficients for the system. If 1b, they should correspond to the C6 of the homodimer. If 2b, the C6 of the dimer itself.
        d6                    - List of all the d6 coefficients for the system. If 1b, they should correspond to the d6 of the homodimer. If 2b, the d6 of the dimer itself.
    """
    hname = "dispersion.h"
    ff = open(hname,'w')
    a = """
#ifndef DISPERSION_H
#define DISPERSION_H

#include "tang-toennies.h"
#include <cmath>
#include <algorithm>
#include <vector>



struct mbnrg_disp {
  mbnrg_disp();
  mbnrg_disp(std::vector<double> c);
  ~mbnrg_disp();

"""
    ff.write(a)

    # Need to select if we have 1b or 2b
    if symmetry_parser.get_num_fragments() == 1:
        pairs = symmetry_parser.get_pairs(vsites=virtual_sites_poly)
    elif symmetry_parser.get_num_fragments() == 2:
        pairs = symmetry_parser.get_intermolecular_pairs(vsites=virtual_sites_poly)
    else:
        pairs = []

    # Write C6 declaration
    a = ""
    b = ""
    for pair_index, (symmetry1, symmetry2) in enumerate(pairs):
        a += "  double " + get_C6_var_name(symmetry1, symmetry2) + " = " + str(c6[pair_index]) + ";\n"
        b += "  double " + get_d6_var_name(symmetry1, symmetry2) + " = " + str(d6[pair_index]) + ";\n"

    ff.write(a + "\n")
    ff.write(b)

    a = """
  const double m_C8 = 0.0;
  const double m_d8 = 0.0;

  const double if6 = 1.0/x2o::factorial<6>();
  const double if8 = 1.0/x2o::factorial<8>();

  std::vector<double> xyz;

  double get_dispersion();

  std::vector<double> get_nonlinear_terms();

  inline void set_nonlinear_parameters(std::vector<double> d6) {
"""
    ff.write(a)

    a = ""
    for pair_index, (symmetry1, symmetry2) in enumerate(pairs):
        a += "    " + get_d6_var_name(symmetry1, symmetry2) + " = d6[" + str(pair_index) + "];\n"

    ff.write(a)

    a = """
  }

  inline void set_linear_parameters(std::vector<double> c6) {
"""
    ff.write(a)

    a = ""
    for pair_index, (symmetry1, symmetry2) in enumerate(pairs):
        a += "    " + get_C6_var_name(symmetry1, symmetry2) + " = c6[" + str(pair_index) + "];\n"

    ff.write(a)

    a = """
  }

  inline double x6(const double& C6, const double& d6,
                  const double& C8, const double& d8,
                  const double* p1, const double* p2,
                        double* g1,       double* g2)
  {
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r;
    const double tt6 = x2o::tang_toennies(6, d6r);

    const double d8r = d8*r;
    const double tt8 = x2o::tang_toennies(8, d8r);


    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;
    const double inv_r8 = inv_r6*inv_rsq;

    const double e6 = C6*tt6*inv_r6;
    const double e8 = C8*tt8*inv_r8;

    const double grd = (6*e6 + 8*e8)*inv_rsq
        - (C6*std::pow(d6, 7)*if6*std::exp(-d6r)
        +  C8*std::pow(d8, 9)*if8*std::exp(-d8r))/r;

    g1[0] += dx*grd;
    g2[0] -= dx*grd;

    g1[1] += dy*grd;
    g2[1] -= dy*grd;

    g1[2] += dz*grd;
    g2[2] -= dz*grd;

    return - (e6 + e8);
  }

  inline double x6(const double& C6, const double& d6,
                  const double& C8, const double& d8,
                  const double* p1, const double* p2)
  {
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r;
    const double tt6 = x2o::tang_toennies(6, d6r);

    const double d8r = d8*r;
    const double tt8 = x2o::tang_toennies(8, d8r);


    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;
    const double inv_r8 = inv_r6*inv_rsq;

    const double e6 = C6*tt6*inv_r6;
    const double e8 = C8*tt8*inv_r8;

    return - (e6 + e8);
  }

};

#endif

"""
    ff.write(a)
    ff.close()


def write_dispersion_cpp(symmetry_parser, virtual_sites_poly, excl12 = None, excl13 = None, excl14 = None):
    """
    Writes the C++ file for the dispersion energy calculation in the fitting

    Args:
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        virtual_sites_poly    - Virtual site labels
        excl12                - List with the pairs that are excluded at 1-2 distance
        excl13                - List with the pairs that are excluded at 1-3 distance
        excl14                - List with the pairs that are excluded at 1-4 distance
    """
    # Need to select if we have 1b or 2b
    if symmetry_parser.get_num_fragments() == 1:
        pairs = symmetry_parser.get_pairs(vsites=virtual_sites_poly)
    elif symmetry_parser.get_num_fragments() == 2:
        pairs = symmetry_parser.get_intermolecular_pairs(vsites=virtual_sites_poly)
    else:
        pairs = []

    cppname = "dispersion.cpp"
    ff = open(cppname,'w')
    a = """#include "dispersion.h"

mbnrg_disp::mbnrg_disp() {}
mbnrg_disp::~mbnrg_disp() {}

mbnrg_disp::mbnrg_disp(std::vector<double> c) {
  xyz = c;
}

std::vector<double> mbnrg_disp::get_nonlinear_terms() {
    std::vector<double> nl_terms(""" + str(len(pairs)) + """,0.0);
    """
    ff.write(a)

    # Get pointer to coordinates only for real sites
    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(symmetry_parser, virtual_sites_poly, "xyz.data()")

    ff.write(pointer_to_coordinates)

    # Need to select if we have 1b or 2b
    if symmetry_parser.get_num_fragments() == 1:
        pairs = symmetry_parser.get_pairs(vsites=virtual_sites_poly)
    elif symmetry_parser.get_num_fragments() == 2:
        pairs = symmetry_parser.get_intermolecular_pairs(vsites=virtual_sites_poly)
    else:
        pairs = []

    a = ""

    excl = []
    for excl_level in [excl12, excl13, excl14]:
        for excl_pair in excl_level:
            excl.append(excl_pair)

    # case where we have a monomer
    if symmetry_parser.get_num_fragments() == 1:
        mon = list(symmetry_parser.get_atoms())
        for atom1_index, (atom1_symmetry, atom1_sym_index, atom1_frag_index) in enumerate(mon):
            for atom2_index, (atom2_symmetry, atom2_sym_index, atom2_frag_index) in enumerate(mon[atom1_index+1:]):
                index_pair = [atom1_index, atom2_index + 1 + atom1_index]
                if index_pair not in excl:
                    pair = tuple(sorted([atom1_symmetry, atom2_symmetry]))
                    vector_index = str(pairs.index(pair))

                    a += "    nl_terms[" + vector_index + "] += x6(1.0, {}, m_C8, m_d8, {}, {});\n".format(get_d6_var_name(pair[0], pair[1]), get_coords_var_name(atom1_symmetry, atom1_sym_index, atom1_frag_index), get_coords_var_name(atom2_symmetry, atom2_sym_index, atom2_frag_index))

    # case in which we have a dimer
    elif symmetry_parser.get_num_fragments() == 2:
        mons = list(symmetry_parser.get_atoms())
        for atom1_index, (atom1_symmetry, atom1_sym_index, atom1_frag_index) in enumerate(mons):
            for atom2_index, (atom2_symmetry, atom2_sym_index, atom2_frag_index) in enumerate(mons[atom1_index+1:]):
                if atom1_frag_index[0] == atom2_frag_index[0]:
                    continue
                if atom1_symmetry in virtual_sites_poly or atom2_symmetry in virtual_sites_poly:
                    continue
                pair = tuple(sorted([atom1_symmetry, atom2_symmetry]))
                vector_index = str(pairs.index(pair))
                a += "    nl_terms[" + vector_index + "] += x6(1.0, {}, m_C8, m_d8, {}, {});\n".format(get_d6_var_name(pair[0], pair[1]), get_coords_var_name(atom1_symmetry, atom1_sym_index, atom1_frag_index), get_coords_var_name(atom2_symmetry, atom2_sym_index, atom2_frag_index))

    ff.write(a)

    a = """

    return nl_terms;
}

double mbnrg_disp::get_dispersion() {

    double disp = 0.0;

"""
    ff.write(a)

    # Get pointer to coordinates only for real sites
    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(symmetry_parser, virtual_sites_poly, "xyz.data()")

    ff.write(pointer_to_coordinates)

    a = ""

    excl = []
    for excl_level in [excl12, excl13, excl14]:
        for excl_pair in excl_level:
            excl.append(excl_pair)

    # Case in which we have a monomer
    if symmetry_parser.get_num_fragments() == 1:
        mon = list(symmetry_parser.get_atoms())
        for atom1_index, (atom1_symmetry, atom1_sym_index, atom1_frag_index) in enumerate(mon):
            for atom2_index, (atom2_symmetry, atom2_sym_index, atom2_frag_index) in enumerate(mon[atom1_index+1:]):
                index_pair = [atom1_index, atom2_index + 1 + atom1_index]
                if index_pair not in excl:
                    pair = sorted([atom1_symmetry, atom2_symmetry])
                    a += "    disp += x6({}, {}, m_C8, m_d8, {}, {});\n".format(get_C6_var_name(pair[0], pair[1]), get_d6_var_name(pair[0], pair[1]), get_coords_var_name(atom1_symmetry, atom1_sym_index, atom1_frag_index), get_coords_var_name(atom2_symmetry, atom2_sym_index, atom2_frag_index))

    # Case in which we have a dimer
    elif symmetry_parser.get_num_fragments() == 2:
        mons = list(symmetry_parser.get_atoms())
        for atom1_index, (atom1_symmetry, atom1_sym_index, atom1_frag_index) in enumerate(mons):
            for atom2_index, (atom2_symmetry, atom2_sym_index, atom2_frag_index) in enumerate(mons[atom1_index+1:]):
                if atom1_frag_index[0] == atom2_frag_index[0]:
                    continue
                if atom1_symmetry in virtual_sites_poly or atom2_symmetry in virtual_sites_poly:
                    continue
                pair = tuple(sorted([atom1_symmetry, atom2_symmetry]))
                a += "    disp += x6({}, {}, m_C8, m_d8, {}, {});\n".format(get_C6_var_name(pair[0], pair[1]), get_d6_var_name(pair[0], pair[1]), get_coords_var_name(atom1_symmetry, atom1_sym_index, atom1_frag_index), get_coords_var_name(atom2_symmetry, atom2_sym_index, atom2_frag_index))

    ff.write(a)

    a = """
    return disp;
}

"""

    ff.write(a)
    ff.close()


def get_nl_params_initialization_string(nl_param_all, k_min, k_max, d_min, d_max, k_min_intra, k_max_intra, d_min_intra, d_max_intra):
    """
    Returns a string that initializes the non-linear parameters in the C+++ code within the proper range.

    Args:
        nl_param_all           - All the non linear parameters
        k_min_intra            - Minimum value allowed for k_intra
        k_max_intra            - Maximum value allowed for k_intra
        k_min                  - Minimum value allowed for k
        k_max                  - Maximum value allowed for k
        d_min_intra            - Minimum value allowed for d_intra
        d_max_intra            - Maximum value allowed for d_intra
        d_min                  - Minimum value allowed for d
        d_max                  - Maximum value allowed for d
    """

    nl_params_string = ""
    for i in range(len(nl_param_all)):
        if nl_param_all[i].startswith('d'):
            if nl_param_all[i].startswith('d_intra'):
                nl_params_string += '    x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(d_max_intra) - float(d_min_intra)) + ' + ' + d_min_intra + ';\n'
            else:
                nl_params_string += '      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(d_max) - float(d_min)) + ' + ' + d_min + ';\n'
        else:
            if nl_param_all[i].startswith('k_intra'):
                nl_params_string += '      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(k_max_intra) - float(k_min_intra)) + ' + ' + k_min_intra + ';\n'
            else:
                nl_params_string += '      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(k_max) - float(k_min)) + ' + ' + k_min + ';\n'

    return nl_params_string


def get_nbody_electrostatics_string(number_of_monomers, number_of_atoms, number_of_sites, ptr_to_coords):
    """
    Returns a string that contains the C++ code with the n-body electrostatics calculation for the n-body system that is being fit.

    Args:
        number_of_monomers     - Number of monomers in the system
        number_of_atoms        - Number of real atoms in the monomer
        number_of_sites        - Number of sites (electrostatic sites) of the monomer
        ptr_to_coords          - Name of the variable that points to the coordinates in the C++ code.
    """
    electrostatics_string = ""

    # generate first index list to locate monomers
    first_index_of_atoms = []
    first_index_of_sites = []
    fi = 0
    fi_s = 0
    for i in range(number_of_monomers):
        first_index_of_atoms.append(fi)
        first_index_of_sites.append(fi_s)
        fi += number_of_atoms[i]
        fi_s += number_of_sites[i]

    mon_indexes = list(range(number_of_monomers))

    # Declare all the monomers:
    for i in range(number_of_monomers):
        electrostatics_string += "        x::mon" + str(i+1) + " m" + str(i+1) + "(" + ptr_to_coords + " + " + str(3*first_index_of_atoms[i]) + ");\n"

    electrostatics_string += "\n"

    nb_elec_string = []
    factors = []
    # Loop over all n-mers: monomers, dimers, trimers
    for i in range(number_of_monomers):
        nmer = i+1
        factors.append((-1)**(number_of_monomers - i + 1))
        nb_elec_string.append([])

        elec_keyword = ""
        for comb in it.combinations(mon_indexes,nmer):
            # set up nmer properties:
            nsites = 0
            elec_keyword = "_".join([str(x) for x in comb])
            fi_of_nmer_sites = []
            for mon in comb:
                fi_of_nmer_sites.append(nsites)
                nsites += number_of_sites[mon]

            # declare the energy outside the scope:
            electrostatics_string += "        double elec_" + elec_keyword + " = 0.0;\n"
            # declare the variables inside a scope
            a = """
        {
            int num_sites = """ + str(nsites) + """;

            int isw_sites[num_sites];

            double xyz_sites[num_sites*3];
            double chg_sites[num_sites];
            double pol_sites[num_sites];
            double polfac_sites[num_sites];
            excluded_set_type excl12;
            excluded_set_type excl13;
            excluded_set_type excl14;
"""
            electrostatics_string += a

            # Now fill in the data
            for mon in comb:
                mon_index = comb.index(mon)
                a = """
            std::fill(isw_sites + {1}, isw_sites + {1} + m{0}.get_nsites(), m{0}.is_w);
            std::copy(m{0}.get_sitecrds(), m{0}.get_sitecrds() + 3*m{0}.get_nsites(), xyz_sites + 3*{1});
            std::copy(m{0}.get_charges(), m{0}.get_charges() + m{0}.get_nsites(), chg_sites + {1});
            std::copy(m{0}.get_pol(), m{0}.get_pol() + m{0}.get_nsites(), pol_sites + {1});
            std::copy(m{0}.get_polfacs(), m{0}.get_polfacs() + m{0}.get_nsites(), polfac_sites + {1});

""".format(mon+1,fi_of_nmer_sites[mon_index],elec_keyword)

                electrostatics_string += a
                a = """
            for (auto i = m{0}.get_begin_12(); i != m{0}.get_end_12(); i++) {{
                std::pair<size_t,size_t> p =
                    std::make_pair(i->first + {1},
                                   i->second + {1});
                excl12.insert(p);
            }}

            for (auto i = m{0}.get_begin_13(); i != m{0}.get_end_13(); i++) {{
                std::pair<size_t,size_t> p =
                    std::make_pair(i->first + {1},
                                   i->second + {1});
                excl13.insert(p);
            }}

            for (auto i = m{0}.get_begin_14(); i != m{0}.get_end_14(); i++) {{
                std::pair<size_t,size_t> p =
                    std::make_pair(i->first + {1},
                                   i->second + {1});
                excl14.insert(p);
            }}

""".format(mon+1,fi_of_nmer_sites[mon_index],elec_keyword)

                electrostatics_string += a
            nb_elec_string[-1].append("elec_" + elec_keyword)

            a = """

            ttm::electrostatics m_electrostatics;

            ttm::smear_ttm4x smr;
            smr.m_aDD_intra_12 = 0.3;
            smr.m_aDD_intra_13 = 0.3;
            smr.m_aDD_intra_14 = 0.055;

            elec_{2} = m_electrostatics(num_sites, chg_sites, polfac_sites, pol_sites,
                                          xyz_sites, excl12, excl13, excl14,
                                          isw_sites, smr, 0);
        }}
""".format(mon+1,fi_of_nmer_sites[mon_index],elec_keyword)

            electrostatics_string += a

    nb_elec = []
    for i in range(number_of_monomers):
         nb_elec.append(str(factors[i]) + "*(")
         nb_elec[-1] += " + ".join(nb_elec_string[i])
         nb_elec[-1] += ")"

    electrostatics_string += "\n        elec_e.push_back(" + " + ".join(nb_elec) + ");\n"

    return electrostatics_string


def write_ttmnrg_fitting_code(symmetry_parser, virtual_sites_poly, number_of_monomers, number_of_atoms, number_of_sites, system_name, b_min, b_max, b_min_init, b_max_init, E_range):

    """
    Writes the fitting code for TTM-nrg"

    Args:
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        virtual_sites_poly     - Virtual site labels
        number_of_monomers     - Number of monomers in the system
        number_of_atoms        - Number of real atoms in the monomer
        number_of_sites        - Number of sites (electrostatic sites) of the monomer
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        b_min                  - Minimum value allowed for b
        b_max                  - Maximum value allowed for b
        b_min_init             - Minimum value allowed for initializing b
        b_max_init             - Maximum value allowed for initializing b
        E_range                - Energy range used for computing weights for each point in the fittings
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    # We need the pairs again to know how many linear and non linear terms we have
    pairs = symmetry_parser.get_intermolecular_pairs(vsites=virtual_sites_poly)

    num_terms = len(pairs)

    cppname = "fit-" + str(number_of_monomers) + "b-ttm.cpp"
    ff = open(cppname,'w')
    a = """#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <vector>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <chrono>

#include <gsl/gsl_multimin.h>

#include "wlsq.h"

#include "fit-utils.h"
#include "training_set.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"
#include "buckingham.h"
#include "constants.h"
"""

    ff.write(a)

    for i in range(1,number_of_monomers+1):
        ff.write("#include \"mon" + str(i) + ".h\"\n")

    a = """
namespace {

double E_range = """ + str(E_range) + """; // kcal/mol

//----------------------------------------------------------------------------//

static std::vector<tset::nb_system> training_set;
static std::vector<double> elec_e;
static std::vector<double> rep_e;
static std::vector<double> disp_e;
static std::vector<double> ts_weights;

size_t num_linear_params;
size_t num_non_linear_params;

namespace linear {

//----------------------------------------------------------------------------//

static double* A(0);
static double* y(0);
static double* params(0);

//----------------------------------------------------------------------------//

void allocate()
{
    A = new double[training_set.size()*num_linear_params];
    y = new double[training_set.size()];

    rep_e = std::vector<double>(training_set.size());
    disp_e = std::vector<double>(training_set.size());

    params = new double[num_linear_params];
}

//----------------------------------------------------------------------------//

double compute_chisq(const gsl_vector* X, void* unused) {

    for (size_t i = 0; i < num_non_linear_params; i++) {
        if (X->data[i] < """ + b_min + """ || X->data[i] > """ + b_max + """) return 1e+06;
    }

    std::vector<double> nlp(X->data, X->data + num_non_linear_params);

    for (size_t n = 0; n < training_set.size(); ++n) {
        // Calculate dispersion
        mbnrg_disp disp(training_set[n].xyz);
        disp.set_nonlinear_parameters(nlp);
        disp_e[n] = disp.get_dispersion();

        // Calculate non-linear terms of buckingham
        mbnrg_buck buck(training_set[n].xyz);
        buck.set_nonlinear_parameters(nlp);
        std::vector<double> nl_terms = buck.get_nonlinear_terms();

        // Update energy
        y[n] = training_set[n].nb_energy - disp_e[n];

        // Update the A matrix
        for (size_t p = 0; p < num_linear_params; ++p) {
            A[num_linear_params*n + p] = nl_terms[p];
        }
    }

    double chisq = 0.0;
    int rank = 0;

    kit::wlsq::solve(training_set.size(), num_linear_params,
                      A, y, ts_weights.data(), params, chisq, rank);

    std::cout << "<#> chisq = " << chisq
              << std::endl;

    for (size_t i = 0; i < num_linear_params; i++)
        if (params[i] < 0)
            chisq += fabs(params[i] * 100000);

    if (!gsl_finite (chisq)) {
      return 10000000.0;
    }

    return chisq;
}

//----------------------------------------------------------------------------//

} // namespace linear

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: ./fit-""" + str(number_of_monomers) + """b-ttm <training_set.xyz> [DE] > fit.log"
                  << std::endl;
        return 0;
    }

    long long int duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch()).count();

    srand(duration);

    double x0[""" + str(num_terms) + """];
"""

    ff.write(a)

    for i in range(num_terms):
        ff.write("    x0[" + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(b_max_init) - float(b_min_init)) + ' + ' + b_min_init + ';\n')

    a = """

    for (size_t i = 0; i < """ + str(num_terms) + """ ; i++) {
        std::cout << std::setprecision(6) << "x0[" << i << "] = " 
                  << x0[i] << std::endl;
    }

    argv++;
    argc--;

    try {
        size_t nsys = tset::load_nb_system(*argv, training_set);
        std::cout << "'" << *(argv++) << "' : "
                      << nsys << " configurations" << std::endl;
        if (--argc > 0) E_range = atof(*(argv++));
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    // Allocating memory
    num_linear_params = """ + str(num_terms) + """;
    num_non_linear_params = """ + str(num_terms) + """;
    linear::allocate();

    ts_weights = std::vector<double>(training_set.size(),1.0);

    double E_min;
    size_t N_eff;
    size_t index_min;

    tset::setup_weights(training_set, E_range, E_min, index_min,
                       ts_weights, N_eff);

    std::cout << "\\n>>   E_min = " << E_min << " kcal/mol; index " << index_min <<
                 "\\n>> E_range = " << E_range << " kcal/mol" <<
                 "\\n\\n>> training set size = " << training_set.size()
              << "\\n>>    effective size = " << N_eff << '\\n'
              << std::endl;

    std::vector<double> ts_energies(training_set.size(),0.0);

        for (size_t n = 0; n < training_set.size(); n++) {
        // Saving reference energies
        ts_energies[n] = training_set[n].nb_energy ;
"""

    ff.write(a)

    # Now write the many body electrostatics calculation
    a = get_nbody_electrostatics_string(number_of_monomers, number_of_atoms, number_of_sites, "training_set[n].xyz.data()")

    ff.write(a)

    a = """
        training_set[n].nb_energy -= elec_e[n];
    }

    gsl_vector* x = gsl_vector_alloc(num_non_linear_params);
    std::copy(x0, x0 + num_non_linear_params, x->data);

    gsl_multimin_function chisq_func;

    chisq_func.n = num_non_linear_params;
    chisq_func.f = linear::compute_chisq;

    gsl_vector* ss = gsl_vector_alloc(num_non_linear_params);
    gsl_vector_set_all(ss, 0.1);

    std::cout << "\\n<> initial simplex sides:\\n";
    for (size_t n = 0; n < num_non_linear_params; ++n)
        std::cout << n << " : " << ss->data[n] << "\\n";

    gsl_multimin_fminimizer* s =
        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2rand,
                                      num_non_linear_params);

    gsl_multimin_fminimizer_set(s, &chisq_func, x, ss);

    int iter(0), status(0);
    //
    // Main optimization loop
    //

        do {
        ++iter;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
             break;

        const double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-4);
        //status = gsl_multimin_test_size(size, 1e-3); //changed it to test convergence

        if (status == GSL_SUCCESS)
            std::cout << "!!! converged !!!" << std::endl;

        std::cout << iter << ' ' << s->fval << ' ' << size << std::endl;
        if (iter > 0 && iter%10 == 0) {
            std::cout << "\\n<> solution:\\n";
            for (size_t n = 0; n <  num_non_linear_params; ++n)
                std::cout << n << " : " << s->x->data[n] << "\\n";
            std::cout << "<>" << std::endl;
        }
    } while (status == GSL_CONTINUE && iter < 5000);

    linear::compute_chisq(s->x, 0);
    std::vector<double> final_l(linear::params, linear::params + num_linear_params);
    std::vector<double> final_nl(s->x->data,s->x->data + num_non_linear_params);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    //
    // report
    //

    for (size_t i = 0; i < training_set.size(); ++i)
        training_set[i].nb_energy += elec_e[i];


    std::ofstream correlation_file;
    correlation_file.open ("correlation.dat");
    
    correlation_file << "#" << std::setw(10) << "index"
                     << std::setw(15)       << "Reference" 
                     << std::setw(15)       << "Calculated"
                     << std::setw(15)       << "Error^2"
                     << "    \\n"; 

    std::ofstream terms_file;
    terms_file.open("individual_terms.dat");
    terms_file << "#" << std::setw(10) << "index"
                     << std::setw(20)       << "Reference" 
                     << std::setw(20)       << "Calculated"
                     << std::setw(20)       << "Repulsion (TTM-nrg)"
                     << std::setw(20)       << "Dispersion"
                     << std::setw(40)       << "Electrostatics (Permanent + Induced)"
                     << "    \\n";
    
    double err_L2(0), err_wL2(0), err_Linf(0);
    double err_L2_lo(0), err_Linf_lo(0), nlo(0);
    double weight_sum(0);

    for (size_t i = 0; i < training_set.size(); ++i) {

        // Calculate dispersion
        mbnrg_disp disp(training_set[i].xyz);
        disp.set_nonlinear_parameters(final_nl);
        disp_e[i] = disp.get_dispersion();

        // Calculate non-linear terms of buckingham
        mbnrg_buck buck(training_set[i].xyz);
        buck.set_nonlinear_parameters(final_nl);
        buck.set_linear_parameters(final_l);
        rep_e[i] = buck.get_buckingham();

        const double E_model = rep_e[i] + elec_e[i] + disp_e[i];
        const double delta = E_model - training_set[i].nb_energy;
        if (std::abs(delta) > err_Linf)
            err_Linf = std::abs(delta);

        correlation_file <<  std::setw(10) << i+1
                         << std::setw(20) << std::scientific
                         << std::setprecision(8)  <<  training_set[i].nb_energy
                         << std::setw(20) <<  E_model
                         << std::setw(20) <<  delta*delta  << "    \\n" ;

        terms_file <<  std::setw(10) << i+1
                         << std::setw(20) << std::scientific
                         << std::setprecision(8)  <<  training_set[i].nb_energy
                         << std::setw(20) <<  E_model
                         << std::setw(20) <<  rep_e[i] 
                         << std::setw(20) <<  disp_e[i]
                         << std::setw(40) <<  elec_e[i]
                         << "    \\n" ;

        err_L2 += delta*delta;
        err_wL2 += ts_weights[i]*delta*delta;
        weight_sum += ts_weights[i];

        if (training_set[i].binding_energy - E_min < E_range) {
            nlo += 1.0;
            err_L2_lo += delta*delta;
            if (std::abs(delta) > err_Linf_lo)
                err_Linf_lo = std::abs(delta);
        }
    }

    correlation_file.close();
    terms_file.close();

    err_L2 /= training_set.size();
    err_wL2 /= weight_sum;

    err_L2_lo /= nlo;

    std::cout << "      err[L2] = " << std::sqrt(err_L2) << "    #rmsd of full ts\\n"
              << "     err[wL2] = " << std::sqrt(err_wL2) << "   #weighted rmsd of full ts\\n"
              << "    err[Linf] = " << err_Linf << "   #highest error in full ts\\n"
              << "  err[L2,low] = " << std::sqrt(err_L2_lo) << "   #rmsd of low-energy ts \\n"
              << "err[Linf,low] = " << err_Linf_lo << "   #highest error in low-energy ts "
              << std::endl;

    std::cout << "\\n>> saving as '" << "ttm-nrg_params.dat" << "'\\n";
    std::ofstream final_params;
    final_params.open ("ttm-nrg_params.dat");
    for (size_t i = 0; i < num_linear_params; i++) {
        final_params << final_l[i] << " ";
    }
    final_params << std::endl;

    for (size_t i = 0; i < num_non_linear_params; i++) {
        final_params << final_nl[i] << " ";
    }
    final_params << std::endl;

    final_params.close();

    return 0;
}

"""

    ff.write(a)
    ff.close()

def write_ttmnrg_eval_code(symmetry_parser, virtual_sites_poly, number_of_monomers, number_of_atoms, number_of_sites, system_name):
    """
    Writes the eval code for TTM-nrg"

    Args:
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        virtual_sites_poly     - Virtual site labels
        number_of_monomers     - Number of monomers in the system
        number_of_atoms        - Number of real atoms in the monomer
        number_of_sites        - Number of sites (electrostatic sites) of the monomer
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    # We need the pairs again to know how many linear and non linear terms we have
    pairs = symmetry_parser.get_intermolecular_pairs(vsites=virtual_sites_poly)

    cppname = "eval-" + str(number_of_monomers) + "b-ttm.cpp"
    ff = open(cppname,'w')

    a = """#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <iterator>

#include "fit-utils.h"
#include "training_set.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"
#include "buckingham.h"
#include "constants.h"
"""

    ff.write(a)

    for i in range(1,number_of_monomers+1):
        ff.write("#include \"mon" + str(i) + ".h\"\n")

    a = """

int main(int argc, char** argv) {

    std::vector<tset::nb_system> training_set;
    std::vector<double> elec_e;
    std::vector<double> disp_e;
    std::vector<double> rep_e;
    std::vector<double> nb_energy;
    std::vector<double> ts_weights;

    std::vector<double> l_params(""" + str(len(pairs)) + """,0.0);
    std::vector<double> nl_params(""" + str(len(pairs)) + """,0.0);

    if (argc < 3) {
        std::cerr << "usage: ./eval-""" + str(number_of_monomers) + """b-ttm  <parameters.dat> <test_set.xyz> "
                  << std::endl;
        return 0;
    }

    argv++;
    argc--;

    try {
        std::ifstream params;
        params.open(*argv);

        // Read the file ensuring that its format is as it is supposed to be
        size_t lineno = 1;
        while (!params.eof()) {
            std::string line;
            getline(params,line);

            // We reached the end of the file
            if (line == "") break;

            // Check that only 2 lines exist in the parameter file (A and b)
            if (lineno == 3) {
                throw std::string("The TTM-nrg parameter file has more than two lines. This does not seem correct...");
            }
            
            // Retrieve the elements of the line and check proper size
            std::istringstream iss(line);
            std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());
            if (nl_params.size() != results.size()) {
                std::ostringstream ss;
                ss << "There is an issue with your TTM-nrg parameter file. "
                   << "The size read for line " << lineno << " is " << results.size()
                   << " when according to our records should be " << nl_params.size();
                throw ss.str();
            }

            // Now convert the strings to float
            if (lineno == 1) {
                for (size_t i = 0; i < results.size(); i++) {
                    l_params[i] = std::stod(results[i]);
                }
            } else if (lineno == 2) {
                for (size_t i = 0; i < results.size(); i++) {
                    nl_params[i] = std::stod(results[i]);
                }
            } else {
                throw "If you see this message, the world is about to end.";
            }
            lineno++;
        }

        ++argv;
        --argc;

        size_t nsys = tset::load_nb_system(*argv, training_set, false);
        std::cout << "'" << *(argv++) << "' : "
                      << nsys << " configurations" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    } catch (std::string s) {
        std::cerr << s << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error while parsing inputs...";
        return 1;
    }

    for (size_t n = 0; n < training_set.size(); n++) {
"""

    ff.write(a)

    a = get_nbody_electrostatics_string(number_of_monomers, number_of_atoms, number_of_sites, "training_set[n].xyz.data()")

    ff.write(a)

    a = """
        mbnrg_disp disp(training_set[n].xyz);
        disp.set_nonlinear_parameters(nl_params);
        disp_e.push_back(disp.get_dispersion());

        mbnrg_buck buck(training_set[n].xyz);
        buck.set_nonlinear_parameters(nl_params);
        buck.set_linear_parameters(l_params);
        rep_e.push_back(buck.get_buckingham());

        nb_energy.push_back(elec_e[n] + disp_e[n] + rep_e[n]);
    }

    std::cerr << std::setw(9)   << "#frame["
              << std::setw(6)   << std::setfill('.') << "###"
              << std::setw(3)   << "]= " << std::setfill(' ')
              << std::setw(20)  << "Calculated"
              << std::setw(20)  << "Repulsion (TTM-nrg)"
              << std::setw(20)  << "Dispersion"
              << std::setw(40)  << "Electrostatics (Permanent + Induced)"
              << std::endl;

    for (size_t n = 0; n < training_set.size(); n++) {
        std::cout << std::scientific << std::setprecision(8)
                  << std::setw(10)   << "frame["
                  << std::setw(6)    << std::setfill('.') << n
                  << std::setw(3)    << "]= " << std::setfill(' ')
                  << std::setw(20)   << nb_energy[n]
                  << std::setw(20)   << rep_e[n]
                  << std::setw(20)   << disp_e[n]
                  << std::setw(40)   << elec_e[n]
                  << std::endl;
    }


    return 0;
}

"""
    ff.write(a)
    ff.close()


def write_mbnrg_fitting_code(number_of_monomers, number_of_atoms, number_of_sites, system_name, nl_param_all, k_min_init, k_max_init, d_min_init, d_max_init, k_min_intra_init, k_max_intra_init, d_min_intra_init, d_max_intra_init, E_range, alpha):
    """
    Writes the fitting code for MB-nrg"

    Args:
        number_of_monomers     - Number of monomers in the system
        number_of_atoms        - Number of real atoms in the monomer
        number_of_sites        - Number of sites (electrostatic sites) of the monomer
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        nl_param_all           - All the non linear parameters
        k_min_init             - Minimum value allowed for initializing k
        k_max_init             - Maximum value allowed for initializing k
        d_min_init             - Minimum value allowed for initializing d
        d_max_init             - Maximum value allowed for initializing d
        k_min_intra_init       - Minimum value allowed for initializing k_intra
        k_max_intra_init       - Maximum value allowed for initializing k_intra
        d_min_intra_init       - Minimum value allowed for initializing d_intra
        d_max_intra_init       - Maximum value allowed for initializing d_intra
        E_range                - Energy range used for computing weights for each point in the fittings
        alpha                  - Parameter used for regularization during the fittings
    """
    system_name = system_name.replace("(", "_o_").replace(")", "_c_")

    cppname = "fit-" + str(number_of_monomers) + "b.cpp"
    ff = open(cppname,'w')
    a = """#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <chrono>

#include <gsl/gsl_multimin.h>

#include "rwlsq.h"

#include "fit-utils.h"
#include "training_set.h"
#include "mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit.h"
#include "electrostatics.h"
#include "coulomb.h"
"""

    ff.write(a)

    a = """#include "dispersion.h"
#ifdef USE_BUCKINGHAM
#include "buckingham.h"
#endif
"""

    if (number_of_monomers<3): ff.write(a)

    a = """#include "constants.h"
"""

    ff.write(a)

    for i in range(1,number_of_monomers+1):
        ff.write("#include \"mon" + str(i) + ".h\"\n")

    a = """
namespace {

static mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit::polyholder_""" + str(number_of_monomers) + "b_" + system_name + """_fit model;

double alpha = """ + str(alpha) + """;

double E_range = """ + str(E_range) + """; // kcal/mol

//----------------------------------------------------------------------------//

static std::vector<tset::nb_system> training_set;
static std::vector<double> elec_e;

"""
    ff.write(a)

    a = """#ifdef USE_BUCKINGHAM
static std::vector<double> buck_e;
#endif
static std::vector<double> disp_e;
"""
    if (number_of_monomers<3): ff.write(a)

    a = """static std::vector<double> ts_weights;

namespace linear {

//----------------------------------------------------------------------------//

static double* A(0);
static double* y(0);
static double* params(0);

//----------------------------------------------------------------------------//

void allocate()
{
    A = new double[training_set.size()*model.nparams()];
    y = new double[training_set.size()];

    params = new double[model.nparams()];
}

//----------------------------------------------------------------------------//

double compute_chisq(const gsl_vector* X, void* unused)
{
    model.set_nonlinear_parameters(X->data);

    if (model.nonlinear_parameters_out_of_range())
        return 1.0e+6;

    for (size_t n = 0; n < training_set.size(); ++n) {
        std::vector<double> mmm(model.nparams(),0.0);
        y[n] = training_set[n].nb_energy
               - model.calculate_nl_part_of_terms(training_set[n].xyz, mmm);
        for (size_t p = 0; p < model.nparams(); ++p) {
            A[p + n*model.nparams()] = mmm[p];
        }
    }

    double chisq;

    double penaltysq;
    kit::rwlsq::solve(training_set.size(), model.nparams(),
                      A, y, ts_weights.data(), alpha, params, chisq, penaltysq);

    std::cout << "<#> chisq = " << chisq
              << " : penaltysq = " << penaltysq
              << std::endl;

    if (!gsl_finite (chisq)) {
      return 10000000.0;
    }

    return chisq;
}

//----------------------------------------------------------------------------//

} // namespace linear

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: ./fit-""" + str(number_of_monomers) + """b <training_set.xyz> [DE] [ridge_alpha] > fit.log"
                  << std::endl;
        return 0;
    }

    long long int duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch()).count();

    srand(duration);

    double x0[""" + str(len(nl_param_all)) + """];

"""

    ff.write(a)

    nl_params_string = get_nl_params_initialization_string(nl_param_all, k_min_init, k_max_init, d_min_init, d_max_init, k_min_intra_init, k_max_intra_init, d_min_intra_init, d_max_intra_init)

    ff.write(nl_params_string)

    a = """
    model.set_nonlinear_parameters(x0);

    std::cout << "<> using ridge regression with alpha = "
              << alpha << std::endl;

    {
        const char fn[] = "fit-""" + str(number_of_monomers) + """b-initial.cdl";
        std::ofstream ofs(fn);
        std::cout << "\\n>> dumping initial model as '" << fn << "' >>\\n\\n";
        model.as_cdl(ofs);
    }

    argv++;
    argc--;

    try {
        size_t nsys = tset::load_nb_system(*argv, training_set);
        std::cout << "'" << *(argv++) << "' : "
                      << nsys << " configurations" << std::endl;
        if (--argc > 0) E_range = atof(*(argv++));

        if (--argc > 0) alpha = atof(*(argv++));

    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    ts_weights = std::vector<double>(training_set.size(),1.0);

    double E_min;
    size_t N_eff;
    size_t index_min;

    tset::setup_weights(training_set, E_range, E_min, index_min,
                       ts_weights, N_eff);

    std::cout << "\\n>>   E_min = " << E_min << " kcal/mol; index " << index_min <<
                 "\\n>> E_range = " << E_range << " kcal/mol" <<
                 "\\n\\n>> training set size = " << training_set.size()
              << "\\n>>    effective size = " << N_eff << '\\n'
              << std::endl;

    std::vector<double> ts_energies(training_set.size(),0.0);
    for (size_t n = 0; n < training_set.size(); n++) {
        // Saving reference energies
        ts_energies[n] = training_set[n].nb_energy ;
"""

    ff.write(a)

    # Now write the many body electrostatics calculation
    a = get_nbody_electrostatics_string(number_of_monomers, number_of_atoms, number_of_sites, "training_set[n].xyz.data()")

    ff.write(a)

    a = """
        training_set[n].nb_energy -= elec_e[n];
"""

    ff.write(a)

    a = """
        mbnrg_disp disp(training_set[n].xyz);
        disp_e.push_back(disp.get_dispersion());
        training_set[n].nb_energy -= disp_e[n];

        #ifdef USE_BUCKINGHAM
        mbnrg_buck buck(training_set[n].xyz);
        buck_e.push_back(buck.get_buckingham());
        training_set[n].nb_energy -= buck_e[n];
        #endif
"""

    if (number_of_monomers<3): ff.write(a)

    a = """
    }

    linear::allocate();

    gsl_vector* x = gsl_vector_alloc(model.get_num_nonlinear_params());
    model.get_nonlinear_parameters(x->data);

    gsl_multimin_function chisq_func;

    chisq_func.n = model.get_num_nonlinear_params();
    chisq_func.f = linear::compute_chisq;

    gsl_vector* ss = gsl_vector_alloc(model.get_num_nonlinear_params());
    gsl_vector_set_all(ss, 0.1);

    std::cout << "\\n<> initial simplex sides:\\n";
    for (size_t n = 0; n < model.get_num_nonlinear_params(); ++n)
        std::cout << n << " : " << ss->data[n] << "\\n";

    gsl_multimin_fminimizer* s =
        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2rand,
                                      model.get_num_nonlinear_params());

    gsl_multimin_fminimizer_set(s, &chisq_func, x, ss);

    int iter(0), status(0);
    //
    // Main optimization loop
    //

        do {
        ++iter;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
             break;

        const double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-4);
        //status = gsl_multimin_test_size(size, 1e-3); //changed it to test convergence

        if (status == GSL_SUCCESS)
            std::cout << "!!! converged !!!" << std::endl;

        std::cout << iter << ' ' << s->fval << ' ' << size << std::endl;
        if (iter > 0 && iter%10 == 0) {
            std::cout << "\\n<> solution:\\n";
            for (size_t n = 0; n <  model.get_num_nonlinear_params(); ++n)
                std::cout << n << " : " << s->x->data[n] << "\\n";
            std::cout << "<>" << std::endl;
        }
    } while (status == GSL_CONTINUE && iter < 5000);

    model.set_nonlinear_parameters(s->x->data);
    linear::compute_chisq(s->x, 0);
    model.set_linear_parameters(linear::params);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    //
    // report
    //

    for (size_t i = 0; i < training_set.size(); ++i)"""

    ff.write(a)

    if (number_of_monomers<3):
        a = """
        training_set[i].nb_energy += elec_e[i] + disp_e[i];

        #ifdef USE_BUCKINGHAM
        for (size_t i = 0; i < training_set.size(); ++i)
            training_set[i].nb_energy += buck_e[i];
        #endif
"""
    else:
        a = """
        training_set[i].nb_energy += elec_e[i];
"""

    ff.write(a)

    a = """

    std::ofstream correlation_file;
    correlation_file.open ("correlation.dat");
    
    correlation_file << "#" << std::setw(10) << "index"
                     << std::setw(15)       << "Reference" 
                     << std::setw(15)       << "Calculated"
                     << std::setw(15)       << "Error^2"
                     << "    \\n"; 
    
    std::ofstream terms_file;
    terms_file.open("individual_terms.dat");
    terms_file << std::setw(10) << "#index"
                     << std::setw(20)       << "Reference" 
                     << std::setw(20)       << "Calculated"
#ifdef USE_BUCKINGHAM
                     << std::setw(20)       << "Repulsion (TTM-nrg)"
#endif
                     << std::setw(20)       << "Polynomials"
"""

    ff.write(a)
    if (number_of_monomers<3):
        a = """                     << std::setw(20)       << "Dispersion"

"""
        ff.write(a)
    a = """                     << std::setw(40)       << "Electrostatics (Permanent + Induced)"
                     << "    \\n";

    double err_L2(0), err_wL2(0), err_Linf(0);
    double err_L2_lo(0), err_Linf_lo(0), nlo(0);
    double weight_sum(0);

    for (size_t i = 0; i < training_set.size(); ++i) {
"""
    ff.write(a)

    if (number_of_monomers<3):
        a = """
        double E_model = model.calculate(training_set[i].xyz) + elec_e[i] + disp_e[i];"""
        a += """
        #ifdef USE_BUCKINGHAM
        E_model += buck_e[i];
        #endif
"""
    else:
        a = """
        double E_model = model.calculate(training_set[i].xyz) + elec_e[i];"""


    ff.write(a)

    a = """
        const double delta = E_model - training_set[i].nb_energy;
        if (std::abs(delta) > err_Linf)
            err_Linf = std::abs(delta);

        correlation_file <<  std::setw(10) << i+1
                         << std::setw(15) << std::scientific
                         << std::setprecision(8)  <<  training_set[i].nb_energy
                         << std::setw(15) <<  E_model
                         << std::setw(15) <<  delta*delta  << "    \\n" ;
 
        terms_file <<  std::setw(10) << i+1
                         << std::setw(20) << std::scientific
                         << std::setprecision(8)  <<  training_set[i].nb_energy
                         << std::setw(20) <<  E_model
#ifdef USE_BUCKINGHAM
                         << std::setw(20) <<  buck_e[i]
#endif
                         << std::setw(20) << model.calculate(training_set[i].xyz)
"""

    ff.write(a)
    if (number_of_monomers<3):
        a = """                         << std::setw(20) <<  disp_e[i] 
"""
        ff.write(a)
    a = """                         << std::setw(40) <<  elec_e[i]
                         << "    \\n" ;

        err_L2 += delta*delta;
        err_wL2 += ts_weights[i]*delta*delta;
        weight_sum += ts_weights[i];

        if (training_set[i].binding_energy - E_min < E_range) {
            nlo += 1.0;
            err_L2_lo += delta*delta;
            if (std::abs(delta) > err_Linf_lo)
                err_Linf_lo = std::abs(delta);
        }
    }

    correlation_file.close();
    terms_file.close();

    err_L2 /= training_set.size();
    err_wL2 /= weight_sum;

    err_L2_lo /= nlo;

    std::cout << "      err[L2] = " << std::sqrt(err_L2) << "    #rmsd of full ts\\n"
              << "     err[wL2] = " << std::sqrt(err_wL2) << "   #weighted rmsd of full ts\\n"
              << "    err[Linf] = " << err_Linf << "   #highest error in full ts\\n"
              << "  err[L2,low] = " << std::sqrt(err_L2_lo) << "   #rmsd of low-energy ts \\n"
              << "err[Linf,low] = " << err_Linf_lo << "   #highest error in low-energy ts "
              << std::endl;

    //
    //  save
    //

    {
        const char fn[] = "fit-""" + str(number_of_monomers) + """b.cdl";
        std::ofstream ofs(fn);
        std::cout << "\\n>> saving as '" << fn << "'\\n";
        model.as_cdl(ofs);
    }


    return 0;
}

"""

    ff.write(a)
    ff.close()


def write_makefile_mbnrg(number_of_monomers, system_name):
    """
    Writes the Makefile

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
    """

    mon_files = []
    for i in range(number_of_monomers):
        mon_files.append("mon" + str(i+1) + ".o")

    fit_exes = "fit-" + str(number_of_monomers) + "b eval-" + str(number_of_monomers) + "b "
    if number_of_monomers == 2:
        fit_exes += "fit-" + str(number_of_monomers) + "b-over-ttm eval-" + str(number_of_monomers) + "b-over-ttm "

    mon_objects = " ".join(mon_files)
    ff = open("Makefile", 'w')

    escaped_name = system_name.replace("(", "_o_").replace(")", "_c_")

    a = """
ifndef INTELHOME
$(info "INTELHOME is not set. Please set it or ignore if the default /opt/intel is OK")
INTELHOME=/opt/intel
endif

ifndef GSLHOME
$(info "GSLHOME is not set. Please set it or ignore if GSL is already in your path")
GSLHOME=""
endif

ifndef NETCDFHOME
$(info "NETCDFHOMEis not set. Please set it or ignore if netcdf is already in your path")
NETCDFHOME=""
endif

CXX=icpc
CXXFLAGS= -g -Wall -std=c++11 -O0 -m64 
LIBS = -lnetcdf -lgsl -lgslcblas -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

AR = /usr/bin/ar
OBJDIR = ../obj
BINDIR = ../bin
LIBDIR = -L./ -L$(GSLHOME)/lib -L$(NETCDFHOME)/lib -L$(INTELHOME)lib/intel64 -L$(INTELHOME)/mkl/lib/intel64 
INCLUDE = -I./ -I$(GSLHOME)/include -I$(NETCDFHOME)/include -I$(INTELHOME)/mkl/include
"""
    ff.write(a)
    if (number_of_monomers<3):
        a = """
FIT_OBJ = fit-utils.o training_set.o io-xyz.o tang-toennies.o dispersion.o\\
training_set.o variable.o vsites.o water_monomer_lp.o gammq.o buckingham.o\\"""
    else:
        a = """
FIT_OBJ = fit-utils.o training_set.o io-xyz.o tang-toennies.o\\
training_set.o variable.o vsites.o water_monomer_lp.o gammq.o\\"""
    ff.write(a)
    a = """
electrostatics.o coulomb.o wlsq.o rwlsq.o ps.o mbnrg_""" + str(number_of_monomers) + "b_" + escaped_name + """_fit.o \\
""" + mon_objects + """ poly_""" + str(number_of_monomers) + "b_" + escaped_name + """_fit.o


all: libfit.a """ + fit_exes + """

libfit.a: $(addprefix $(OBJDIR)/, $(FIT_OBJ))
\t$(AR) cru libfit.a $(addprefix $(OBJDIR)/, $(FIT_OBJ))

fit-""" + str(number_of_monomers) + """b: fit-""" + str(number_of_monomers) + """b.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBDIR) $< $(LIBS) -lfit -o $@
\tmkdir -p $(BINDIR)
\tmv fit-""" + str(number_of_monomers) + """b $(BINDIR)

eval-""" + str(number_of_monomers) + """b: eval-""" + str(number_of_monomers) + """b.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBDIR) $< $(LIBS) -lfit -o $@
\tmkdir -p $(BINDIR)
\tmv eval-""" + str(number_of_monomers) + """b $(BINDIR)

fit-""" + str(number_of_monomers) + """b-over-ttm: fit-""" + str(number_of_monomers) + """b.cpp
\t$(CXX) -DUSE_BUCKINGHAM $(CXXFLAGS) $(INCLUDE) $(LIBDIR) $< $(LIBS) -lfit -o $@
\tmkdir -p $(BINDIR)
\tmv fit-""" + str(number_of_monomers) + """b-over-ttm $(BINDIR)

eval-""" + str(number_of_monomers) + """b-over-ttm: eval-""" + str(number_of_monomers) + """b.cpp
\t$(CXX) -DUSE_BUCKINGHAM $(CXXFLAGS) $(INCLUDE) $(LIBDIR) $< $(LIBS) -lfit -o $@
\tmkdir -p $(BINDIR)
\tmv eval-""" + str(number_of_monomers) + """b-over-ttm $(BINDIR)
"""
    ff.write(a)

    a = """
$(OBJDIR)/%.o: %.cpp $(OBJDIR)/.sentinel
\t$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBDIR) -c $< $(LIBS) -o $@

$(OBJDIR)/.sentinel:
\tmkdir -p $(OBJDIR)
\ttouch $@

clean:
\trm -rf $(addprefix $(OBJDIR)/, $(FIT_OBJ)) libfit*.a $(BINDIR)
"""
    ff.write(a)
    ff.close()

def write_makefile_ttmnrg(number_of_monomers, system_name):
    """
    Writes the Makefile

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
    """
    
    if (number_of_monomers != 2):
        raise PotentialFittingError("You are trying to generate a TTM-nrg fitting code for more than 2 monomers. That is not possible.")

    mon_files = []
    for i in range(number_of_monomers):
        mon_files.append("mon" + str(i+1) + ".o")

    fit_exes = "fit-" + str(number_of_monomers) + "b-ttm eval-" + str(number_of_monomers) + "b-ttm "

    mon_objects = " ".join(mon_files)
    ff = open("Makefile", 'w')

    escaped_name = system_name.replace("(", "_o_").replace(")", "_c_")

    a = """
ifndef INTELHOME
$(info "INTELHOME is not set. Please set it or ignore if the default /opt/intel is OK")
INTELHOME=/opt/intel
endif

ifndef GSLHOME
$(info "GSLHOME is not set. Please set it or ignore if GSL is already in your path")
GSLHOME=""
endif

ifndef NETCDFHOME
$(info "NETCDFHOMEis not set. Please set it or ignore if netcdf is already in your path")
NETCDFHOME=""
endif

CXX=icpc
CXXFLAGS= -g -Wall -std=c++11 -O0 -m64 
LIBS = -lnetcdf -lgsl -lgslcblas -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

AR = /usr/bin/ar
OBJDIR = ../obj
BINDIR = ../bin
LIBDIR = -L./ -L$(GSLHOME)/lib -L$(NETCDFHOME)/lib -L$(INTELHOME)lib/intel64 -L$(INTELHOME)/mkl/lib/intel64 
INCLUDE = -I./ -I$(GSLHOME)/include -I$(NETCDFHOME)/include -I$(INTELHOME)/mkl/include

FIT_OBJ = fit-utils.o training_set.o io-xyz.o tang-toennies.o dispersion.o\\
training_set.o variable.o vsites.o water_monomer_lp.o gammq.o buckingham.o\\
electrostatics.o coulomb.o wlsq.o rwlsq.o ps.o """ + mon_objects + """


all: libfit.a """ + fit_exes + """

libfit.a: $(addprefix $(OBJDIR)/, $(FIT_OBJ))
\t$(AR) cru libfit.a $(addprefix $(OBJDIR)/, $(FIT_OBJ))

fit-""" + str(number_of_monomers) + """b-ttm: fit-""" + str(number_of_monomers) + """b-ttm.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBDIR) $< $(LIBS) -lfit -o $@
\tmkdir -p $(BINDIR)
\tmv fit-""" + str(number_of_monomers) + """b-ttm $(BINDIR)

eval-""" + str(number_of_monomers) + """b-ttm: eval-""" + str(number_of_monomers) + """b-ttm.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBDIR) $< $(LIBS) -lfit -o $@
\tmkdir -p $(BINDIR)
\tmv eval-""" + str(number_of_monomers) + """b-ttm $(BINDIR)

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/.sentinel
\t$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBDIR) -c $< $(LIBS) -o $@

$(OBJDIR)/.sentinel:
\tmkdir -p $(OBJDIR)
\ttouch $@

clean:
\trm -rf $(addprefix $(OBJDIR)/, $(FIT_OBJ)) libfit*.a $(BINDIR)
"""
    ff.write(a)
    ff.close()


def write_poly_fit_header(number_of_monomers, system_name, degree, nvars, npoly):
    """
    Writes the header file for the polynomial structure

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        degree                 - The degree of the polynomial
        nvars                  - Number of variables in the polynomial
        npoly                  - Number of terms in the polynomial
    """
    system_name = system_name.replace("(", "_o_").replace(")", "_c_")

    fname = "poly_" + str(number_of_monomers) + "b_" + system_name + "_fit.h"
    ff = open(fname,'w')
    a = """
#ifndef POLY_""" + str(number_of_monomers) + "B_" + system_name + """_H
#define POLY_""" + str(number_of_monomers) + "B_" + system_name + """_H

namespace mb_system_fit {

struct poly_model{
    static const unsigned degree = """ + str(degree) + """;
    static const unsigned n_vars = """ + str(nvars) + """;

    static const unsigned size = """ + str(npoly) + """;

    void eval(const double x[n_vars], double*);
};

} // namespace

#endif

"""
    ff.write(a)
    ff.close()


def write_poly_fit_cpp(number_of_monomers, system_name, nvars, npoly, directcpp):
    """
    Writes the C++ file for the polynomial structure

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        nvars                  - Number of variables in the polynomial
        npoly                  - Number of terms in the polynomial
        directcpp              - Path to the poly-direct.cpp file
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")

    fnamecpp = "poly_" + str(number_of_monomers) + "b_" + system_name + "_fit.cpp"
    ff = open(fnamecpp,'w')
    fnameh = "poly_" + str(number_of_monomers) + "b_" + system_name + "_fit.h"
    a = """
#include \"""" + fnameh + """\"

namespace mb_system_fit {

void poly_model::eval(const double x[""" + str(nvars) + """], double a[""" + str(npoly) + """])
{
    double p[""" + str(npoly) + """];
"""
    ff.write(a)

    with open(directcpp, 'r') as fdirect:
        for line in fdirect.readlines():
            if line.startswith('    p['):
                ff.write(line)

    a = """
for(int i = 0; i < """ + str(npoly) + """; ++i)
        a[i] = p[i];

}
} // namespace mb_system
"""

    ff.write(a)
    ff.close()


def write_mbnrg_eval_code(number_of_monomers, number_of_atoms, number_of_sites, system_name):
    """
    Writes the header file for the polynomial structure

    Args:
        number_of_monomers     - Number of monomers in the system
        number_of_atoms        - Number of real atoms in the monomer
        number_of_sites        - Number of sites (electrostatic sites) of the monomer
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
    """
    system_name = system_name.replace("(", "_o_").replace(")", "_c_")

    cppname = "eval-" + str(number_of_monomers) + "b.cpp"
    ff = open(cppname,'w')
    a = """#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <vector>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <chrono>

#include "fit-utils.h"
#include "training_set.h"
#include "mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit.h"
#include "electrostatics.h"
#include "coulomb.h"
"""

    ff.write(a)

    a = """#include "dispersion.h"
#ifdef USE_BUCKINGHAM
#include "buckingham.h"
#endif
"""
    if (number_of_monomers<3): 
        ff.write(a)

    a = """
#include "constants.h"
"""

    ff.write(a)

    for i in range(1,number_of_monomers+1):
        ff.write("#include \"mon" + str(i) + ".h\"\n")

    a = """
namespace {

static mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit::polyholder_""" + str(number_of_monomers) + "b_" + system_name + """_fit model;

static std::vector<tset::nb_system> training_set;

} // namespace


int main(int argc, char** argv) {

    mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit::polyholder_""" + str(number_of_monomers) + "b_" + system_name + """_fit model;

    std::vector<tset::nb_system> training_set;
    std::vector<double> elec_e;
"""

    ff.write(a)

    a = """    std::vector<double> disp_e;
"""

    if (number_of_monomers<3): ff.write(a)

    a = """    std::vector<double> poly_e;
    std::vector<double> nb_energy;
    std::vector<double> ts_weights;


    if (argc < 3) {
        std::cerr << "usage: ./eval-""" + str(number_of_monomers) + """b  <parameters.nc> <test_set.xyz> "
                  << std::endl;
        return 0;
    }

    argv++;
    argc--;

    try {
        model.load_netcdf(*argv);

        ++argv;
        --argc;

        size_t nsys = tset::load_nb_system(*argv, training_set);
        std::cout << "'" << *(argv++) << "' : "
                      << nsys << " configurations" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    for (size_t n = 0; n < training_set.size(); n++) {
"""

    ff.write(a)

    a = get_nbody_electrostatics_string(number_of_monomers, number_of_atoms, number_of_sites, "training_set[n].xyz.data()")

    ff.write(a)

    a = """
        mbnrg_disp disp(training_set[n].xyz);
        disp_e.push_back(disp.get_dispersion());
        #ifdef USE_BUCKINGHAM
        mbnrg_buck buck(training_set[n].xyz);
        buck_e.push_back(buck.get_buckingham());
        #endif
"""

    if (number_of_monomers<3): 
        ff.write(a)

    a = """
        poly_e.push_back(model.calculate(training_set[n].xyz));
"""

    ff.write(a)

    if (number_of_monomers<3):
        a = """
        nb_energy.push_back(elec_e[n] + disp_e[n] + poly_e[n]);
        #ifdef USE_BUCKINGHAM
        nb_energy[n] += buck_e[n];
        #endif"""
    else:
        a = """
        nb_energy.push_back(elec_e[n] + poly_e[n]);"""

    ff.write(a)

    a = """
    }

    std::cerr << std::setw(9)   << "#frame["
              << std::setw(6)   << std::setfill('.') << "###"
              << std::setw(3)   << "]= " << std::setfill(' ')
              << std::setw(20)  << "Calculated"
#ifdef USE_BUCKINGHAM
              << std::setw(20)  << "Repulsion (TTM-nrg)"
#endif
              << std::setw(20)  << "Polynomials"
"""
    ff.write(a)

    if (number_of_monomers<3):
        a="""              << std::setw(20)  << "Dispersion"
"""
        ff.write(a)
    a = """              << std::setw(40)  << "Electrostatics (Permanent + Induced)"
              << std::endl;

    
    for (size_t n = 0; n < training_set.size(); n++) {
        std::cout << std::scientific << std::setprecision(8)
                  << std::setw(10)   << "frame["
                  << std::setw(6)    << std::setfill('.') << n
                  << std::setw(3)    << "]= " << std::setfill(' ')
                  << std::setw(20)   << nb_energy[n]
#ifdef USE_BUCKINGHAM
                  << std::setw(20)   << buck_e[n]
#endif
                  << std::setw(20)   << poly_e[n]
"""

    ff.write(a)

    if (number_of_monomers<3):
        a="""                  << std::setw(20)   << disp_e[n]
"""
        ff.write(a)
    a = """                  << std::setw(40)   << elec_e[n]
                  << std::endl;
    }

    return 0;
}

"""
    ff.write(a)
    ff.close()


def write_poly_header_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_in, version = "v1"):
    """
    Writes the polynomial header for MBX

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        degree                 - The degree of the polynomial
        nvars                  - Number of variables in the polynomial
        npoly                  - Number of terms in the polynomial
        version                - Will be appended to the class and files to differentiate multiple versions of the same system
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    namespace = "mbnrg_" + system_name + "_deg" + str(degree)
    struct_name = "poly_" + system_name + "_deg" + str(degree) + "_" + version
    fname = "poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + ".h"
    ff = open(fname,'w')
    a = """
#ifndef POLY_""" + str(number_of_monomers) + "B_" + namespace.upper() + """_H
#define POLY_""" + str(number_of_monomers) + "B_" + namespace.upper() + """_H

namespace """ + namespace + """ {

struct """ + struct_name + """ {
    static const unsigned degree = """ + str(degree) + """;
    static const unsigned n_vars = """ + str(nvars) + """;

    static const unsigned size = """ + str(npoly) + """;

    double eval(const double x[""" + str(nvars) + """],
              const double a[""" + str(npoly) + """]);
    double eval_direct(const double x[""" + str(nvars) + """],
                     const double a[""" + str(npoly) + """]);
    double eval(const double x[""" + str(nvars) + """],
              const double a[""" + str(npoly) + """],
                    double g[""" + str(nvars) + """]);
    double eval_direct(const double x[""" + str(nvars) + """],
                     const double a[""" + str(npoly) + """],
                           double g[""" + str(nvars) + """]);
};

} // namespace """ + namespace + """

#endif // POLY_""" + str(number_of_monomers) + "B_" + namespace.upper() + """_H

"""
    ff.write(a)

    with open(poly_in, 'r') as polyinp:
        lines = polyinp.readlines()
        a = "\n\n//Polynomial input used to generate these files:\n\n"
        for line in lines:
            a += "//  " + line
    ff.write(a)
    ff.close()
    

def retrieve_polynomial_lines(keyword_start, poly_file):
    """
    From a given polynomial file it returns a string with just the polynomial terms that start with a given keyword

    Args:
        keyword_start          - The keyword that will state if a line belongs to a polynomial term
        poly_file              - The file containing the polynomial expression
    """

    poly_lines = ""
    with open(poly_file,'r') as poly:
        line = poly.readline()
        while line != "":
            if any(line.startswith(key) for key in keyword_start):
                while True:
                    poly_lines += line
                    if ";" in line:
                        break
                    line = poly.readline()
            line = poly.readline()

    return poly_lines


def write_poly_cpp_grad_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_directory, version = "v1"):
    """
    Writes the polynomial C++ file with gradients for MBX

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        degree                 - The degree of the polynomial
        nvars                  - Number of variables in the polynomial
        npoly                  - Number of terms in the polynomial
        poly_directory         - The directory where the polynomial generation has been performed
        version                - Will be appended to the class and files to differentiate multiple versions of the same system
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    namespace = "mbnrg_" + system_name + "_deg" + str(degree)
    struct_name = "poly_" + system_name + "_deg" + str(degree) + "_" + version
    fname = "poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_grad_" + version + ".cpp"
    ff = open(fname,'w')
    a = "#include \"poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + """.h"

namespace """ + namespace + """ {

double """ + struct_name + """::eval(const double x[""" + str(nvars) + """],
            const double a[""" + str(npoly) + """],
                  double g[""" + str(nvars) + """]) {
"""

    ff.write(a)

    poly_lines = retrieve_polynomial_lines(["    const double t", "    g[", "    return"],poly_directory + "/poly-grd.cpp")

    ff.write(poly_lines)

    a = """}

} // namespace """ + namespace + """

"""

    ff.write(a)
    ff.close()

def write_direct_poly_cpp_grad_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_directory, version = "v1"):
    """
    Writes the polynomial C++ file with gradients for MBX

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        degree                 - The degree of the polynomial
        nvars                  - Number of variables in the polynomial
        npoly                  - Number of terms in the polynomial
        poly_directory         - The directory where the polynomial generation has been performed
        version                - Will be appended to the class and files to differentiate multiple versions of the same system
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    namespace = "mbnrg_" + system_name + "_deg" + str(degree)
    struct_name = "poly_" + system_name + "_deg" + str(degree) + "_" + version
    fname = "poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_grad_" + version + ".cpp"
    ff = open(fname,'w')
    a = "#include \"poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + """.h"

namespace """ + namespace + """ {

double """ + struct_name + """::eval(const double x[""" + str(nvars) + """],
            const double a[""" + str(npoly) + """],
                  double g[""" + str(nvars) + """]) {
"""

    ff.write(a)

    poly_lines = retrieve_polynomial_lines(["    double p[", "    p[", "    const double t", "    g["],poly_directory + "/poly-grd-direct.cpp")

    ff.write(poly_lines)

    a = """
    double energy(0);
    for(int i = 0; i < """ + str(npoly) + """; ++i)
        energy += p[i]*a[i];

    return energy;
}

} // namespace """ + namespace + """

"""

    ff.write(a)
    ff.close()


def write_poly_cpp_nograd_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_directory, version = "v1"):
    """
    Writes the polynomial C++ file without gradients for MBX

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        degree                 - The degree of the polynomial
        nvars                  - Number of variables in the polynomial
        npoly                  - Number of terms in the polynomial
        poly_directory         - The directory where the polynomial generation has been perfor
med
        version                - Will be appended to the class and files to differentiate multiple versions of the same system
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    namespace = "mbnrg_" + system_name + "_deg" + str(degree)
    struct_name = "poly_" + system_name + "_deg" + str(degree) + "_" + version
    fname = "poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_nograd_" + version + ".cpp"
    ff = open(fname,'w')
    a = "#include \"poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + """.h"

namespace """ + namespace + """ {

double """ + struct_name + """::eval(const double x[""" + str(nvars) + """],
            const double a[""" + str(npoly) + """]) {
"""

    ff.write(a)

    poly_lines = retrieve_polynomial_lines(["    const double t", "    g[", "    return"],poly_directory + "/poly-nogrd.cpp")

    ff.write(poly_lines)

    a = """}

} // namespace """ + namespace + """

"""

    ff.write(a)
    ff.close()

def write_direct_poly_cpp_nograd_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_directory, version = "v1"):
    """
    Writes the polynomial C++ file without gradients for MBX

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        degree                 - The degree of the polynomial
        nvars                  - Number of variables in the polynomial
        npoly                  - Number of terms in the polynomial
        poly_directory         - The directory where the polynomial generation has been perfor
med
        version                - Will be appended to the class and files to differentiate multiple versions of the same system
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    namespace = "mbnrg_" + system_name + "_deg" + str(degree)
    struct_name = "poly_" + system_name + "_deg" + str(degree) + "_" + version
    fname = "poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_nograd_" + version + ".cpp"
    ff = open(fname,'w')
    a = "#include \"poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + """.h"

namespace """ + namespace + """ {

double """ + struct_name + """::eval(const double x[""" + str(nvars) + """],
            const double a[""" + str(npoly) + """]) {
"""

    ff.write(a)

    poly_lines = retrieve_polynomial_lines(["    double p[", "    p["],poly_directory + "/poly-direct.cpp")

    ff.write(poly_lines)

    a = """
    double energy(0);
    for(int i = 0; i < """ + str(npoly) + """; ++i)
        energy += p[i]*a[i];

    return energy;

}

} // namespace """ + namespace + """

"""

    ff.write(a)
    ff.close()

def get_arguments_for_functions(arg_string, number_of_monomers):
    """
    Helper function that returns a set of arguments in a function for multiple monomers

    Args:
        arg_string             - The argument that depends on the monomer (xyz, mon, grad...)
        number_of_monomers     - The number of monomers
    """

    arg_text = ""
    for i in range(number_of_monomers):
        arg_text += arg_string + str(i+1)
        if i+1 != number_of_monomers:
            arg_text += ", "
    return arg_text


def write_mbx_polynomial_holder_header(number_of_monomers, system_name, degree, nvars, npoly, poly_directory, non_linear_parameters, ri, ro, vsites, version = "v1"):
    """
    Writes the polynomial holder header file

    Args:
        number_of_monomers     - Number of monomers in the system
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        degree                 - The degree of the polynomial
        nvars                  - Number of variables in the polynomial
        npoly                  - Number of terms in the polynomial
        poly_directory         - The directory where the polynomial generation has been perfor
med
        non_linear_parameters  - All the non linear parameters
        ri                     - Inner cutoff (not used if 1b)
        ro                     - Outer cutoff (not used in 2b)
        vsites                 - Virtual site labels
        version                - Will be appended to the class and files to differentiate multiple versions of the same system
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    namespace = "mbnrg_" + system_name + "_deg" + str(degree)
    struct_name = "mbnrg_" + system_name + "_deg" + str(degree) + "_" + version
    fname = "mbnrg_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + ".h"
    defname = fname.replace(".h","_H").upper()
    poly_header = "poly_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + ".h"
    poly_name = "poly_" + system_name + "_deg" + str(degree) + "_" + version
    ff = open(fname,'w')

    a = """#ifndef {0}
#define {0}

#include <cmath>
#include <string>
#include <vector>

#include "tools/constants.h"
#include "tools/variable.h"
#include "tools/water_monomer_lp.h"
#include \"{1}\" \n""".format(defname,poly_header)

    ff.write(a)

    # Get the monomer arguments for the constructor
    arg_constr = get_arguments_for_functions("const std::string mon",number_of_monomers)

    # Get the monomer arguments for the evals (xyz and grad)
    arg_xyz = get_arguments_for_functions("const double *xyz",number_of_monomers)
    arg_grad = get_arguments_for_functions("double *grad",number_of_monomers)

    if number_of_monomers == 1:
       declare_key = "std::vector<double> "
    else:
       declare_key = "double "

    a = """
////////////////////////////////////////////////////////////////////////////////

namespace """ + namespace + """ {

//----------------------------------------------------------------------------//

struct """ + struct_name + """ {
    """ + struct_name + """() {};
    """ + struct_name + "(" +  arg_constr + """);

    ~""" + struct_name + """() {};

    typedef """ + poly_name + """ polynomial;

    """ + declare_key + """eval(""" + arg_xyz + """, const size_t n);
    """ + declare_key + """eval(""" + arg_xyz + ", " + arg_grad + """ , const size_t n, std::vector<double> *virial=0);

  private:
"""
    ff.write(a)

    for nl in non_linear_parameters:
        ff.write("    double m_" + nl + ";\n")

    a = """
    double m_ri = """ + str(ri) + """;
    double m_ro = """ + str(ro) + """;

    double f_switch(const double, double&);

    std::vector<double> coefficients;
};

//----------------------------------------------------------------------------//

} // namespace """ + namespace + """

////////////////////////////////////////////////////////////////////////////////

#endif
"""
    ff.write(a)
    ff.close()


def write_mbx_polynomial_holder_cpp(system_name, symmetry_parser, number_of_monomers, number_of_atoms, vsites, use_lonepairs, non_linear_parameters, variables, number_of_variables, degree, version = "v1"):
    """
    Writes the polynomial holder header file

    Args:
        system_name            - The name of the system in symmetry language. Expects n fragments separated by an underscore "_" such as A1B2_C2D4
        symmetry_parser        - The SymmetryParser object that gives the symmetry of the system.
        number_of_monomers     - Number of monomers in the system
        number_of_atoms        - Number of real atoms in the monomer
        vsites                 - Virtual site labels
        use_lonepairs          - List with 0 and 1 of the same size as monomers. If 0, no lone pairs will be used or declared. If 1, lone pairs will be used. NOTE: TODO As for 09/25/2019, only water can have lone pairs.
        non_linear_parameters  - All the non linear parameters
        variables              - List of lists with the information in the poly.in file
        number_of_variables    - Number of variables
        degree                 - The degree of the polynomial
        version                - Will be appended to the class and files to differentiate multiple versions of the same system
    """

    system_name = system_name.replace("(", "_o_").replace(")", "_c_")
    namespace = "mbnrg_" + system_name + "_deg" + str(degree)
    struct_name = "mbnrg_" + system_name + "_deg" + str(degree) + "_" + version
    fname = "mbnrg_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + ".cpp"
    hname = "mbnrg_" + str(number_of_monomers) + "b_" + system_name + "_deg" + str(degree) + "_" + version + ".h"

    ff = open(fname,'w')

    # Get the monomer arguments for the constructor
    arg_constr = get_arguments_for_functions("const std::string mon",number_of_monomers)

    # Get the monomer arguments for the evals (xyz and grad)
    arg_xyz = get_arguments_for_functions("const double *xyz",number_of_monomers)
    arg_grad = get_arguments_for_functions("double *grad",number_of_monomers)

    if number_of_monomers == 1:
       declare_key = "std::vector<double> "
    else:
       declare_key = "double "

    a = '''#include "''' + hname + '''"

////////////////////////////////////////////////////////////////////////////////

namespace ''' + namespace + """ {

""" + struct_name + "::" + struct_name + "(" + arg_constr + """) {

    // =====>> BEGIN SECTION CONSTRUCTOR <<=====
    // =>> PASTE RIGHT BELOW THIS LINE <==


    // =====>> END SECTION CONSTRUCTOR <<=====
}

//----------------------------------------------------------------------------//

double """ + struct_name + """::f_switch(const double r, double& g)
{
    if (r > m_ro) {
        g = 0.0;
        return 0.0;
    } else if (r > m_ri) {
        const double t1 = M_PI/(m_ro - m_ri);
        const double x = (r - m_ri)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}

//----------------------------------------------------------------------------//

 """ + declare_key + struct_name + """::eval(""" + arg_xyz + """, const size_t n) {
    std::vector<double> energies(n,0.0);
    std::vector<double> energies_sw(n,0.0);

    std::vector<double> xyz(""" + str(sum(number_of_atoms)*3) + """);
    double sw = 0.0;
    polynomial my_poly;

    for (size_t j = 0; j < n; j++) {
"""
    ff.write(a)


    for i in range(len(number_of_atoms)):
        ff.write("        const double *mon" + str(i+1) + " = xyz" + str(i+1) + " + " + str(3*number_of_atoms[i]) + "*j;\n")
    ff.write("\n")

    all_distances = utils_nb_fitting.get_list_of_numeric_pairs("d",number_of_monomers)
    all_switches = utils_nb_fitting.get_list_of_numeric_pairs("sw",number_of_monomers)

    # Write the distances
    for d in all_distances:
        a = """
        const double """ + d[0] + """[3] =
                         {mon""" + str(d[1]) + """[0] - mon""" + str(d[2]) + """[0],
                          mon""" + str(d[1]) + """[1] - mon""" + str(d[2]) + """[1],
                          mon""" + str(d[1]) + """[2] - mon""" + str(d[2]) + """[2]};

        const double """ + d[0] + """rsq = """ + d[0] + """[0]*""" + d[0] + """[0] + """ + d[0] + """[1]*""" + d[0] + """[1] + """ + d[0] + """[2]*""" + d[0] + """[2];
        const double """ + d[0] + """r = std::sqrt(""" + d[0] + """rsq);
"""
        ff.write(a)

    # Write Check for distances but skip if distances is empty
    if len(all_distances) > 0:
        condition = "true "
        for d in all_distances:
            condition += " && " + d[0] + "r > m_ro "
    else:
        condition = "false "

    a = """
        if (""" + condition + """) {
             continue;
        }
"""
    ff.write(a)

    counter = 0
    for i in range(len(number_of_atoms)):
        ff.write("\n        std::copy(mon" + str(i+1) + ", mon" + str(i+1) + " + " + str(number_of_atoms[i]*3) + ", xyz.begin() + " + str(counter) + ");\n")
        counter += number_of_atoms[i]*3

    # Get the pointers to the atoms
    pointer_to_coordinates, pointer_to_vsites= get_pointer_setup_string(symmetry_parser, vsites, "xyz.data()", nspaces = 8)

    ff.write("\n\n" + pointer_to_coordinates)
    ff.write(pointer_to_vsites)

    a = """
        double w12 =     -9.721486914088159e-02;  //from MBpol
        double w13 =     -9.721486914088159e-02;
        double wcross =   9.859272078406150e-02;

    """
    ff.write(a)

    fragments = symmetry_parser.get_sub_parsers()

    # FIXME Only monomer that accepts lone pairs, for now, is MBpol water.
    char_code = 'a'
    for i in range(len(use_lonepairs)):
        if use_lonepairs[i] != 0:
            a = """
        monomer m""" + str(i + 1) + """;
        m""" + str(i + 1) + """.setup(coords_""" + list(fragments[i].get_atoms())[0][0] + "_1_" + char_code + """, w12, wcross, coords_""" + list(fragments[i].get_atoms())[3][0] + "_1_" + char_code + ", coords_" + list(fragments[i].get_atoms())[4][0] + "_2_" + char_code + """);
"""
            ff.write(a)
        char_code = chr(ord(char_code) + 1)

    ff.write("\n        variable vs[" + str(number_of_variables) + "];\n")
    ff.write("\n        double xs[" + str(number_of_variables) + "];\n\n")
    string_vars = get_variables_string(variables, "xs", "vs", nspaces = 8)
    ff.write(string_vars)

    # write the switches
    # for now, will be set as the sum of the products of each n-1 combination
    # for a 3b, it will be the sum of the product of each pair
    sw_string = "\n"
    for sw, d in zip(all_switches, all_distances):
        sw_string += "        double g{} = 0.0;\n".format(sw[0])
        sw_string += "        double {} = f_switch({}, g{});\n".format(sw[0],d[0] + "r" ,sw[0])

    ff.write(sw_string)

    all_switches_only = []
    for sw in all_switches:
        all_switches_only.append(sw[0])
    sw_comb = list(it.combinations(all_switches_only,number_of_monomers - 1))
    sw_string = " + "
    terms = []
    if len(sw_comb[0]) == 0:
        sw_string = "1.0"
    else:
        for swc in sw_comb:
            term = "*"
            term = term.join(swc)
            terms.append(term)

        sw_string = sw_string.join(terms)

    a = """
        sw = """ + sw_string + """;

        energies[j] = my_poly.eval(xs,coefficients.data());
        energies_sw[j] = energies[j]*sw;

    }
"""
    ff.write(a)

    if number_of_monomers == 1:
        ff.write("    return energies_sw;\n")
    else:
        a = """
    double energy = 0.0;
    for (size_t i = 0; i < n; i++)
        energy += energies_sw[i];

    return energy;
"""
        ff.write(a)


    a = """
}

//----------------------------------------------------------------------------//

""" + declare_key + struct_name + """::eval(""" + arg_xyz + ", " + arg_grad + """ , const size_t n, std::vector<double> *virial) {
    std::vector<double> energies(n,0.0);
    std::vector<double> energies_sw(n,0.0);

    std::vector<double> xyz(""" + str(sum(number_of_atoms)*3) + """);
    double sw = 0.0;
    polynomial my_poly;

    for (size_t j = 0; j < n; j++) {
"""
    ff.write(a)

    for i in range(len(number_of_atoms)):
        ff.write("        const double *mon" + str(i+1) + " = xyz" + str(i+1) + " + " + str(3*number_of_atoms[i]) + "*j;\n")

    all_distances = utils_nb_fitting.get_list_of_numeric_pairs("d",number_of_monomers)
    all_switches = utils_nb_fitting.get_list_of_numeric_pairs("sw",number_of_monomers)

    ff.write("\n")

    # Write the distances
    for d in all_distances:
        a = """
        const double """ + d[0] + """[3] =
                         {mon""" + str(d[1]) + """[0] - mon""" + str(d[2]) + """[0],
                          mon""" + str(d[1]) + """[1] - mon""" + str(d[2]) + """[1],
                          mon""" + str(d[1]) + """[2] - mon""" + str(d[2]) + """[2]};

        const double """ + d[0] + """rsq = """ + d[0] + """[0]*""" + d[0] + """[0] + """ + d[0] + """[1]*""" + d[0] + """[1] + """ + d[0] + """[2]*""" + d[0] + """[2];
        const double """ + d[0] + """r = std::sqrt(""" + d[0] + """rsq);
"""
        ff.write(a)

    # Write Check for distances but skip if distances is empty
    if len(all_distances) > 0:
        condition = "true "
        for d in all_distances:
            condition += " && " + d[0] + "r > m_ro "
    else:
        condition = "false "

    a = """
        if (""" + condition + """) {
             continue;
        }

        std::vector<double> gradients(""" + str(sum(number_of_atoms)*3) + """,0.0);
"""
    ff.write(a)

    counter = 0
    for i in range(len(number_of_atoms)):
        ff.write("\n        std::copy(mon" + str(i+1) + ", mon" + str(i+1) + " + " + str(number_of_atoms[i]*3) + ", xyz.begin() + " + str(counter) + ");\n")
        counter += number_of_atoms[i]*3

    # Get the pointers to the atoms
    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(symmetry_parser, vsites, "xyz.data()", nspaces = 8)

    ff.write(pointer_to_coordinates)
    ff.write(pointer_to_vsites)
    ff.write("\n")

    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(symmetry_parser, vsites, "gradients.data()", "g", nspaces = 8, is_const = False)

    ff.write(pointer_to_coordinates)
    ff.write(pointer_to_vsites)
    ff.write("\n")


    a = """
        double w12 =     -9.721486914088159e-02;  //from MBpol
        double w13 =     -9.721486914088159e-02;
        double wcross =   9.859272078406150e-02;

    """
    ff.write(a)

    fragments = symmetry_parser.get_sub_parsers()

    # FIXME Only monomer that accepts lone pairs, for now, is MBpol water.
    char_code = 'a'
    for i in range(len(use_lonepairs)):
        if use_lonepairs[i] != 0:
            a = """
        monomer m""" + str(i + 1) + """;
        m""" + str(i + 1) + """.setup(coords_""" + list(fragments[i].get_atoms())[0][0] + "_1_" + char_code + """, w12, wcross, coords_""" + list(fragments[i].get_atoms())[3][0] + "_1_" + char_code + ", coords_" + list(fragments[i].get_atoms())[4][0] + "_2_" + char_code + """);
"""
            ff.write(a)
        char_code = chr(ord(char_code) + 1)

    ff.write("\n        variable vs[" + str(number_of_variables) + "];\n")
    ff.write("\n        double xs[" + str(number_of_variables) + "];\n\n")
    ff.write("\n        double gxs[" + str(number_of_variables) + "];\n\n")
    string_vars = get_variables_string(variables, "xs", "vs",nspaces = 8)
    ff.write(string_vars)

    # write the switches
    # for now, will be set as the sum of the products of each n-1 combination
    # for a 3b, it will be the sum of the product of each pair
    sw_string = "\n"
    for sw, d in zip(all_switches, all_distances):
        sw_string += "        double g{} = 0.0;\n".format(sw[0])
        sw_string += "        double {} = f_switch({}, g{});\n".format(sw[0],d[0] + "r" ,sw[0])

    ff.write(sw_string)

    all_switches_only = []
    for sw in all_switches:
        all_switches_only.append(sw[0])
    sw_comb = list(it.combinations(all_switches_only,number_of_monomers - 1))
    sw_string = " + "
    terms = []
    if len(sw_comb[0]) == 0:
        sw_string = "1.0"
    else:
        for swc in sw_comb:
            term = "*"
            term = term.join(swc)
            terms.append(term)

        sw_string = sw_string.join(terms)

    a = """
        sw = """ + sw_string + """;

        energies[j] = my_poly.eval(xs,coefficients.data(),gxs);
        energies_sw[j] = energies[j]*sw;

        for (size_t i = 0; i < """ + str(number_of_variables) + """; i++) {
            gxs[i] *= sw;
        }

"""
    ff.write(a)

    string_vars = get_grad_var_string(variables, "vs", "gxs", nspaces = 8)
    ff.write(string_vars)
    
    

    fragments = symmetry_parser.get_sub_parsers()

    # FIXME Only monomer that accepts lone pairs, for now, is MBpol water.
    char_code = 'a'
    for i in range(len(use_lonepairs)):
        if use_lonepairs[i] != 0:
            a = """
        m""" + str(i+1) + """.grads(coords_""" + list(fragments[i].get_atoms())[3][0] + "_1_" + char_code + "_g, coords_" + list(fragments[i].get_atoms())[4][0] + "_2_" + char_code + """_g, w12, wcross, coords_""" + list(fragments[i].get_atoms())[0][0] + "_1_" + char_code + """_g);
"""
            ff.write(a)
        char_code = chr(ord(char_code) + 1)

    # Now finalize switch gradients
    all_switch_summands = sw_string.split(" + ")
    for sw, d in zip(all_switches, all_distances):
        switch_grad = []
        for summand in all_switch_summands:
            if sw[0] in summand:
                switch_grad.append(summand.replace(sw[0],"1.0"))

        switch_text = " + ".join(switch_grad)
        ff.write("        g" + sw[0] + " *= (" + switch_text + ")*energies[j]/" + d[0] + "r;\n")

    a = """

        for (size_t i = 0; i < 3; i++) {
"""

    ff.write(a)

    # And put the switch grads into the atoms
    counter = 0
    for i in range(len(number_of_atoms)):
        string_to_add = "            gradients[" + str(counter) + " + i] += 0.0 "
        current_mon = i + 1
        for sw, d in zip(all_switches, all_distances):
            if current_mon == sw[1]:
                string_to_add += "+ (g" + sw[0] + "*" + d[0] + "[i])"
            elif current_mon == sw[2]:
                string_to_add += "- (g" + sw[0] + "*" + d[0] + "[i])"

        string_to_add += ";\n"
        ff.write(string_to_add)

        counter += number_of_atoms[i]*3

    ff.write("        }\n")
    count = 0
    for i in range(len(number_of_atoms)):
        a = """

        for (size_t i = 0; i < """ + str(3*number_of_atoms[i]) + """; i++) {
            grad""" + str(i+1) + """[i + j*""" + str(3*number_of_atoms[i]) + """] += gradients[""" + str(count) + """ + i];
        }
"""
        count += 3*number_of_atoms[i]
        ff.write(a)
## EL put virial Here
        a = """
        
        if (virial != 0) {"""
    
    ff.write(a)

    pairdict = {  0 : [0,0],
                  1 : [0,1],
                  2 : [0,2],
                  4 : [1,1],
                  5 : [1,2],
                  8 : [2,2] }

    for virial_index in [0,1,2,4,5,8]:



        a = """
        
            (*virial)[""" + str(virial_index) + """] += """

        ff.write(a)


        parser_list=list(symmetry_parser.get_atoms())
        mod_parser_list=[]
        for param in parser_list:
            if param[0] not in ['X','Y','Z']:
               
                mod_parser_list.append(param)


        for labelindex in list(range(0,len(mod_parser_list))):
            pos_name=get_coords_var_name(mod_parser_list[labelindex][0], mod_parser_list[labelindex][1], mod_parser_list[labelindex][2])
            grad_name=get_coords_var_name(mod_parser_list[labelindex][0], mod_parser_list[labelindex][1], mod_parser_list[labelindex][2], "g")
             
            if labelindex==0:
                a = """-"""+pos_name + """[""" + str(pairdict[virial_index][0]) + """]*""" + grad_name + """[""" + str(pairdict[virial_index][1]) + """]\n"""
            elif labelindex==len(mod_parser_list)-1:
                a = """                        -""" + pos_name + """[""" + str(pairdict[virial_index][0]) + """]*""" + grad_name + """[""" + str(pairdict[virial_index][1]) + """];\n"""
            else:
                a = """                        -""" + pos_name + """[""" + str(pairdict[virial_index][0]) + """]*""" + grad_name + """[""" + str(pairdict[virial_index][1]) + """]\n"""   
            ff.write(a)
    a = """
            (*virial)[3] = (*virial)[1];
            (*virial)[6] = (*virial)[2];
            (*virial)[7] = (*virial)[5];
"""
    
    ff.write(a)


    a = """

        }
"""
    ff.write(a)
          
    a = """

    }
"""
    ff.write(a)

    if number_of_monomers == 1:
        ff.write("    return energies_sw;\n")
    else:
        a = """
    double energy = 0.0;
    for (size_t i = 0; i < n; i++)
        energy += energies_sw[i];

    return energy;
"""
        ff.write(a)

    a = """
}

//----------------------------------------------------------------------------//
} // namespace """ + namespace + """
"""

    ff.write(a)
    ff.close()
