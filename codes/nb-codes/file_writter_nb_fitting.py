import utils_nb_fitting
import itertools as it

def write_monomer_class_header(mon_index):
    filename = "mon" + str(mon_index) + ".h"
    ff = open(filename,'w')

    a="""#ifndef MON""" + str(mon_index) + """_H
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

    ff.write(a)
    ff.close()

def write_monomer_class_cpp(mon_index, nsites, natoms, excl12, excl13, excl14, chg, pol, polfac):
    filename = "mon" + str(mon_index) + ".cpp"
    ff = open(filename,'w')

    a = """#include "mon""" + str(mon_index) + """.h"

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

    ff.write(a)
    for p in excl12:
        ff.write('    excluded12.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
    for p in excl13:
        ff.write('    excluded13.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
    for p in excl14:
        ff.write('    excluded14.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
    a = """

  }

  double* mon""" + str(mon_index) + """::set_charges(double* atmcrds) {
    charge = memory;
"""
    ff.write(a)
    for i in range(len(chg)):
        ff.write('    charge[' + str(i) + '] = ' + chg[i] + '*constants::CHARGECON;\n')
    a = """

    return charge;
  }
  
  double* mon""" + str(mon_index) + """::set_pol() {
    atmpolar = memory + nsites + nsites*3;
"""
    ff.write(a)
    for i in range(len(pol)):
        ff.write('    atmpolar[' + str(i) + '] = ' + pol[i] + ';\n')
    a = """
    return atmpolar;
  }

  double* mon""" + str(mon_index) + """::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;
    
"""
    ff.write(a)
    for i in range(len(polfac)):
        ff.write('    polfac[' + str(i) + '] = ' + polfac[int(i)] + ';\n')
    a = """
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

    ff.write(a)
    ff.close()

def write_mbpol_monomer(mon_index):
    filename = "mon" + str(mon_index) + ".cpp"
    ff = open(filename,'w')

    a = """
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
    ff.write(a)
    ff.close()

def write_fit_polynomial_holder_header(system_name, number_of_monomers, non_linear_parameters, ri, ro):
    defflag = "MBNRG_" + str(number_of_monomers) + "B_" + system_name + "_FIT"
    system_keyword_mbnrg = "mbnrg_" + str(number_of_monomers) + "b_" + system_name
    system_keyword_poly = "poly_" + str(number_of_monomers) + "b_" + system_name
    system_keyword_polyholder = "polyholder_" + str(number_of_monomers) + "b_" + system_name
    header_mbnrg_fit = system_keyword_mbnrg + "_fit.h"
    header_poly_fit = system_keyword_poly + "_fit.h"

    # Generate arguments for the eval
    args_eval = "const std::vector<double> mon1"
    for i in range(2,number_of_monomers + 1):
        args_eval += ", const std::vector<double> mon" + str(i)

    ff = open(header_mbnrg_fit,'w')
    a = """#ifndef """ + defflag + """
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

    ff.write(a)
    a=""
    for nl_param in non_linear_parameters:
        a += "    double m_" + nl_param + ";\n"

    ff.write(a)

    a = """protected:
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

    ff.write(a)
    ff.close()

def get_individual_atoms_with_type(monomer_atom_types, vsites):
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

def get_pointer_setup_string(monomer_atom_types, vsites, xyz_var_name):
    crd_shift = 0
    mon_id = 'a'
    string_pointers = ""
    string_pointers_vs = ""
    for monomer in monomer_atom_types:
        num_ats = int(len(monomer)/2 + 0.49999)
        for i in range(num_ats):
            atom = monomer[2*i]
            if not atom in vsites:
                for j in range(monomer[2*i + 1]):
                    string_pointers += "    const double* " + atom + "_" + str(j+1) + "_" + mon_id + " = " + xyz_var_name + " + " + str(crd_shift) + ";\n"
                    crd_shift += 3
            else:
                for j in range(monomer[2*i + 1]):
                    string_pointers_vs += "    double " + atom + "_" + str(j+1) + "_" + mon_id + "[3];\n"
        string_pointers += "\n"
        string_pointers_vs += "\n"
        mon_id = chr(ord(mon_id)+1)    
    return string_pointers, string_pointers_vs
            

def get_variables_string(variables, name_vout, name_vstruct):
    variables_string = ""

    for i in range(len(variables)):
        var = variables[i]
        types_1 = utils_nb_fitting.get_atom_types(var[0])
        types_2 = utils_nb_fitting.get_atom_types(var[2])
        atom1_name = "{}_{}_{}".format(types_1[0], types_1[1], var[1])
        atom2_name = "{}_{}_{}".format(types_2[0], types_2[1], var[3])
        
        sorted_atoms = "".join(sorted([types_1[0], types_2[0]]))

        kconst = "m_k_"
        dconst = "m_d_"
        if var[1] == var[3]:
            kconst += "intra_"
            dconst += "intra_"
 
        argument_constants = kconst + sorted_atoms
        if "0" in var[-1]:
            argument_constants = dconst + sorted_atoms + ', ' + kconst

        variables_string += "    {}[{}] = {}[{}].v_{}({}, {}, {});\n".format(name_vout,i,name_vstruct,i,var[-1],argument_constants,atom1_name,atom2_name)

    return variables_string     


def write_fit_polynomial_holder_cpp(system_name, monomer_atom_types, number_of_monomers, number_of_atoms, vsites, use_lonepairs, non_linear_parameters, variables, number_of_variables, ri, ro, k_min_intra, k_max_intra, k_min, k_max, d_min_intra, d_max_intra, d_min, d_max):
    system_keyword_mbnrg = "mbnrg_" + str(number_of_monomers) + "b_" + system_name
    system_keyword_poly = "poly_" + str(number_of_monomers) + "b_" + system_name
    system_keyword_polyholder = "polyholder_" + str(number_of_monomers) + "b_" + system_name
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
    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(monomer_atom_types, vsites, "xyz.data()")
    
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

    # FIXME Only monomer that accepts lone pairs, for now, is MBpol water.
    char_code = 'a'
    for i in range(len(use_lonepairs)):
        if use_lonepairs[i] != 0:
            a = """
    monomer m""" + str(i+1) + """;
    m""" + str(i+1) + """.setup(""" + monomer_atom_types[i][0] + "_1_" + char_code + """, w12, wcross, """ + monomer_atom_types[i][4] + "_1_" + char_code + ", " + monomer_atom_types[i][4] + "_2_" + char_code + """);
"""
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


def write_dispersion_header(monomer_atom_types, virtual_sites_poly, c6, d6):
    hname = "dispersion.h"
    ff = open(hname,'w')
    a = """
#ifndef DISPERSION_H
#define DISPERSION_H

#include "tang-toennies.h"
#include <cmath>
#include <algorithm>



struct mbnrg_disp {
  mbnrg_disp();
  mbnrg_disp(std::vector<double> c);
  ~mbnrg_disp();

"""
    ff.write(a)

    # Need to select if we have 1b or 2b
    if len(monomer_atom_types) == 1:
        pairs = utils_nb_fitting.get_nonbonded_pairs(virtual_sites_poly,monomer_atom_types[0])
    else:
        pairs = utils_nb_fitting.get_nonbonded_pairs(virtual_sites_poly,monomer_atom_types[0],monomer_atom_types[1])

    # Write C6 declaration
    a = ""
    b = ""
    for pair in pairs:
        a += "  const double m_C6_" + pair + " = " + str(c6[pairs.index(pair)]) + ";\n"
        b += "  const double m_d6_" + pair + " = " + str(d6[pairs.index(pair)]) + ";\n"

    ff.write(a + "\n")
    ff.write(b)

    a = """
  const double m_C8 = 0.0;
  const double m_d8 = 0.0;

  const double if6 = 1.0/x2o::factorial<6>();
  const double if8 = 1.0/x2o::factorial<8>();

  std::vector<double> xyz;

  double get_dispersion();


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


def write_dispersion_cpp(monomer_atom_types, vsites, excl12 = None, excl13 = None, excl14 = None):
    cppname = "dispersion.cpp"
    ff = open(cppname,'w')
    a = """#include "dispersion.h"

mbnrg_disp::mbnrg_disp() {}
mbnrg_disp::~mbnrg_disp() {}

mbnrg_disp::mbnrg_disp(std::vector<double> c) {
  xyz = c;
}

double mbnrg_disp::get_dispersion() {

    double disp = 0.0;
  
"""
    ff.write(a)

    # Get pointer to coordinates only for real sites
    pointer_to_coordinates, pointer_to_vsites = get_pointer_setup_string(monomer_atom_types, vsites, "xyz.data()")

    ff.write(pointer_to_coordinates)

    # Write the dispersion calculations that are excluded
    atom_types = get_individual_atoms_with_type(monomer_atom_types, vsites)
    # Case in which we have a monomer
    a = ""
    if len(monomer_atom_types) == 1:
        mon = atom_types[0]
        for i in range(len(mon)-1):
            for j in range(i+1,len(mon)):
                if [i,j] not in excl12 and [i,j] not in excl13 and [i,j] not in excl14:
                     a += "    disp += x6(m_C6_{}{}, m_d6_{}{}, m_C8, m_d8, {}_{}_{}, {}_{}_{});\n".format(mon[i][0], mon[j][0], mon[i][0], mon[j][0], mon[i][0], mon[i][1], mon[i][2], mon[j][0], mon[j][1], mon[j][2])

    elif len(monomer_atom_types) == 2:
        mon1 = atom_types[0]
        mon2 = atom_types[1]
        for i in range(len(mon1)):
            for j in range(len(mon2)):
                 a += "    disp += x6(m_C6_{}{}, m_d6_{}{}, m_C8, m_d8, {}_{}_{}, {}_{}_{});\n".format(mon1[i][0], mon2[j][0], mon1[i][0], mon2[j][0], mon1[i][0], mon1[i][1], mon1[i][2], mon2[j][0], mon2[j][1], mon2[j][2])

    ff.write(a)

    a = """
    return disp;
}
        
"""

    ff.write(a)
    ff.close()

def get_nl_params_initialization_string(nl_param_all, k_min, k_max, d_min, d_max, k_min_intra, k_max_intra, d_min_intra, d_max_intra):
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
    electrostatics_string = ""

    # generate first index list to locate monomers
    first_index_of_atoms = []
    first_index_of_sites = []
    fi = 0
    fi_s = 0
    for i in range(number_of_monomers):
        first_index_of_atoms.append(fi)
        first_index_of_sites.append(fi)
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

def write_fitting_code(number_of_monomers, number_of_atoms, number_of_sites, system_name, nl_param_all, k_min, k_max, d_min, d_max, k_min_intra, k_max_intra, d_min_intra, d_max_intra):
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
#include "dispersion.h"
#include "constants.h"
"""

    ff.write(a)

    for i in range(1,number_of_monomers+1):
        ff.write("#include \"mon" + str(i) + ".h\"\n")

    a = """
namespace {

static mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit::polyholder_""" + str(number_of_monomers) + "b_" + system_name + """_fit model;

double alpha = 0.0005;

double E_range = 20.0; // kcal/mol

//----------------------------------------------------------------------------//

static std::vector<tset::nb_system> training_set;
static std::vector<double> elec_e;
static std::vector<double> disp_e;
static std::vector<double> ts_weights;

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

    nl_params_string = get_nl_params_initialization_string(nl_param_all, k_min, k_max, d_min, d_max, k_min_intra, k_max_intra, d_min_intra, d_max_intra)

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

        mbnrg_disp disp(training_set[n].xyz);
        disp_e.push_back(disp.get_dispersion());

        training_set[n].nb_energy -= disp_e[n];
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

    for (size_t i = 0; i < training_set.size(); ++i)
        training_set[i].nb_energy += elec_e[i] + disp_e[i];
        

    std::ofstream correlation_file;
    correlation_file.open ("correlation.dat");
    double err_L2(0), err_wL2(0), err_Linf(0);
    double err_L2_lo(0), err_Linf_lo(0), nlo(0);

    for (size_t i = 0; i < training_set.size(); ++i) {
        
        const double E_model = model.calculate(training_set[i].xyz) + elec_e[i] + disp_e[i];
        const double delta = E_model - training_set[i].nb_energy;
        if (std::abs(delta) > err_Linf)
            err_Linf = std::abs(delta);

        correlation_file <<  std::setw(10) << i+1   
                         << std::setw(15) << std::scientific 
                         << std::setprecision(4)  <<  training_set[i].nb_energy
                         << std::setw(15) <<  E_model  
                         << std::setw(15) <<  delta*delta  << "    \\n" ;

        err_L2 += delta*delta;
        err_wL2 += ts_weights[i]*delta*delta;

        if (training_set[i].binding_energy - E_min < E_range) {
            nlo += 1.0;
            err_L2_lo += delta*delta;
            if (std::abs(delta) > err_Linf_lo)
                err_Linf_lo = std::abs(delta);
        }
    }

    correlation_file.close();

    err_L2 /= training_set.size();
    err_wL2 /= training_set.size();

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

def write_makefile(number_of_monomers, system_name):

    mon_files = []
    for i in range(number_of_monomers):
        mon_files.append("mon" + str(i+1) + ".o")

    mon_objects = " ".join(mon_files)
    ff = open("Makefile", 'w')
    a = """
CXX=icpc
CXXFLAGS= -g -Wall -std=c++11 -O0 -m64 -I/opt/intel/mkl/include
LIBS = -lnetcdf -lgsl -lgslcblas -L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

AR = /usr/bin/ar
OBJDIR = .
LIBDIR = ./
INCLUDE = -I./

FIT_OBJ = fit-utils.o training_set.o io-xyz.o tang-toennies.o dispersion.o\\
training_set.o variable.o vsites.o water_monomer_lp.o gammq.o \\
electrostatics.o coulomb.o rwlsq.o ps.o mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit.o \\
""" + mon_objects + """ poly_""" + str(number_of_monomers) + "b_" + system_name + """_fit.o

all: libfit.a fit-""" + str(number_of_monomers) + """b eval-""" + str(number_of_monomers) + """b

libfit.a: $(addprefix $(OBJDIR)/, $(FIT_OBJ))
\t$(AR) cru libfit.a $(addprefix $(OBJDIR)/, $(FIT_OBJ))

fit-""" + str(number_of_monomers) + """b: fit-""" + str(number_of_monomers) + """b.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) $< $(LIBS) -lfit -o $@

eval-""" + str(number_of_monomers) + """b: eval-""" + str(number_of_monomers) + """b.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) $< $(LIBS) -lfit -o $@

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/.sentinel
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) -c $< $(LIBS) -o $@

$(OBJDIR)/.sentinel:
\tmkdir -p $(OBJDIR)
\ttouch $@

clean:
\trm -rf $(addprefix $(OBJDIR)/, $(FIT_OBJ)) libfit*.a fit-""" + str(number_of_monomers) + """b
"""
    ff.write(a)
    ff.close()


def write_poly_fit_header(number_of_monomers, system_name, degree, nvars, npoly):
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

def write_eval_code(number_of_monomers, number_of_atoms, number_of_sites, system_name):
    cppname = "eval-" + str(number_of_monomers) + "b.cpp"
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

#include "fit-utils.h"
#include "training_set.h"
#include "mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"
#include "constants.h"
"""

    ff.write(a)

    for i in range(1,number_of_monomers+1):
        ff.write("#include \"mon" + str(i) + ".h\"\n")

    a = """
namespace {

static mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit::polyholder_""" + str(number_of_monomers) + "b_" + system_name + """_fit model;

static std::vector<tset::nb_system> training_set;
static std::vector<double> elec_e;
static std::vector<double> disp_e;
static std::vector<double> nb_energy;

} // namespace


int main(int argc, char** argv) {

    mbnrg_""" + str(number_of_monomers) + "b_" + system_name + """_fit::polyholder_""" + str(number_of_monomers) + "b_" + system_name + """_fit model;

    std::vector<tset::nb_system> training_set;
    std::vector<double> elec_e;
    std::vector<double> disp_e;
    std::vector<double> poly_e;
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

        poly_e.push_back(model.calculate(training_set[n].xyz));

        nb_energy.push_back(elec_e[n] + disp_e[n] + poly_e[n]);
    }
    
    for (size_t n = 0; n < training_set.size(); n++) {
        std::cout << std::scientific << std::setprecision(8)
                  << std::setw(10) << "frame[" 
                  << std::setw(6) << std::setfill('.') << n
                  << std::setw(3) << "]= " << std::setfill(' ') 
                  << std::setw(25) << nb_energy[n] << std::endl;
    }

    return 0;
}

"""
    ff.write(a)
    ff.close()












