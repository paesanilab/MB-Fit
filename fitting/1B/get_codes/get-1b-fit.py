#!/usr/bin/env python
# coding: utf-8

# In[78]:

import sys
import os
import json


# In[79]:

# Check proper if input is provided


# In[80]:

if len(sys.argv) != 4:
    print("Usage: ./script <input.in> <poly-direct.cpp_with_path> <config.ini path>\n ")
    sys.exit()
else:
    name = sys.argv[1]
    directcpp = sys.argv[2]
    config_filename = sys.argv[3]

import configparser

config = configparser.ConfigParser()
config.read(config_filename)


# In[81]:

# This should be the commandline argument
#name = "A1B4.in"
f = open(name, 'r')
mon1 = f.readline().split('\'')[1]


# In[82]:

# This should be the second command line argument
#directcpp = 'poly-direct.cpp'


# In[83]:

# For Andrea:
# Find a way to find the number of atoms in each monomer
# Store them in nat1
nat = config.getint("fitting", "number_of_atoms")

# Find the number of sites (electrostatic sites)
nsites = config.getint("fitting", "number_of_electrostatic_sites")

# Obtain the lists with the excluded pairs
excl12 = json.loads(config.get("fitting", "excluded_pairs_12"))
excl13 = json.loads(config.get("fitting", "excluded_pairs_13"))
excl14 = json.loads(config.get("fitting", "excluded_pairs_14"))


#Obtain charges (in the order of input), pols and polfacs
chg = json.loads(config.get("fitting", "charges"))
pol = json.loads(config.get("fitting", "polarizabilities"))
polfac = json.loads(config.get("fitting", "polarizability_fractions"))

#Ask the user for the max value of k and d
k_min = config.getfloat("fitting", "k_min")
k_max = config.getfloat("fitting", "k_max")

d_min = config.getfloat("fitting", "d_min")
d_max = config.getfloat("fitting", "d_max")


# Obtain C6 and d6 from user in the same order as the given pairs AA, AB ...:
# Look at the input to retrieve this
#These are the INTERMOLECULAR C6 for same m
C6 = json.loads(config.get("fitting", "C6"))
d6 = json.loads(config.get("fitting", "d6"))

# FInd a way to get the degree
degree = config.getint("common", "polynomial_order")
# Find a way to get the number ov variables
nvars = 10
# Find a way to get the size of the polynomial
npoly = 82

# Define kind of variables 
# 'coul' is e^(-k(d-d0)/d)
# 'exp'  is e^(-k(d-d0))
var = config.get("fitting", "var")

# Define Energy Range for the fitting
E_range = config.getfloat("fitting", "energy_range")

# Define list of variables that are fictitious
vsites = json.loads(config.get("fitting", "virtual_sites_labels"))


# In[84]:

types = list(mon1)


# In[85]:

# Generating the non linear parameter list
nlparam = []
t1 = []
# Monomer 1 parameters
for i in range(0,len(types),2):
    for j in range(int(types[i+1])):
        t1.append(types[i])

print(t1)


# In[86]:

# Appending for mon1
for i in range(len(t1)):
    if t1[i] in vsites:
        continue
    for j in range(1,len(t1)):
        if t1[j] in vsites:
            continue
        for nlp in ['d','k']:
            const_intra = nlp + '_' + t1[i] + t1[j]
            if not const_intra in nlparam:
                nlparam.append(const_intra)

# Getting the pairs
pairs = []
real_pairs = []

# First, the excluded pairs (12,13,14)
excluded_pairs = []
for i in range(len(excl12)):
    x = t1[excl12[i][0]]
    y = t1[excl12[i][1]]
    if not x in vsites and not y in vsites:
        p = sorted([x,y])
        ps = p[0] + p[1]
        if not ps in excluded_pairs:
            excluded_pairs.append(ps)
            
for i in range(len(excl13)):
    x = t1[excl13[i][0]]
    y = t1[excl13[i][1]]
    if not x in vsites and not y in vsites:
        p = sorted([x,y])
        ps = p[0] + p[1]
        if not ps in excluded_pairs:
            excluded_pairs.append(ps)
            
for i in range(len(excl14)):
    x = t1[excl14[i][0]]
    y = t1[excl14[i][1]]
    if not x in vsites and not y in vsites:
        p = sorted([x,y])
        ps = p[0] + p[1]
        if not ps in excluded_pairs:
            excluded_pairs.append(ps)
            
for i in range(len(t1)):
    for j in range(i+1,len(t1)):
        p = sorted([t1[i],t1[j]])
        ps = p[0] + p[1]
        if not ps in pairs:
            pairs.append(ps)
        if not p[0] in vsites and not p[1] in vsites and not ps in real_pairs:
            real_pairs.append(ps)

print(nlparam)
print(pairs)
print(real_pairs)
print(excluded_pairs)


# In[87]:

# Save number in num_nonlinear
num_nonlinear = len(nlparam)


# ## Creating mon1.h

# In[88]:

mon1_class = open('mon1.h','w')

a="""
#ifndef MON1_H
#define MON1_H
#include "molecule.h"

namespace x   {

   class mon1 : public molecule {
      public :
      mon1();
      ~mon1();
      
      mon1(double* crd);
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

mon1_class.write(str(a))
mon1_class.close()


# ## Creating their cpp files

# In[89]:

ff = open('mon1.cpp','w')

a = """
#include "mon1.h"
#include <algorithm>
#include <cstddef>
#include "fit-constants.h"

//#define VERBOSE
#ifdef VERBOSE 
#include <iostream>
#endif

const double CHARGECON = constants::CHARGECON;

// NOTE: ASSUMING ABB

namespace x  {

  mon1::mon1() { }

  mon1::~mon1() {
    delete[] memory;
  }

  mon1::mon1( double* crd) {
    
    nsites = """ + str(nsites) + """; 
    realsites = """ + str(nat) + """;
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

  double* mon1::set_charges(double* atmcrds) {
    charge = memory;
"""
ff.write(a)
for i in range(len(chg)):
    ff.write('    charge[' + str(i) + '] = ' + str(chg[i]) + '*CHARGECON;\n')
a = """

    return charge;
  }
  
  double* mon1::set_pol() {
    atmpolar = memory + nsites + nsites*3;
"""
ff.write(a)
for i in range(len(pol)):
    ff.write('    atmpolar[' + str(i) + '] = ' + str(pol[i]) + ';\n')
a = """
    return atmpolar;
  }

  double* mon1::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;
    
"""
ff.write(a)
for i in range(len(polfac)):
    ff.write('    polfac[' + str(i) + '] = ' + str(polfac[i]) + ';\n')
a = """
    return polfac;
  }

  double* mon1::set_sitecrds(double* atmcrds) {
    sitecrds = memory + nsites;
    
// ##DEFINE HERE## ficticious site coordinates
    return atmcrds;  
  }

  void mon1::allocate() {
    memory = new double [nsites  //charge
      + nsites*3  //site coordinates
      + nsites  //polarizabilities 
      + nsites];  //polfac
    }

    int mon1::get_nsites() { return nsites; }
    int mon1::get_realsites() { return realsites; }
    double* mon1::get_charges() { return charge; }
    double* mon1::get_sitecrds() { return sitecrds; }
    double* mon1::get_pol() { return atmpolar; }
    double* mon1::get_polfacs() { return polfac; }

    excluded_set_type::iterator mon1::get_begin_12() { return excluded12.begin(); }
    excluded_set_type::iterator mon1::get_begin_13() { return excluded13.begin(); }
    excluded_set_type::iterator mon1::get_begin_14() { return excluded14.begin(); }
    excluded_set_type::iterator mon1::get_end_12() { return excluded12.end(); }
    excluded_set_type::iterator mon1::get_end_13() { return excluded13.end(); }
    excluded_set_type::iterator mon1::get_end_14() { return excluded14.end(); }

} // namespace x
"""

ff.write(a)
ff.close()


# ## Create training_set.h/cpp files

# In[90]:

ff = open('training_set.h','w')
a = """
#ifndef TRAINING_SET_H
#define TRAINING_SET_H

#include <cstddef>
#include <vector>

namespace tset {

struct monomer{
    double energy_onebody;
    double xyz[""" + str((nat)*3) + """];
};

size_t load_monomers(const char* filename, std::vector<monomer>& ts);

} // namespace tset

#endif // TRAINING_SET_H
"""
ff.write(a)
ff.close()


# ## X1B h file

# In[91]:

hname = "x1b_" + mon1 + "_v1.h"
polyhname = "poly_1b_" + mon1 + "_v1x.h"
polyhname2 = "poly_1b_" + mon1 + ".h"
defname = "X1B_" + mon1 + "_V1_H"
ff = open(hname,'w')
ff.write('#ifndef ' + defname + '\n')
ff.write('#define ' + defname + '\n \n')
ff.write('#include "' + polyhname + '" \n')
ff.write('#include "' + polyhname2 + '" \n \n')

a = """
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace x1b_""" + mon1 + """ {

//----------------------------------------------------------------------------//

"""
ff.write(a)
ff.write('struct x1b_' + mon1 + "_v1x { \n")
a = """
    typedef mb_system::poly_model poly_type;
    typedef mb_system_fit::poly_model poly_type_fit;


    static std::string name();
    void load_netcdf(const char*);

    // returns 1B contribution only
    // XYZ is for the real sites
"""
ff.write(a)
ff.write('    double operator()(const double xyz[' + str(3*(nat)) + ']) const; \n')
a = """
    double eval(const double* mon1) const;
    double eval(const double* mon1, double* g1) const;
    void cart_to_vars(const double xyz[""" + str(3*(nat)) + """], double* vars) const;
    
    // For fitting purposes
    size_t get_nvars() {return poly_type::n_vars;}
    size_t get_num_nonlinear_params() {return """ + str(num_nonlinear) + """;}
    void set_nonlinear_parameters(const double*);
    void set(const double* xxx);
    void get_nonlinear_parameters(double*);
    bool nonlinear_parameters_out_of_range() const;
    inline double basis(const double xyz[""" + str(3*(nat)) + """], double*) const; //basis polynomials
    inline size_t nparams() {
      return poly_type::size;
    }
    

    inline void as_cdl(std::ostream&) const;
    void write_cdl(std::ostream&, unsigned, unsigned, const double*) const;
    
private:

"""
ff.write(a)
for nl in nlparam:
    ff.write("    double m_" + nl + ";\n")
a = """

private:
    double m_poly[poly_type::size];
};

// For fitting
inline void x1b_""" + mon1 + """_v1x::as_cdl(std::ostream& os) const
{
    write_cdl(os, """ + str(degree) + """, poly_type::size, m_poly);
}


inline double
x1b_""" + mon1 + """_v1x::basis(const double xyz[""" + str(3*(nat)) + """], double* mmm) const 
{
    double v[poly_type::n_vars];

    cart_to_vars(xyz, v);

    poly_type_fit polyn;
    polyn.eval(v, mmm);

    return 0; // we will substract electrostatics and dispersion in the main
}

//----------------------------------------------------------------------------//

} // namespace x1b_""" + mon1 + """

////////////////////////////////////////////////////////////////////////////////

#endif 
"""
ff.write(a)
ff.close()


# ## CPP file

# In[92]:

cppname = "x1b_" + mon1 + "_v1.cpp"
ff = open(cppname,'w')
a = """
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>

#include <netcdf.h>

#include "stuff.h"
#include "constants.h"
"""
hname = "x1b_" + mon1 + "_v1.h"
ff.write(a)
ff.write('#include "' + hname + '" \n \n')
a = """
////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int kode) {

    std::cerr << " ** Fatal Error in x1b_""" + mon1 + """_v1x::load_netcdf() ** :" 
              << nc_strerror(kode) << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

struct variable {
    double v_exp(const double& r0, const double& k,
                 const double * p1, const double * p2 );

    double v_coul(const double& r0, const double& k,
                  const double * p1, const double * p2 );

    double g[3]; // diff(value, p1 - p2)
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_exp(const double& r0, const double& k,
                       const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);

    const double exp1 = std::exp(k*(r0 - r));
    const double gg = - k*exp1/r;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return exp1;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_coul(const double& r0, const double& k,
                        const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    const double r = std::sqrt(rsq);

    const double exp1 = std::exp(k*(r0 - r));
    const double rinv = 1.0/r;
    const double val = exp1*rinv;

    const double gg = - (k + rinv)*val*rinv;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return val;
}



//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

//struct vsites {
//    //void TwoParticleAverageSite() {}
//    //void ThreeParticleAverageSite() {}
//    void OutOfPlaneSite(const double& w12, const double& w13,
//                        const double& wcross, const double x1[3],
//                        const double y1[3], const double y2[3],
//                        double vs[3]);
//    //void LocalCoordinatesSite{}
//};
//
//void vsites::OutOfPlaneSite(const double& w12,
//                            const double& w13,
//                            const double& wcross,
//                            const double x1[3],
//                            const double y1[3],
//                            const double y2[3],
//                            double vs[3]) {
//    double r12[3], r13[3];
//
//    for (int i = 0; i < 3; ++i) {
//        r12[i] = y1[i] - x1[i];
//        r13[i] = y2[i] - x1[i];
//    }
//                            
//    double rc[3];
//    rc[0] = r12[1]*r13[2] - r12[2]*r13[1];
//    rc[1] = r12[2]*r13[0] - r12[0]*r13[2];
//    rc[2] = r12[0]*r13[1] - r12[1]*r13[0];
//    
//    vs[0] = x1[0] + w12 * r12[0] + w13 * r13[0] + wcross * rc[0];
//    vs[1] = x1[1] + w12 * r12[1] + w13 * r13[1] + wcross * rc[1];
//    vs[2] = x1[2] + w12 * r12[2] + w13 * r13[2] + wcross * rc[2];
//}

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace x1b_""" + mon1 + """ {

//----------------------------------------------------------------------------//


std::string x1b_""" + mon1 + """_v1x::name() {
    return "x1b_""" + mon1 + """_v1x";
}
"""
ff.write(a)
ff.write('void x1b_' + mon1 + '_v1x::set_nonlinear_parameters(const double* xxx) \n { \n ')
for nl in nlparam:
    ff.write("    m_" + nl + " = *xxx++; \n")
a = """
}

//----------------------------------------------------------------------------//

"""
ff.write(a)
ff.write('void x1b_' + mon1 + '_v1x::get_nonlinear_parameters(double* xxx) \n { \n ')
for nl in nlparam:
    ff.write("   *xxx++ = m_" + nl + "; \n")
a = """
}

//----------------------------------------------------------------------------//

void x1b_""" + mon1 + """_v1x::set(const double* xxx) {
  std::copy(xxx, xxx + nparams(), m_poly);
}

//----------------------------------------------------------------------------//


void x1b_""" + mon1 + """_v1x::write_cdl(std::ostream& os, unsigned DEG,
                              unsigned npoly, const double* poly) const
{
    using namespace std;

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific
       << "netcdf x1b_""" + mon1 + """_v1x {" << endl
       << "  // global attributes " << endl
       << "  :name = \\"x1b_""" + mon1 + """_v1x<" << DEG << ">\\";" << endl;

    // x1b_""" + mon1 + """_v1x::as_cdl(os);
"""
ff.write(a)
ff.write('    os ')
for nl in nlparam:
    ff.write('       << "  :' + nl + ' = " << setw(22) << m_' + nl +  ' << "; // A^(-1))" << endl \n')
a = """
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
ff.write('bool x1b_' + mon1 + '_v1x::nonlinear_parameters_out_of_range() const { \n')
ff.write('    const double k_min =  ' + str(k_min) + ' ;\n')
ff.write('    const double k_max =  ' + str(k_max) + ' ;\n')

ff.write('    const double d_min =  ' + str(d_min) + ' ;\n')
ff.write('    const double d_max =  ' + str(d_max) + ' ;\n')

ff.write('return false')
for nl in nlparam:
    if nl.startswith('d'):
        ff.write('\n       || m_' + nl + ' < d_min ')
        ff.write('\n       || m_' + nl + ' > d_max ')
    else:
        ff.write('\n       || m_' + nl + ' < k_min ')
        ff.write('\n       || m_' + nl + ' > k_max ')
    
ff.write('; \n')            
a = """

}


"""
ff.write(a)
ff.write('void  x1b_' + mon1 + '_v1x::cart_to_vars(const double* xyz, double* v) const { \n')
ff.write('    // NOTE: XYZ contains ONLY the real sites. The lone pairs etc are calculated here \n')

nc = 0
set_m1 = []
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if not types[i] in vsites:
            ff.write('    const double* ' + types[i] + '_' + str(n) + '_a' + '= xyz + ' + str(3 * nc) + ';\n')
            set_m1.append(types[i] + '_' + str(n) + '_a')
            n = n + 1
            nc = nc + 1
ff.write('\n')
            
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if types[i] in vsites:
            ff.write('    double ' + types[i] + '_' + str(n) + '_a[3]' + ';\n')
            set_m1.append(types[i] + '_' + str(n) + '_a')
            n = n + 1
ff.write('\n')

a = """
    
    variable vr[""" + str(nvars) + """];
    using x2o::distance;
    
"""

ff.write(a)
nv = 0
# Intramolecular distances:
for i in range(0,len(set_m1) - 1):
    for j in range(i + 1,len(set_m1)):
        ti = set_m1[i].split('_')[0]
        tj = set_m1[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        if not ti in vsites and not tj in vsites:
            ff.write('    v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var + '(m_d_' + t + ', m_k_' + t + ', ' + set_m1[i] + ', ' + set_m1[j] + ');\n')
            nv = nv + 1
ff.write('\n')


a = """
    
#define PR(x)
"""
ff.write(a)
for i in range(nv):
    ff.write('  PR(v[' + str(i) + ']);\n')
a = """
} 

//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

void x1b_""" + mon1 + """_v1x::load_netcdf(const char* fn)
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
for nl in nlparam:
    ff.write("    RETRIEVE(" + nl + ")\n")
a = """

#   undef RETRIEVE

    int varid;

    if ((rc = nc_inq_varid(ncid, "poly", &varid)))
        error(rc);

    for (size_t n = 0; n < poly_type::size; ++n) {
        if ((rc = nc_get_var1_double(ncid, varid, &n, m_poly + n)))
            error(rc);
    }

    if ((rc = nc_close(ncid)))
        error(rc);
}


//----------------------------------------------------------------------------//

double x1b_""" + mon1 + """_v1x::eval(const double* mon1) const
{

    double xcrd[""" + str(3*(nat)) + """]; // coordinates of real sites ONLY

    std::copy(mon1, mon1 + """ + str(3*(nat)) + """, xcrd);
    
    double v[""" + str(nvars) + """]; 

    cart_to_vars(xcrd, v); 
    
    const double E_poly = mb_system::poly_model::eval(m_poly, v);

    return E_poly;
}

double x1b_""" + mon1 + """_v1x::operator()(const double crd[""" + str(3*(nat)) + """]) const
{
    const double E_poly = eval(crd);

    return E_poly;
}

} // namespace x1b_""" + mon1 + """

////////////////////////////////////////////////////////////////////////////////
"""
ff.write(a)
ff.close()


# ## Dispersion.h

# In[93]:

hname = "dispersion.h"
ff = open(hname,'w')
a = """
#ifndef DISPERSION_H
#define DISPERSION_H

#include "tang-toennies.h"
#include <cmath>
#include <algorithm>



struct x1b_disp {
  x1b_disp();
  x1b_disp(double *, size_t);
  ~x1b_disp();

"""
ff.write(a)
for i in range(len(real_pairs)):
    if not real_pairs[i] in excluded_pairs:
        ff.write('  const double m_C6_' + real_pairs[i] + ' = ' + str(C6[i]) + ' ; \n')
for i in range(len(real_pairs)):
    if not real_pairs[i] in excluded_pairs:
        ff.write('  const double m_d6_' + real_pairs[i] + ' = ' + str(d6[i]) + ' ; \n')
    
a = """

  const double m_C8 = 0.0;
  const double m_d8 = 0.0;

  const double if6 = 1.0/x2o::factorial<6>();
  const double if8 = 1.0/x2o::factorial<8>();

  double * xyz1;

  double get_dispersion();
  double get_dispersion(double * grdx);


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


# ## dispersion.cpp

# In[94]:

cppname = "dispersion.cpp"
ff = open(cppname,'w')
a = """
#include "dispersion.h"

x1b_disp::x1b_disp() {
  xyz1 = new double[3];
}
x1b_disp::~x1b_disp() {
  delete[] xyz1;
}

x1b_disp::x1b_disp(double * c1, size_t n1) {
  xyz1 = new double[3*n1];
  std::copy(c1, c1 + 3*n1, xyz1);
}

double x1b_disp::get_dispersion() {

  double disp = 0.0;
"""
ff.write(a)

nc = 0
set_m1 = []
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if not types[i] in vsites:
            ff.write('  const double* ' + types[i] + '_' + str(n) + ' = xyz1 + ' + str(3 * nc) + ';\n')
            set_m1.append(types[i] + '_' + str(n))
            n = n + 1
            nc = nc + 1
ff.write('\n')
   
for i in range(0,len(set_m1) - 1):
    for j in range(i + 1,len(set_m1)):
        ti = set_m1[i].split('_')[0]
        tj = set_m1[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        if not t in excluded_pairs:
            ff.write('  disp += x6(m_C6_' + t + ', m_d6_' + t + ', m_C8, m_d8, ' + set_m1[i] + ', ' + set_m1[j] + ');\n')
            
    ff.write('\n')
a = """
  return disp;
}

double x1b_disp::get_dispersion(double * grd) {

  double disp = 0.0;
"""
ff.write(a)

nc = 0
set_m1 = []
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if not types[i] in vsites:
            ff.write('  const double* ' + types[i] + '_' + str(n) + ' = xyz1 + ' + str(3 * nc) + ';\n')
            set_m1.append(types[i] + '_' + str(n))
            n = n + 1
            nc = nc + 1
ff.write('\n')

nc = 0
set_m1 = []
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if not types[i] in vsites:
            ff.write('  const double* ' + types[i] + '_' + str(n) + '_g = grd + ' + str(3 * nc) + ';\n')
            set_m1.append(types[i] + '_' + str(n))
            n = n + 1
            nc = nc + 1
ff.write('\n')



   
for i in range(0,len(set_m1) - 1):
    for j in range(i + 1,len(set_m1)):
        ti = set_m1[i].split('_')[0]
        tj = set_m1[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        if not t in excluded_pairs:
            ff.write('  disp += x6(m_C6_' + t + ', m_d6_' + t + ', m_C8, m_d8, \n             '                 + set_m1[i] + ', ' + set_m1[j] + set_m1[i] + '_g, ' + set_m1[j] + '_g, ' + ');\n')
        
    ff.write('\n')
a = """
  return disp;
}
"""
ff.write(a)
ff.close()


# ## Fitting routine

# In[95]:

cppname = "fit-1b.cpp"
ff = open(cppname,'w')
a = """
#include <cmath>
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

#define RIDGE_REGRESSION 1

#ifdef RIDGE_REGRESSION
#   include "rwlsq.h"
#else
#   include "wlsq.h"
#endif

#include "mon1.h"
#include "fit-utils.h"
#include "training_set.h"
#include "x1b_""" + mon1 + """_v1.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"


#define dont_be_VERBOSE yes

//#define DEBUG 

#ifdef DEBUG
#define PR(x) std::cout << #x << ": " << (x) << std::endl;
#else
#define PR(x)
#endif /* DEBUG */

using namespace std;

namespace {

static x1b_""" + mon1 + """::x1b_""" + mon1 + """_v1x model;


#ifdef RIDGE_REGRESSION
const double alpha = 0.0005;
#endif

// ##DEFINE HERE## energy range
const double E_range = """ + str(E_range) + """; // kcal/mol

//----------------------------------------------------------------------------//

static std::vector<tset::monomer> training_set;
static std::vector<double> elec_e;
static std::vector<double> disp_e;
static double* ts_weights = 0;

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

    PR(model_type::nparams())
    for (size_t n = 0; n < training_set.size(); ++n) {
        double mmm[model.nparams()];
        y[n] = training_set[n].energy_onebody
               - model.basis(training_set[n].xyz, mmm);
        for (size_t p = 0; p < model.nparams(); ++p) {
            A[p + n*model.nparams()] = mmm[p];
        }
    }

#   ifdef VERBOSE
    std::cout << "=== calling wlsq::solve() ["
              << kit::wlsq::implementation()
              << "] ===" << std::endl;
#   endif

    double chisq;

//    try {
#     ifdef RIDGE_REGRESSION
      double penaltysq;
      kit::rwlsq::solve(training_set.size(), model.nparams(),
                        A, y, ts_weights, alpha, params, chisq, penaltysq);

      std::cout << "<#> chisq = " << chisq
                << " : penaltysq = " << penaltysq
                << std::endl;
#     else
      int rank;
      kit::wlsq::solve(training_set.size(), model_type::nparams(),
                       A, y, ts_weights, params, chisq, rank);
      std::cout << "<#> chisq = " << chisq
                << " : rank = " << rank
                << std::endl;
#     endif

      if (!gsl_finite (chisq)) {
        return 10000000.0;
      }
//    } catch (const std::exception& e) {
//      return 1000000.0;
//    }
#   ifdef VERBOSE
    std::cout << "\\n--> chisq = " << chisq
              << "\\n-->  rank = " << rank
              << '\\n' << std::endl;
#   endif

    return chisq;
}

//----------------------------------------------------------------------------//

} // namespace linear

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: fit-1b ts1 ..."
                  << std::endl;
        return 0;
    }
      
    long long int duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch()).count();

    srand(duration);

    double x0[""" + str(len(nlparam)) + """];
"""
ff.write(a)
for i in range(len(nlparam)):
    if nlparam[i].startswith('d'):
        ff.write('      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(d_max) - float(d_min)) + ' + ' + str(d_min) + ';\n')
    else:
        ff.write('      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(k_max) - float(k_min)) + ' + ' + str(k_min) + ';\n')

            
a = """
      

    model.set_nonlinear_parameters(x0);

    #   ifdef RIDGE_REGRESSION
    std::cout << "<> using ridge regression with alpha = "
              << alpha << std::endl;
#   endif

    std::cout << "\\n<><><> model type = '" << model.name() << "'\\n";

    {
        const char fn[] = "fit-1b-initial.cdl";
        std::ofstream ofs(fn);
        std::cout << "\\n>> dumping initial model as '" << fn << "' >>\\n\\n";
        model.as_cdl(ofs);
    }

    ++argv;

    try {
        while (--argc != 0) {
            size_t nd = tset::load_monomers(*argv, training_set);
            std::cout << "'" << *(argv++) << "' : "
                      << nd << " monomers loaded" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    //
    // assign weights
    //

    ts_weights = new double[training_set.size()];

    double E_min;
    size_t N_eff;

    tset::setup_weights(training_set, E_range, E_min,
                       ts_weights, N_eff, true);

    std::cout << "\\n>>   E_min = " << E_min << " kcal/mol"
                 "\\n>> E_range = " << E_range << " kcal/mol"
                 "\\n\\n>> training set size = " << training_set.size()
              << "\\n>>    effective size = " << N_eff << '\\n'
              << std::endl;


    // electrostatics
    double tb_ref[training_set.size()];
    for (size_t n = 0; n < training_set.size(); n++) {
      // Saving reference 1b energies  
      tb_ref[n] = training_set[n].energy_onebody ;

      x::mon1 m1(training_set[n].xyz);

      int system_nsites = m1.get_nsites();
      int * system_is_w;
      double* system_sitecrds;
      double* system_charge;
      double* system_polfac;
      double* system_pol;

      system_sitecrds = new double [system_nsites*3];
      system_charge = new double [system_nsites];
      system_polfac = new double [system_nsites];
      system_pol = new double [system_nsites]; //allocates memory to the pointers
      system_is_w = new int[system_nsites];
      
      std::fill(system_is_w, system_is_w + m1.get_nsites(), m1.is_w);

      std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
                system_sitecrds);

      std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
                system_charge);

      std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
                system_pol);

      std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
                system_polfac);


      excluded_set_type exclude12;
      excluded_set_type exclude13;
      excluded_set_type exclude14;

      for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
        exclude12.insert(*i);
      }
      

      for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
        exclude13.insert(*i);
      }
      

      for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
        exclude14.insert(*i);
      }
      

      ttm::electrostatics m_electrostatics;

      ttm::smear_ttm4x smr; 
      smr.m_aDD_intra_12 = 0.3;
      smr.m_aDD_intra_13 = 0.3;
      smr.m_aDD_intra_14 = 0.055;

      double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
                                    system_sitecrds, exclude12, exclude13, exclude14, 
                                    system_is_w, smr, 0);

      // Take out electrostatic energy:
      elec_e.push_back(ener);
      training_set[n].energy_onebody -= ener;
      // std::cerr << "Conf " << n << " : Elec= " << ener ;

      // Now need to take out dispersion
      x1b_disp disp(m1.get_sitecrds(), m1.get_realsites());
      ener = disp.get_dispersion();
      disp_e.push_back(ener);
      training_set[n].energy_onebody -= ener ;

      // std::cerr << " , Disp= " << ener << std::endl;
      delete[] system_sitecrds ;
      delete[] system_charge  ;
      delete[] system_polfac  ;
      delete[] system_pol ;
      delete[] system_is_w;

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
    model.set(linear::params);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    //
    // report
    //

    for (size_t i = 0; i < training_set.size(); ++i)
      training_set[i].energy_onebody += elec_e[i] + disp_e[i];

    ofstream correlation_file;
    correlation_file.open ("correlation.dat");
    double err_L2(0), err_wL2(0), err_Linf(0);
    double err_L2_lo(0), err_Linf_lo(0), nlo(0);

    for (size_t i = 0; i < training_set.size(); ++i) {

        const double E_model = model(training_set[i].xyz) + elec_e[i] + disp_e[i];
        const double delta = E_model - tb_ref[i];
        if (std::abs(delta) > err_Linf)
            err_Linf = std::abs(delta);

  correlation_file <<  i+1   << "   "
    << E_model  << "   "
    <<  tb_ref[i]  << "    "
          <<  delta*delta  << "    \\n" ;

        err_L2 += delta*delta;
        err_wL2 += ts_weights[i]*delta*delta;

        if (training_set[i].energy_onebody - E_min < E_range) {
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
        const char fn[] = "fit-1b.cdl";
        std::ofstream ofs(fn);
        std::cout << "\\n>> saving as '" << fn << "'\\n";
        model.as_cdl(ofs);
    }


    return 0;
}

"""
ff.write(a)
ff.close()


# ## Makefile

# In[103]:

fname = "Makefile"
ff = open(fname,'w')
a = """
CXX=g++
CXXFLAGS= -g -Wall -std=c++11 -O0 -m64 -I/opt/intel/mkl/include
LIBS = -lnetcdf -lgsl -lgslcblas -L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

AR = /usr/bin/ar
OBJDIR = .
LIBDIR = ./
INCLUDE = -I./

FIT_OBJ = fit-utils.o coulomb.o gammq.o electrostatics.o io-xyz.o \\
kvstring.o mon1.o rwlsq.o wlsq.o stuff.o tang-toennies.o \\
training_set.o poly_1b_""" + mon1 + """_v1x.o \\
x1b_""" + mon1 + """_v1.o poly_1b_""" + mon1 + """_v1.o \\
dispersion.o poly_1b_""" + mon1 + """.o 

EVAL_OBJ = fit-utils.o coulomb.o electrostatics.o gammq.o io-xyz.o \\
kvstring.o mon1.o  rwlsq.o wlsq.o stuff.o tang-toennies.o \\
training_set.o poly_1b_""" + mon1 + """_v1x.o \\
x1b_""" + mon1 + """_v1x.o poly_1b_""" + mon1 + """_v1.o \\
dispersion.o poly_1b_""" + mon1 + """.o

all: libfit.a libeval.a fit-1b eval-1b

libfit.a: $(addprefix $(OBJDIR)/, $(FIT_OBJ))
\t$(AR) cru libfit.a $(addprefix $(OBJDIR)/, $(FIT_OBJ))

libeval.a: $(addprefix $(OBJDIR)/, $(EVAL_OBJ))
\t$(AR) cru libeval.a $(addprefix $(OBJDIR)/, $(EVAL_OBJ))

fit-1b: fit-1b.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) $< $(LIBS) -lfit -o $@

eval-1b: eval-1b.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) $< $(LIBS) -leval -o $@

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/.sentinel
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) -c $< $(LIBS) -o $@

$(OBJDIR)/.sentinel:
\tmkdir -p $(OBJDIR)
\ttouch $@

clean:
\trm -rf $(addprefix $(OBJDIR)/, $(FIT_OBJ)) libfit.a libeval fit-1b eval-1b
"""
ff.write(a)
ff.close()


# ## Poly_fit_header

# In[97]:

fname = "poly_1b_" + mon1 + ".h"
ff = open(fname,'w')
a = """
#ifndef POLY_1B_""" + mon1 + """_H
#define POLY_1B_""" + mon1 + """_H

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


# ## Modifi poly-direct.cpp file to give just the non linear terms

# In[98]:

fdirect = open(directcpp, 'r')
fnamecpp = "poly_1b_" + mon1 + ".cpp"
fpolycpp = open(fnamecpp,'w')
fnameh = "poly_1b_" + mon1 + ".h"
a = """
#include \"""" + fnameh + """\"

namespace mb_system_fit {

void poly_model::eval(const double x[""" + str(nvars) + """], double a[""" + str(npoly) + """])
{
    double p[""" + str(npoly) + """];
"""
fpolycpp.write(a)

for line in fdirect.readlines():
    if line.startswith('    p['):
        fpolycpp.write(line)
    
a = """
for(int i = 0; i < """ + str(npoly) + """; ++i)
        a[i] = p[i];

}
} // namespace mb_system
"""

fpolycpp.write(a)
fpolycpp.close()


# ## Evaluation code

# In[108]:

ff = open('eval-1b.cpp','w')
a = """
#include <cmath>
#include <cassert>
#include <cstdlib>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "mon1.h"
#include "training_set.h"
#include "x1b_""" + mon1 + """_v1x.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"
#include "io-xyz.h"

#define GRADIENTS

static std::vector<double> elec_e;
static std::vector<double> disp_e;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: eval fit-1b.nc monomer.xyz"
                  << std::endl;
        return 0;
    }
    std::cout << std::scientific << std::setprecision(9);
    x1b_""" + mon1 + """::x1b_""" + mon1 + """_v1x pot;
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
    
    double xyz[""" + str(3*(nat)) + """];
    std::copy(crd.begin(), crd.end(), xyz);
    
    x::mon1 m1(xyz);
    
    int system_nsites = m1.get_nsites();
    int * system_is_w;
    double* system_sitecrds;
    double* system_charge;
    double* system_polfac;
    double* system_pol;

    system_sitecrds = new double [system_nsites*3];
    system_charge = new double [system_nsites];
    system_polfac = new double [system_nsites];
    system_pol = new double [system_nsites]; //allocates memory to the pointers
    system_is_w = new int[system_nsites];
      
    std::fill(system_is_w, system_is_w + m1.get_nsites(), m1.is_w);

    std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
              system_sitecrds);

    std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
              system_charge);

    std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
              system_pol);

    std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
              system_polfac);

    excluded_set_type exclude12;
    excluded_set_type exclude13;
    excluded_set_type exclude14;

    for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
        exclude12.insert(*i);
    }

    for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
        exclude13.insert(*i);
    }

    for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
        exclude14.insert(*i);
    }

    ttm::electrostatics m_electrostatics;
    ttm::smear_ttm4x smr; 
    smr.m_aDD_intra_12 = 0.3;
    smr.m_aDD_intra_13 = 0.3;
    smr.m_aDD_intra_14 = 0.055;

    double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
                                   system_sitecrds, exclude12, exclude13, exclude14, 
                                   system_is_w, smr, 0);

    elec_e.push_back(ener);
      
    // Now need to take out dispersion
    x1b_disp disp(m1.get_sitecrds(), m1.get_realsites());
    ener = disp.get_dispersion();
    disp_e.push_back(ener);
    
    double Epoly = pot(xyz);
    std::cout << "E_nograd = " << Epoly + elec_e[0] + disp_e[0] << std::endl;
    std::cout << "E_poly = " << Epoly << std::endl;
    std::cout << "E_elec = " << elec_e[0] << std::endl;
    std::cout << "E_disp = " << disp_e[0] << std::endl;
    
#ifdef GRADIENTS
    const double eps = 1.0e-5;
    double grd[""" + str(3*(nat)) + """];
    std::fill(grd, grd + """ + str(3*(nat)) + """, 0.0);
    Epoly = pot(xyz, grd);
    std::cout << "E_grad = " << Epoly << std::endl;
    for(size_t n = 0; n < """ + str(3*(nat)) + """; ++n){
      const double x_orig = xyz[n];

      xyz[n] = x_orig + eps;
      const double Ep = pot(xyz) ;
      std::cout << "E_p = " << Ep << std::endl;

      xyz[n] = x_orig + 2*eps;
      const double E2p = pot(xyz) ;
      std::cout << "E_2p = " << E2p << std::endl;

      xyz[n] = x_orig - 2*eps;
      const double E2m = pot(xyz) ;
      std::cout << "E_2m = " << E2m << std::endl;

      xyz[n] = x_orig - eps;
      const double Em = pot(xyz) ;
      std::cout << "E_m = " << Em << std::endl;

      const double gfd = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
      xyz[n] = x_orig;

      std::cout << elements[n/3] << "   "  << "Analit: " << grd[n] << " Numerical: " << gfd
                     << " Diff: " << std::fabs(grd[n] - gfd) << '\\n';
    }
#endif

    // Free memory
    delete[] system_sitecrds ;
    delete[] system_charge  ;
    delete[] system_polfac  ;
    delete[] system_pol ;
    delete[] system_is_w;
    
    return 0;
}
"""
ff.write(a)
ff.close()


# ## X1B.h for software

# In[106]:

hname = "x1b_" + mon1 + "_v1x.h"
polyhname = "poly_1b_" + mon1 + "_v1x.h"
defname = "X1b_" + mon1 + "_V1X_H"
ff = open(hname,'w')
ff.write('#ifndef ' + defname + '\n')
ff.write('#define ' + defname + '\n \n')
ff.write('#include "' + polyhname + '" \n')

a = """
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace x1b_""" + mon1 + """ {

//----------------------------------------------------------------------------//

"""
ff.write(a)
ff.write('struct x1b_' + mon1 + "_v1x { \n")
a = """
    typedef mb_system::poly_model poly_type;


    static std::string name();
    void load_netcdf(const char*);

    // returns 1b contribution only
    // XYZ is for the real sites
"""
ff.write(a)
ff.write('    double operator()(const double xyz[' + str(3*(nat)) + '], double grd[' + str(3*(nat)) + ']) const; \n')
ff.write('    double operator()(const double xyz[' + str(3*(nat)) + ']) const; \n')
a = """
    double eval(const double* mon1) const;
    double eval(const double* mon1, double* g1) const;
    
private:

"""
ff.write(a)
for nl in nlparam:
    ff.write("    double m_" + nl + ";\n")
a = """

private:
    double m_poly[poly_type::size];
};

//----------------------------------------------------------------------------//

} // namespace x1b_""" + mon1 + """

////////////////////////////////////////////////////////////////////////////////

#endif 
"""
ff.write(a)
ff.close()


# ## X1B.cpp for software

# In[111]:

cppname = "x1b_" + mon1 + "_v1x.cpp"
ff = open(cppname,'w')
a = """
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>

#include <netcdf.h>

"""
hname = "x1b_" + mon1 + "_v1x.h"
ff.write(a)
ff.write('#include "' + hname + '" \n \n')
a = """
////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int kode) {

    std::cerr << " ** Fatal Error in x1b_""" + mon1 + """_v1x::load_netcdf() ** :" 
              << nc_strerror(kode) << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

struct variable {
    double v_exp(const double& r0, const double& k,
                 const double * p1, const double * p2 );

    double v_coul(const double& r0, const double& k,
                  const double * p1, const double * p2 );
                  
    void grads(const double& gg, double * grd1, double * grd2,
               const double * p1, const double * p2);

    double g[3]; // diff(value, p1 - p2)
};

//----------------------------------------------------------------------------//

void variable::grads(const double& gg, double * grd1, double * grd2, 
                     const double * p1, const double * p2) {
    for (size_t i = 0; i < 3 ; i++) {
        double d = gg*g[i];
        grd1[i] += d;
        grd2[i] -= d;
    }
}

//----------------------------------------------------------------------------//

double variable::v_exp(const double& r0, const double& k,
                       const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);

    const double exp1 = std::exp(k*(r0 - r));
    const double gg = - k*exp1/r;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return exp1;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_coul(const double& r0, const double& k,
                        const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    const double r = std::sqrt(rsq);

    const double exp1 = std::exp(k*(r0 - r));
    const double rinv = 1.0/r;
    const double val = exp1*rinv;

    const double gg = - (k + rinv)*val*rinv;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return val;
}

struct monomer {
    double oh1[3];
    double oh2[3];

    void setup(const double* ohh,
               const double& in_plane_g, const double& out_of_plane_g,
               double x1[3], double x2[3]);

    void grads(const double* g1, const double* g2,
               const double& in_plane_g, const double& out_of_plane_g,
               double* grd) const;
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::setup(const double* ohh,
                    const double& in_plane_g, const double& out_of_plane_g,
                    double* x1, double* x2)
{
    for (int i = 0; i < 3; ++i) {
        oh1[i] = ohh[i + 3] - ohh[i];
        oh2[i] = ohh[i + 6] - ohh[i];
    }

    const double v[3] = {
        oh1[1]*oh2[2] - oh1[2]*oh2[1],
        oh1[2]*oh2[0] - oh1[0]*oh2[2],
        oh1[0]*oh2[1] - oh1[1]*oh2[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double in_plane = ohh[i] + 0.5*in_plane_g*(oh1[i] + oh2[i]);
        const double out_of_plane = out_of_plane_g*v[i];

        x1[i] = in_plane + out_of_plane;
        x2[i] = in_plane - out_of_plane;
    }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::grads(const double* g1, const double* g2,
                    const double& in_plane_g, const double& out_of_plane_g,
                    double* grd) const
{
    const double gm[3] = {
        g1[0] - g2[0],
        g1[1] - g2[1],
        g1[2] - g2[2]
    };

    const double t1[3] = {
        oh2[1]*gm[2] - oh2[2]*gm[1],
        oh2[2]*gm[0] - oh2[0]*gm[2],
        oh2[0]*gm[1] - oh2[1]*gm[0]
    };

    const double t2[3] = {
        oh1[1]*gm[2] - oh1[2]*gm[1],
        oh1[2]*gm[0] - oh1[0]*gm[2],
        oh1[0]*gm[1] - oh1[1]*gm[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double gsum = g1[i] + g2[i];
        const double in_plane = 0.5*in_plane_g*gsum;

        const double gh1 = in_plane + out_of_plane_g*t1[i];
        const double gh2 = in_plane - out_of_plane_g*t2[i];

        grd[i + 0] += gsum - (gh1 + gh2); // O
        grd[i + 3] += gh1; // H1
        grd[i + 6] += gh2; // H2
    }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

//struct vsites {
//    //void TwoParticleAverageSite() {}
//    //void ThreeParticleAverageSite() {}
//    void OutOfPlaneSite(const double& w12, const double& w13,
//                        const double& wcross, const double x1[3],
//                        const double y1[3], const double y2[3],
//                        double vs[3]);
//    //void LocalCoordinatesSite{}
//};
//
//void vsites::OutOfPlaneSite(const double& w12,
//                            const double& w13,
//                            const double& wcross,
//                            const double x1[3],
//                            const double y1[3],
//                            const double y2[3],
//                            double vs[3]) {
//    double r12[3], r13[3];
//
//    for (int i = 0; i < 3; ++i) {
//        r12[i] = y1[i] - x1[i];
//        r13[i] = y2[i] - x1[i];
//    }
//                            
//    double rc[3];
//    rc[0] = r12[1]*r13[2] - r12[2]*r13[1];
//    rc[1] = r12[2]*r13[0] - r12[0]*r13[2];
//    rc[2] = r12[0]*r13[1] - r12[1]*r13[0];
//    
//    vs[0] = x1[0] + w12 * r12[0] + w13 * r13[0] + wcross * rc[0];
//    vs[1] = x1[1] + w12 * r12[1] + w13 * r13[1] + wcross * rc[1];
//    vs[2] = x1[2] + w12 * r12[2] + w13 * r13[2] + wcross * rc[2];
//}

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace x1b_""" + mon1 + """ {

//----------------------------------------------------------------------------//


std::string x1b_""" + mon1 + """_v1x::name() {
    return "x1b_""" + mon1 + """_v1x";
}

"""
ff.write(a)
a = """
//----------------------------------------------------------------------------//

void x1b_""" + mon1 + """_v1x::load_netcdf(const char* fn)
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
for nl in nlparam:
    ff.write("    RETRIEVE(" + nl + ")\n")
a = """



#   undef RETRIEVE

    int varid;

    if ((rc = nc_inq_varid(ncid, "poly", &varid)))
        error(rc);

    for (size_t n = 0; n < poly_type::size; ++n) {
        if ((rc = nc_get_var1_double(ncid, varid, &n, m_poly + n)))
            error(rc);
    }

    if ((rc = nc_close(ncid)))
        error(rc);
}

//----------------------------------------------------------------------------//

double x1b_""" + mon1 + """_v1x::eval(const double* mon1) const
{

    double xcrd[""" + str(3*(nat)) + """]; // coordinates of real sites ONLY

    std::copy(mon1, mon1 + """ + str(3*(nat)) + """, xcrd);
    
    double v[""" + str(nvars) + """];
    
"""
ff.write(a)
nc = 0
set_m1 = []
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if not types[i] in vsites:
            ff.write('    const double* ' + types[i] + '_' + str(n) + '= xcrd + ' + str(3 * nc) + ';\n')
            set_m1.append(types[i] + '_' + str(n))
            n = n + 1
            nc = nc + 1

ff.write('\n')            
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if types[i] in vsites:
            ff.write('    double ' + types_a[i] + '_' + str(n) + '[3]' + ';\n')
            set_m1.append(types[i] + '_' + str(n))
            n = n + 1
ff.write('\n')

a = """    
    variable vr[""" + str(nvars) + """];
    
"""

ff.write(a)
nv = 0
# Intramolecular distances:
for i in range(0,len(set_m1) - 1):
    for j in range(i + 1,len(set_m1)):
        ti = set_m1[i].split('_')[0]
        tj = set_m1[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        if not ti in vsites and not tj in vsites:
            ff.write('    v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var + '(m_d_' + t + ', m_k_' + t + ', ' + set_m1[i] + ', ' + set_m1[j] + ');\n')
            nv = nv + 1
ff.write('\n')

a = """     
    
    const double E_poly = mb_system::poly_model::eval(m_poly, v);
    
    return E_poly;
}

double x1b_""" + mon1 + """_v1x::eval(const double* mon1, 
                double * grd1) const
{

    double xcrd[""" + str(3*(nat)) + """]; // coordinates of real sites ONLY

    std::copy(mon1, mon1 + """ + str(3*(nat)) + """, xcrd);
    
    double v[""" + str(nvars) + """];
    
"""
ff.write(a)
nc = 0
set_m1 = []
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if not types[i] in vsites:
            ff.write('    const double* ' + types[i] + '_' + str(n) + '= xcrd + ' + str(3 * nc) + ';\n')
            set_m1.append(types[i] + '_' + str(n))
            n = n + 1
            nc = nc + 1

ff.write('\n')            
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if types[i] in vsites:
            ff.write('    double ' + types_a[i] + '_' + str(n) + '[3]' + ';\n')
            set_m1.append(types[i] + '_' + str(n))
            n = n + 1
ff.write('\n')


 
a = """    
    variable vr[""" + str(nvars) + """];
    
"""
ff.write(a)
nv = 0
# Intramolecular distances:
for i in range(0,len(set_m1) - 1):
    for j in range(i + 1,len(set_m1)):
        ti = set_m1[i].split('_')[0]
        tj = set_m1[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        if not ti in vsites and not tj in vsites:
            ff.write('    v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var + '(m_d_' + t + ', m_k_' + t + ', ' + set_m1[i] + ', ' + set_m1[j] + ');\n')
            nv = nv + 1
ff.write('\n')
a = """     
    
    double g[""" + str(nvars) + """];
    
    const double E_poly = mb_system::poly_model::eval(m_poly, v, g);
    
    double xgrd[""" + str(3*(len(t1))) + """];
    std::fill(xgrd, xgrd + """ + str(3*(len(t1))) + """, 0.0);

""" 
ff.write(a)   
nc = 0

for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if not types[i] in vsites:
            ff.write('    double* ' + types[i] + '_' + str(n) + '_g' + ' = xgrd + ' + str(3 * nc) + ';\n')
            n = n + 1
            nc = nc + 1
ff.write('\n')
          
for i in range(0,len(types),2):
    n = 1
    for j in range(int(types[i+1])):
        if types[i] in vsites:
            ff.write('    double* ' + types[i] + '_' + str(n) + '_g' + ' = xgrd + ' + str(3 * nc) + ';\n')
            n = n + 1
            nc = nc + 1
ff.write('\n')
    
nv = 0
# Intramolecular distances:
for i in range(0,len(set_m1) - 1):
    for j in range(i + 1,len(set_m1)):
        ti = set_m1[i].split('_')[0]
        tj = set_m1[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        if not ti in vsites and not tj in vsites:
            ff.write('    vr[' + str(nv) + '].grads(g[' + str(nv) + '], ' + set_m1[i] + '_g, ' + set_m1[j] + '_g, ' + set_m1[i] + ', ' + set_m1[j] + ');\n')
            nv = nv + 1
ff.write('\n') 

a = """
    for (size_t i = 0; i < """ + str(3*nat) + """; i++)
        grd1[i] = xgrd[i];

    return E_poly;
}

double x1b_""" + mon1 + """_v1x::operator()(const double crd[""" + str(3*(nat)) + """]) const
{
    const double E_poly = eval(crd);

    return E_poly;
}

double x1b_""" + mon1 + """_v1x::operator()(const double crd[""" + str(3*(nat)) + """],
                            double grd[""" + str(3*(nat)) + """]) const
{
    const double E_poly = eval(crd, grd);

    return E_poly;
}

} // namespace x1b_""" + mon1 + """

////////////////////////////////////////////////////////////////////////////////
"""
ff.write(a)
ff.close()


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



