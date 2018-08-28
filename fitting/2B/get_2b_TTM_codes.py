
# coding: utf-8

# In[ ]:


import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import settings_reader

# In[ ]:


# Check proper if input is provided


# In[ ]:


if len(sys.argv) != 4:
    print("Usage: ./script <settings.ini> <config.ini> <input.in>")
    sys.exit()
else:
    settings_path = sys.argv[1]
    config_path = sys.argv[2]
    name = sys.argv[3]

settings = settings_reader.SettingsReader(settings_path)
config = settings_reader.SettingsReader(config_path)

# In[ ]:


mon1 = name.split('.')[0].split('_')[0]
mon2 = name.split('.')[0].split('_')[1]


# In[ ]:


# For Andrea:
# Find a way to find the number of atoms in each monomer
# Store them in nat1 and nat2
nat1, nat2 = (int(num_atoms) for num_atoms in settings.get("molecule", "fragments").split(","))

# Find the number of sites
#nsites1 = 4
nsites1 = nat1
nsites2 = nat2

# Is one of the molecules water? 0 = none, 1 = mon1; 2=mon2
#is_w = 1
is_w = 0
# Obtain the lists with the excluded pairs
#excl12_a = [[0,1],[0,2],[0,3],[1,3],[2,3]]
excluded_pairs_12 = config.getlist("fitting", "excluded_pairs_12", int)
excluded_pairs_13 = config.getlist("fitting", "excluded_pairs_13", int)
excluded_pairs_14 = config.getlist("fitting", "excluded_pairs_14", int)
excl12_a = excluded_pairs_12[0]
excl13_a = excluded_pairs_13[0]
excl14_a = excluded_pairs_14[0]

excl12_b = excluded_pairs_12[1]
excl13_b = excluded_pairs_13[1]
excl14_b = excluded_pairs_14[1]

#Obtain charges (in the order of input), pols and polfacs
charges = config.getlist("fitting", "charges", float)
polarizabilities = config.getlist("fitting", "polarizabilities", float)
polarizability_fractions = config.getlist("fitting", "polarizability_fractions", float)
chg_a = charges[0]
chg_b = charges[1]
pol_a = polarizabilities[0]
pol_b = polarizabilities[1]
polfac_a = polarizability_fractions[0]
polfac_b = polarizability_fractions[1]

#Ask the user for the max value of k and d
k_min = config.getfloat("fitting", "k_min")
k_max = config.getfloat("fitting", "k_max")

# Obtain C6 from user in the same order as the given pairs AA, AB ...:
# ALL pairs are required (intra and inter)
c6_constants = config.getlist("fitting", "c6", float)
C6 = c6_constants[len(c6_constants) - 1]

# Define Energy Range for the fitting
E_range = config.getfloat("fitting", "energy_range")


# In[ ]:


types_a = list(mon1)
types_b = list(mon2)

#Concatenating excluded pairs:
excl_a = excl12_a + excl13_a + excl14_a
excl_b = excl12_b + excl13_b + excl14_b


# In[ ]:


# Generating the non linear parameter list
t1 = []
# Monomer 1 parameters
for i in range(0,len(types_a),2):
    for j in range(int(types_a[i+1])):
        t1.append(types_a[i])
t2 = []
# Monomer 2 parameters
for i in range(0,len(types_b),2):
    for j in range(int(types_b[i+1])):
        t2.append(types_b[i])


# In[ ]:


# Get non linear parameters and pairs involved
pairs = []
nlparam = []

# No intramolecular pairs. Those contributions are expected to be removed by interaction energy
# Intermolecular    
for i in range(len(t1)):
    for j in range(len(t2)):
        p = sorted([t1[i],t2[j]])
        ps = p[0] + p[1]
        const_inter = 'k_' + p[0] + p[1]
        if not const_inter in nlparam:
            nlparam.append(const_inter)
            pairs.append(ps)


# In[ ]:


print(pairs)
print(nlparam)


# In[ ]:


# Save number in num_nonlinear
num_nonlinear = len(nlparam)


# ## Creating mon1.h and mon2.h

# In[ ]:


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


# In[ ]:


mon1_class = open('mon2.h','w')

a="""
#ifndef MON2_H
#define MON2_H
#include "molecule.h"

namespace x   {

   class mon2 : public molecule {
      public :
      mon2();
      ~mon2();
      
      mon2(double* crd);
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

mon1_class.write(a)
mon1_class.close()


# ## Creating their cpp files

# In[ ]:


ff = open('mon1.cpp','w')

a = """
#include "mon1.h"
#include "ps.h"
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
    
    nsites = """ + str(nsites1) + """; 
    realsites = """ + str(nat1) + """;
    is_w = 0;
    allocate();
    
    sitecrds = set_sitecrds(crd);
    atmpolar = set_pol();
    charge = set_charges(crd);
    polfac = set_polfacs(atmpolar);

    // Excluded pairs
"""
ff.write(a)
for p in excl12_a:
    ff.write('    excluded12.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
for p in excl13_a:
    ff.write('    excluded13.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
for p in excl14_a:
    ff.write('    excluded14.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
a = """

  }

  double* mon1::set_charges(double* atmcrds) {
    charge = memory;
"""
ff.write(a)
for i in range(len(chg_a)):
    ff.write('    charge[' + str(i) + '] = ' + chg_a[i] + '*CHARGECON;\n')
a = """

    return charge;
  }
  
  double* mon1::set_pol() {
    atmpolar = memory + nsites + nsites*3;
"""
ff.write(a)
for i in range(len(pol_a)):
    ff.write('    atmpolar[' + str(i) + '] = ' + pol_a[i] + ';\n')
a = """
    return atmpolar;
  }

  double* mon1::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;
    
"""
ff.write(a)
for i in range(len(polfac_a)):
    ff.write('    polfac[' + str(i) + '] = ' + polfac_a[int(i)] + ';\n')
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


# In[ ]:


ff = open('mon2.cpp','w')

a = """
#include "mon2.h"
#include "ps.h"
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

  mon2::mon2() { }

  mon2::~mon2() {
    delete[] memory;
  }

  mon2::mon2( double* crd) {
    
    nsites = """ + str(nsites2) + """;  
    realsites = """ + str(nat2) + """;
    
    is_w = 0;
    
    allocate();
    
    sitecrds = set_sitecrds(crd);
    atmpolar = set_pol();
    charge = set_charges(crd);
    polfac = set_polfacs(atmpolar);

    // Excluded pairs
"""
ff.write(a)
for p in excl12_b:
    ff.write('    excluded12.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
for p in excl13_b:
    ff.write('    excluded13.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
for p in excl14_b:
    ff.write('    excluded14.insert(std::make_pair(' + str(p[0]) + ',' + str(p[1]) + '));\n')
a = """

  }

  double* mon2::set_charges(double* atmcrds) {
    charge = memory;
"""
ff.write(a)
for i in range(len(chg_b)):
    ff.write('    charge[' + str(i) + '] = ' + chg_b[i] + '*CHARGECON;\n')
a = """

    return charge;
  }
  
  double* mon2::set_pol() {
    atmpolar = memory + nsites + nsites*3;
"""
ff.write(a)
for i in range(len(pol_b)):
    ff.write('    atmpolar[' + str(i) + '] = ' + pol_b[i] + ';\n')
a = """
    return atmpolar;
  }

  double* mon2::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;
    
"""
ff.write(a)
for i in range(len(polfac_b)):
    ff.write('    polfac[' + str(i) + '] = ' + polfac_b[i] + ';\n')
a = """
    return polfac;
  }

  double* mon2::set_sitecrds(double* atmcrds) {
    sitecrds = memory + nsites;
    
// ##DEFINE HERE## ficticious site coordinates
    return atmcrds;  
  }

  void mon2::allocate() {
    memory = new double [nsites  //charge
      + nsites*3  //site coordinates
      + nsites  //polarizabilities 
      + nsites];  //polfac
    }

    int mon2::get_nsites() { return nsites; }
    int mon2::get_realsites() { return realsites; }
    double* mon2::get_charges() { return charge; }
    double* mon2::get_sitecrds() { return sitecrds; }
    double* mon2::get_pol() { return atmpolar; }
    double* mon2::get_polfacs() { return polfac; }

    excluded_set_type::iterator mon2::get_begin_12() { return excluded12.begin(); }
    excluded_set_type::iterator mon2::get_begin_13() { return excluded13.begin(); }
    excluded_set_type::iterator mon2::get_begin_14() { return excluded14.begin(); }
    excluded_set_type::iterator mon2::get_end_12() { return excluded12.end(); }
    excluded_set_type::iterator mon2::get_end_13() { return excluded13.end(); }
    excluded_set_type::iterator mon2::get_end_14() { return excluded14.end(); }

} // namespace x
"""

ff.write(a)
ff.close()


# ## Creating water monomer
# If applicable...

# In[ ]:


if is_w != 0:
    ff = open('mon' + str(is_w) + '.cpp','w')
    a = """
#include "mon""" + str(is_w) + """.h"

#include <iostream>

//TODO: merge constants.h and fit-constants.h appropriately
#include "fit-constants.h"

namespace {

const double gammaM = 0.426706882;
const double gamma1 = 1.0 - gammaM;
const double gamma2 = gammaM/2;

const double CHARGECON = constants::CHARGECON;

inline void compute_M_site_crd
    (const double O[3], const double H1[3], const double H2[3], double M[3])
{
    for (size_t i = 0; i < 3; ++i)
        M[i] = gamma1*O[i] + gamma2*(H1[i] + H2[i]);
}

} // namespace

namespace x {

mon""" + str(is_w) + """::mon""" + str(is_w) + """() {

}

mon""" + str(is_w) + """::~mon""" + str(is_w) + """() {
  delete[] memory;
}

mon""" + str(is_w) + """::mon""" + str(is_w) + """(double* crd) {
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

double* mon""" + str(is_w) + """::set_sitecrds(double* atmcrds) {
    sitecrds = memory + nsites;
    // assumes O H H 
    compute_M_site_crd(atmcrds, atmcrds + 3, atmcrds + 6, sitecrds + 9);
    std::copy(atmcrds, atmcrds + 9, sitecrds);
    return sitecrds;
}

double* mon""" + str(is_w) + """::set_charges(double* atmcrds) {
    double chgtmp[3];
    h2o::ps::dms_nasa(0.0, 0.0, 0.0, atmcrds, chgtmp, 0, false);
    const double tmp = 0.5*gammaM/(1.0 - gammaM);
    charge = memory;

    charge[0] = 0.0;                        // O
    charge[1] = CHARGECON*(chgtmp[1] + tmp*(chgtmp[1] + chgtmp[2])); // H1
    charge[2] = CHARGECON*(chgtmp[2] + tmp*(chgtmp[1] + chgtmp[2])); // H2
    charge[3] = CHARGECON*(chgtmp[0]/(1.0 - gammaM));       // M

    //std::cout << "charge[1] = " << charge[1] << std::endl;

    return charge;
}

double* mon""" + str(is_w) + """::set_pol() {
    atmpolar = memory + nsites + nsites*3;
    atmpolar[0] = 1.310; // polarO
    atmpolar[1] = 0.294; // polarH
    atmpolar[2] = 0.294; // polarH
    atmpolar[3] = 0.0;   // polarM
    return atmpolar;
}

double* mon""" + str(is_w) + """::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;
    polfac[0] = 1.310; // polarO
    polfac[1] = 0.294; // polarH
    polfac[2] = 0.294; // polarH
    polfac[3] = 1.31;   // polarM
    return polfac;
}

void mon""" + str(is_w) + """::allocate() {
    memory = new double [nsites // charges
  + nsites*3              // sitecrds 
  + nsites                // polarizabilities
  + nsites];              // polfacs
}

int mon""" + str(is_w) + """::get_nsites() { return nsites; }
int mon""" + str(is_w) + """::get_realsites() { return realsites; }
double* mon""" + str(is_w) + """::get_sitecrds() { return sitecrds; }
double* mon""" + str(is_w) + """::get_charges() { return charge; }
double* mon""" + str(is_w) + """::get_polfacs() { return polfac; }
double* mon""" + str(is_w) + """::get_pol() { return atmpolar; }

excluded_set_type::iterator mon""" + str(is_w) + """::get_begin_12() { return excluded12.begin(); }
excluded_set_type::iterator mon""" + str(is_w) + """::get_begin_13() { return excluded13.begin(); }
excluded_set_type::iterator mon""" + str(is_w) + """::get_begin_14() { return excluded14.begin(); }
excluded_set_type::iterator mon""" + str(is_w) + """::get_end_12() { return excluded12.end(); }
excluded_set_type::iterator mon""" + str(is_w) + """::get_end_13() { return excluded13.end(); }
excluded_set_type::iterator mon""" + str(is_w) + """::get_end_14() { return excluded14.end(); }


} // namespace x

"""
    ff.write(a)
    ff.close()


# ## Create training_set.h/cpp files

# In[ ]:


ff = open('training_set.h','w')
a = """
#ifndef TRAINING_SET_H
#define TRAINING_SET_H

#include <cstdlib>
#include <vector>

namespace tset {

struct dimer {
    double energy_total;
    double energy_twobody;
    double energy_onebody[2];
    
    // Note: XYZ is for the real sites
"""
ff.write(a)
ff.writelines("    double xyz[" + str((nat1 + nat2)*3) + "];\n")
a = """
    double mb_total_energy() const
    {
        return energy_total;
    }
};

size_t load_dimers(const char* filename, std::vector<dimer>& ts);

} // namespace tset

#endif // TRAINING_SET_H
"""
ff.write(a)
ff.close()


# ## Buckingham.h

# In[ ]:


hname = "buckingham.h"
ff = open(hname,'w')
a = """
#ifndef BUCKINGHAM_H
#define BUCKINGHAM_H

#include <cmath>
#include <algorithm>



struct x2b_buck {
  x2b_buck();
  x2b_buck(double *, double * , size_t, size_t);
  ~x2b_buck();
  
  void set_parameters(double * a, double * b);
  void get_basis(double * mmm);
  
"""
ff.write(a)
for i in range(len(pairs)):
    ff.write('  double m_A_' + pairs[i] + '; \n')
for i in range(len(pairs)):
    ff.write('  double m_b_' + pairs[i] + '; \n')
    
a = """

  double get_buckingham();
  double get_buckingham(double * grdx);
  
  double * xyz1;
  double * xyz2;

  inline double buck(const double a, const double b,
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

  inline double buck(const double a, const double b,
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


# ## Buckingham.cpp

# In[ ]:


cppname = "buckingham.cpp"
ff = open(cppname,'w')
a = """
#include "buckingham.h"

x2b_buck::x2b_buck() {
  xyz1 = new double[3];
  xyz2 = new double[3];
}
x2b_buck::~x2b_buck() {
  delete[] xyz1;
  delete[] xyz2;
}

x2b_buck::x2b_buck(double * c1, double * c2, size_t n1, size_t n2) {
  xyz1 = new double[3*n1];
  xyz2 = new double[3*n2];
  std::copy(c1, c1 + 3*n1, xyz1);
  std::copy(c2, c2 + 3*n2, xyz2);
}

void x2b_buck::set_parameters(double * a, double * b) {
"""
ff.write(a)
for i in range(len(pairs)):
    ff.write('  m_A_' + pairs[i] + ' = a[' + str(i) + ']; \n')
for i in range(len(pairs)):
    ff.write('  m_b_' + pairs[i] + ' = b[' + str(i) + ']; \n')
    
a = """

}

double x2b_buck::get_buckingham() {

  double ebuck = 0.0;
"""
ff.write(a)

nc = 0
set_m1 = []
set_m2 = []
for i in range(0,len(types_a),2):
    n = 1
    for j in range(int(types_a[i+1])):
        ff.write('  const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xyz1 + ' + str(3 * nc) + ';\n')
        set_m1.append(types_a[i] + '_' + str(n) + '_a')
        n = n + 1
        nc = nc + 1
ff.write('\n')

nc = 0
for i in range(0,len(types_b),2):
    n = 1
    for j in range(int(types_b[i+1])):
        ff.write('  const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xyz2 + ' + str(3 * nc) + ';\n')
        set_m2.append(types_b[i] + '_' + str(n) + '_b')
        n = n + 1
        nc = nc + 1
ff.write('\n')

for i in range(0,len(set_m1)):
    for j in range(0,len(set_m2)):
        ti = set_m1[i].split('_')[0]
        tj = set_m2[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        ff.write('  ebuck += buck(m_A_' + t + ', m_b_' + t + ', ' + set_m1[i] + ', ' + set_m2[j] + ');\n')
        
    ff.write('\n')
a = """
  return ebuck;
}

void x2b_buck::get_basis(double * mmm) {

  for (int i = 0; i < """ + str(len(pairs)) + """; i++) {
    mmm[i] = 0.0;
  }
  
"""
ff.write(a)

nc = 0
set_m1 = []
set_m2 = []
for i in range(0,len(types_a),2):
    n = 1
    for j in range(int(types_a[i+1])):
        ff.write('  const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xyz1 + ' + str(3 * nc) + ';\n')
        set_m1.append(types_a[i] + '_' + str(n) + '_a')
        n = n + 1
        nc = nc + 1
ff.write('\n')

nc = 0
for i in range(0,len(types_b),2):
    n = 1
    for j in range(int(types_b[i+1])):
        ff.write('  const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xyz2 + ' + str(3 * nc) + ';\n')
        set_m2.append(types_b[i] + '_' + str(n) + '_b')
        n = n + 1
        nc = nc + 1
ff.write('\n')

for i in range(0,len(set_m1)):
    for j in range(0,len(set_m2)):
        ti = set_m1[i].split('_')[0]
        tj = set_m2[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        index_t = pairs.index(t)
        ff.write('  mmm[' + str(index_t) + '] += buck(1.0, m_b_' + t + ', ' + set_m1[i] + ', ' + set_m2[j] + ');\n')
        
    ff.write('\n')
a = """
  
}

double x2b_buck::get_buckingham(double * grd) {

  double ebuck = 0.0;
"""
ff.write(a)

nc = 0
set_m1 = []
set_m2 = []
for i in range(0,len(types_a),2):
    n = 1
    for j in range(int(types_a[i+1])):
        ff.write('  const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xyz1 + ' + str(3 * nc) + ';\n')
        set_m1.append(types_a[i] + '_' + str(n) + '_a')
        n = n + 1
        nc = nc + 1
ff.write('\n')

nc = 0
for i in range(0,len(types_b),2):
    n = 1
    for j in range(int(types_b[i+1])):
        ff.write('  const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xyz2 + ' + str(3 * nc) + ';\n')
        set_m2.append(types_b[i] + '_' + str(n) + '_b')
        n = n + 1
        nc = nc + 1
ff.write('\n')

nc = 0
set_m1 = []
set_m2 = []
for i in range(0,len(types_a),2):
    n = 1
    for j in range(int(types_a[i+1])):
        ff.write('  double* ' + types_a[i] + '_' + str(n) + '_a_g' + '= grd + ' + str(3 * nc) + ';\n')
        set_m1.append(types_a[i] + '_' + str(n) + '_a')
        n = n + 1
        nc = nc + 1
ff.write('\n')

for i in range(0,len(types_b),2):
    n = 1
    for j in range(int(types_b[i+1])):
        ff.write('  double* ' + types_b[i] + '_' + str(n) + '_b_g' + '= grd + ' + str(3 * nc) + ';\n')
        set_m2.append(types_b[i] + '_' + str(n) + '_b')
        n = n + 1
        nc = nc + 1
ff.write('\n')

for i in range(0,len(set_m1)):
    for j in range(0,len(set_m2)):
        ti = set_m1[i].split('_')[0]
        tj = set_m2[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        ff.write('  ebuck += buck(m_A_' + t + ', m_b_' + t + ',  \n             '                  + set_m1[i] + ', ' + set_m2[j] + ', ' + set_m1[i] + '_g, ' + set_m2[j] +  '_g);\n')
        
    ff.write('\n')
a = """
  return ebuck;
}
"""
ff.write(a)
ff.close()


# ## Dispersion.h

# In[ ]:


hname = "dispersion.h"
ff = open(hname,'w')
a = """
#ifndef DISPERSION_H
#define DISPERSION_H

#include "tang-toennies.h"
#include <cmath>
#include <algorithm>



struct x2b_disp {
  x2b_disp();
  x2b_disp(double *, double * , size_t, size_t);
  ~x2b_disp();
  
  void set_parameters(double * b);

"""
ff.write(a)
for i in range(len(pairs)):
    ff.write('  const double m_C6_' + pairs[i] + ' = ' + C6[i] + ' ; \n')
for i in range(len(pairs)):
    ff.write('  double m_d6_' + pairs[i] + '; \n')
    
a = """

  const double m_C8 = 0.0;
  const double m_d8 = 0.0;

  const double if6 = 1.0/x2o::factorial<6>();
  const double if8 = 1.0/x2o::factorial<8>();

  double * xyz1;
  double * xyz2;

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

# In[ ]:


cppname = "dispersion.cpp"
ff = open(cppname,'w')
a = """
#include "dispersion.h"

x2b_disp::x2b_disp() {
  xyz1 = new double[3];
  xyz2 = new double[3];
}
x2b_disp::~x2b_disp() {
  delete[] xyz1;
  delete[] xyz2;
}

x2b_disp::x2b_disp(double * c1, double * c2, size_t n1, size_t n2) {
  xyz1 = new double[3*n1];
  xyz2 = new double[3*n2];
  std::copy(c1, c1 + 3*n1, xyz1);
  std::copy(c2, c2 + 3*n2, xyz2);
}

void x2b_disp::set_parameters(double * b) {
"""
ff.write(a)
for i in range(len(pairs)):
    ff.write('  m_d6_' + pairs[i] + ' = b[' + str(i) + ']; \n')
    
a = """
}

double x2b_disp::get_dispersion() {

  double disp = 0.0;
"""
ff.write(a)

nc = 0
set_m1 = []
set_m2 = []
for i in range(0,len(types_a),2):
    n = 1
    for j in range(int(types_a[i+1])):
        ff.write('  const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xyz1 + ' + str(3 * nc) + ';\n')
        set_m1.append(types_a[i] + '_' + str(n) + '_a')
        n = n + 1
        nc = nc + 1
ff.write('\n')

nc = 0
for i in range(0,len(types_b),2):
    n = 1
    for j in range(int(types_b[i+1])):
        ff.write('  const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xyz2 + ' + str(3 * nc) + ';\n')
        set_m2.append(types_b[i] + '_' + str(n) + '_b')
        n = n + 1
        nc = nc + 1
ff.write('\n')

for i in range(0,len(set_m1)):
    for j in range(0,len(set_m2)):
        ti = set_m1[i].split('_')[0]
        tj = set_m2[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        ff.write('  disp += x6(m_C6_' + t + ', m_d6_' + t + ', m_C8, m_d8, ' + set_m1[i] + ', ' + set_m2[j] + ');\n')
        
    ff.write('\n')
a = """
  return disp;
}

double x2b_disp::get_dispersion(double * grd) {

  double disp = 0.0;
"""
ff.write(a)

nc = 0
set_m1 = []
set_m2 = []
for i in range(0,len(types_a),2):
    n = 1
    for j in range(int(types_a[i+1])):
        ff.write('  const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xyz1 + ' + str(3 * nc) + ';\n')
        set_m1.append(types_a[i] + '_' + str(n) + '_a')
        n = n + 1
        nc = nc + 1
ff.write('\n')

nc = 0
for i in range(0,len(types_b),2):
    n = 1
    for j in range(int(types_b[i+1])):
        ff.write('  const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xyz2 + ' + str(3 * nc) + ';\n')
        set_m2.append(types_b[i] + '_' + str(n) + '_b')
        n = n + 1
        nc = nc + 1
ff.write('\n')

nc = 0
set_m1 = []
set_m2 = []
for i in range(0,len(types_a),2):
    n = 1
    for j in range(int(types_a[i+1])):
        ff.write('  double* ' + types_a[i] + '_' + str(n) + '_a_g' + '= grd + ' + str(3 * nc) + ';\n')
        set_m1.append(types_a[i] + '_' + str(n) + '_a')
        n = n + 1
        nc = nc + 1
ff.write('\n')

for i in range(0,len(types_b),2):
    n = 1
    for j in range(int(types_b[i+1])):
        ff.write('  double* ' + types_b[i] + '_' + str(n) + '_b_g' + '= grd + ' + str(3 * nc) + ';\n')
        set_m2.append(types_b[i] + '_' + str(n) + '_b')
        n = n + 1
        nc = nc + 1
ff.write('\n')

for i in range(0,len(set_m1)):
    for j in range(0,len(set_m2)):
        ti = set_m1[i].split('_')[0]
        tj = set_m2[j].split('_')[0]
        t = ''.join(sorted(ti + tj))
        ff.write('  disp += x6(m_C6_' + t + ', m_d6_' + t + ', m_C8, m_d8, \n             '                  + set_m1[i] + ', ' + set_m2[j] + ', ' + set_m1[i] + '_g, ' + set_m2[j] +  '_g);\n')
        
    ff.write('\n')
a = """
  return disp;
}
"""
ff.write(a)
ff.close()


# ## Fitting routine

# In[ ]:


cppname = "fit-2b-ttm.cpp"
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

#include "wlsq.h"
#include "mon1.h"
#include "mon2.h"
#include "fit-utils.h"
#include "training_set.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"
#include "buckingham.h"

//#define VERBOSE

using namespace std;

namespace {

const double E_range = """ + E_range + """; // kcal/mol
const int nparams = """ + str(len(nlparam)) + """; 

//----------------------------------------------------------------------------//

static std::vector<tset::dimer> training_set;
static std::vector<double> elec_e;
static std::vector<double> disp_e;
static std::vector<double> rep_e;
static double* ts_weights = 0;

namespace linear {

//----------------------------------------------------------------------------//

static double* A(0);
static double* y(0);
static double* params(0);

//----------------------------------------------------------------------------//

void allocate()
{
    A = new double[training_set.size()*nparams];
    y = new double[training_set.size()];

    params = new double[nparams];
}

bool nonlinear_parameters_out_of_range(double * nlpar) {
    return false """

ff.write(a)
for p in range(len(nlparam)):
    ff.write("\n              || nlpar[" + str(p) + "] < " + k_min + " || nlpar[" + str(p) + "] > " + k_max)

a=""";

}

//----------------------------------------------------------------------------//

double compute_chisq(const gsl_vector* X, void* unused)
{

    if (nonlinear_parameters_out_of_range(X->data))
        return 1.0e+6;

    for (size_t n = 0; n < training_set.size(); ++n) {
        double mmm[nparams];
        y[n] = training_set[n].energy_twobody;
        
        x::mon1 m1(training_set[n].xyz);
        size_t nat1 = m1.get_realsites();
        
        x::mon2 m2(training_set[n].xyz + 3*m1.get_realsites());
        size_t nat2 = m2.get_realsites();
        
        x2b_buck buck(training_set[n].xyz, training_set[n].xyz + 3*nat1, nat1, nat2);
        buck.set_parameters(X->data,X->data);
        buck.get_basis(mmm);
        for (size_t p = 0; p < nparams; ++p) {
            A[p + n*nparams] = mmm[p];
        }
        
        x2b_disp disp(training_set[n].xyz, training_set[n].xyz + 3*nat1, nat1, nat2);
        disp.set_parameters(X->data);
        y[n] -= disp.get_dispersion();
    }

#   ifdef VERBOSE
    std::cout << "=== calling wlsq::solve() ["
              << kit::wlsq::implementation()
              << "] ===" << std::endl;
#   endif

    double chisq;


    int rank;
    kit::wlsq::solve(training_set.size(), nparams,
                   A, y, ts_weights, params, chisq, rank);
    std::cout << "<#> chisq = " << chisq
            << " : rank = " << rank
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
        std::cerr << "usage: fit-2b ts1 ..."
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
    ff.write('      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(k_max) - float(k_min)) + ' + ' + k_min + ';\n')

            
a = """

    std::cout << "\\n<><><> model type = TTM-nrg \\n";

    ++argv;

    try {
        while (--argc != 0) {
            size_t nd = tset::load_dimers(*argv, training_set);
            std::cout << "'" << *(argv++) << "' : "
                      << nd << " dimers" << std::endl;
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
      // Saving reference 2b energies  
      tb_ref[n] = training_set[n].energy_twobody ;

      x::mon1 m1(training_set[n].xyz);
      x::mon2 m2(training_set[n].xyz + 3*m1.get_realsites());

      int system_nsites = m1.get_nsites() + m2.get_nsites();
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
      std::fill(system_is_w + m1.get_nsites(), system_is_w + system_nsites, m2.is_w);
      
      int * is_w_a = system_is_w;
      int * is_w_b = system_is_w + m1.get_nsites();

      std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
                system_sitecrds);
      std::copy(m2.get_sitecrds(), m2.get_sitecrds() + 3 * m2.get_nsites(),
                system_sitecrds + 3 * m1.get_nsites());

      std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
                system_charge);
      std::copy(m2.get_charges(), m2.get_charges() + m2.get_nsites(),
                system_charge + m1.get_nsites());

      std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
                system_pol);
      std::copy(m2.get_pol(), m2.get_pol() + m2.get_nsites(),
                system_pol + m1.get_nsites());

      std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
                system_polfac);
      std::copy(m2.get_polfacs(), m2.get_polfacs() + m2.get_nsites(),
                system_polfac + m1.get_nsites());


      excluded_set_type exclude12;
      excluded_set_type exclude13;
      excluded_set_type exclude14;
      
      excluded_set_type exclude12_a;
      excluded_set_type exclude13_a;
      excluded_set_type exclude14_a;
      
      excluded_set_type exclude12_b;
      excluded_set_type exclude13_b;
      excluded_set_type exclude14_b;

      for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
        exclude12.insert(*i);
        exclude12_a.insert(*i);
      }
      for (auto i = m2.get_begin_12(); i != m2.get_end_12(); i++) {
        std::pair<size_t,size_t> p =
                      std::make_pair(i->first + m1.get_nsites() ,
                                     i->second + m1.get_nsites());
        exclude12.insert(p);
        exclude12_b.insert(*i);
      }

      for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
        exclude13.insert(*i);
        exclude13_a.insert(*i);
      }
      for (auto i = m2.get_begin_13(); i != m2.get_end_13(); i++) {
        std::pair<size_t,size_t> p =
                      std::make_pair(i->first + m1.get_nsites() ,
                                     i->second + m1.get_nsites());
        exclude13.insert(p);
        exclude13_b.insert(*i);
      }

      for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
        exclude14.insert(*i);
        exclude14_a.insert(*i);
      }
      for (auto i = m2.get_begin_14(); i != m2.get_end_14(); i++) {
        std::pair<size_t,size_t> p =
                      std::make_pair(i->first + m1.get_nsites() ,
                                     i->second + m1.get_nsites());
        exclude14.insert(p);
        exclude14_b.insert(*i);
      }

      ttm::electrostatics m_electrostatics;

      ttm::smear_ttm4x smr; 
      smr.m_aDD_intra_12 = 0.3;
      smr.m_aDD_intra_13 = 0.3;
      smr.m_aDD_intra_14 = 0.055;

      double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
                                    system_sitecrds, exclude12, exclude13, exclude14, 
                                    system_is_w, smr, 0);
      double ener_a = m_electrostatics(m1.get_nsites(), m1.get_charges(), m1.get_polfacs(), m1.get_pol(),
                                    m1.get_sitecrds(), exclude12_a, exclude13_a, exclude14_a, 
                                    is_w_a, smr, 0);
      double ener_b = m_electrostatics(m2.get_nsites(), m2.get_charges(), m2.get_polfacs(), m2.get_pol(),
                                    m2.get_sitecrds(), exclude12_b, exclude13_b, exclude14_b, 
                                    is_w_b, smr, 0);

      // Take out electrostatic energy:
      elec_e.push_back(ener - ener_a - ener_b);
      training_set[n].energy_twobody -= elec_e[n];

      
      delete[] system_sitecrds ;
      delete[] system_charge  ;
      delete[] system_polfac  ;
      delete[] system_pol ;
      delete[] system_is_w;

    }

        linear::allocate();

    gsl_vector* x = gsl_vector_alloc(nparams);
    std::copy(x0, x0 + nparams,x->data);

    gsl_multimin_function chisq_func;

    chisq_func.n = nparams;
    chisq_func.f = linear::compute_chisq;

    gsl_vector* ss = gsl_vector_alloc(nparams);
    gsl_vector_set_all(ss, 0.1);

    std::cout << "\\n<> initial simplex sides:\\n";
    for (size_t n = 0; n < nparams; ++n)
        std::cout << n << " : " << ss->data[n] << "\\n";

    gsl_multimin_fminimizer* s =
        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2rand,
                                      nparams);

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
            for (size_t n = 0; n <  nparams; ++n)
                std::cout << n << " : " << s->x->data[n] << "\\n";
            std::cout << "<>" << std::endl;
        }
    } while (status == GSL_CONTINUE && iter < 5000);

    double nl_sol[nparams];
    std::copy(s->x->data, s->x->data + nparams, nl_sol);
    linear::compute_chisq(s->x, 0);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    //
    // report
    //

    for (size_t i = 0; i < training_set.size(); ++i) {
      training_set[i].energy_twobody += elec_e[i];
      x::mon1 m1(training_set[i].xyz);
      size_t nat1 = m1.get_realsites();
        
      x::mon2 m2(training_set[i].xyz + 3*m1.get_realsites());
      size_t nat2 = m2.get_realsites();
        
      x2b_buck buck(training_set[i].xyz, training_set[i].xyz + 3*nat1, nat1, nat2);
      buck.set_parameters(linear::params,nl_sol);
      rep_e.push_back(buck.get_buckingham());
        
      x2b_disp disp(training_set[i].xyz, training_set[i].xyz + 3*nat1, nat1, nat2);
      disp.set_parameters(nl_sol);
      disp_e.push_back(disp.get_dispersion());
    }
      

    ofstream correlation_file;
    ofstream indiv_terms;
    correlation_file.open ("correlation.dat");
    indiv_terms.open("individual_terms.dat");
    
    indiv_terms << "E_repulsion     E_dispersion    E_electrostatics   E_total         E_ref\\n";
    double err_L2(0), err_wL2(0), err_Linf(0);
    double err_L2_lo(0), err_Linf_lo(0), nlo(0);

    for (size_t i = 0; i < training_set.size(); ++i) {

        const double E_model = rep_e[i] + elec_e[i] + disp_e[i];
        const double delta = E_model - tb_ref[i];
        if (std::abs(delta) > err_Linf)
            err_Linf = std::abs(delta);

        correlation_file <<  i+1   << "   "
                         <<  E_model  << "   "
                         <<  tb_ref[i]  << "    "
                         <<  delta*delta  << "    \\n" ;

        err_L2 += delta*delta;
        err_wL2 += ts_weights[i]*delta*delta;

        if (training_set[i].energy_total - E_min < E_range) {
            nlo += 1.0;
            err_L2_lo += delta*delta;
            if (std::abs(delta) > err_Linf_lo)
                err_Linf_lo = std::abs(delta);
        }
        
        indiv_terms << rep_e[i] << "  " << disp_e[i] <<  "  " << elec_e[i] << "  " << E_model << "  " << tb_ref[i] << "  " << std::endl;
    }

    correlation_file.close();
    indiv_terms.close();

    err_L2 /= training_set.size();
    err_wL2 /= training_set.size();

    err_L2_lo /= nlo;

    

    //
    //  save
    //

    ofstream result;
    result.open("ttm-params.txt");
    for (size_t i = 0; i < nparams; i++) {
      result << linear::params[i] << " ";
    }
    result << std::endl;
    
    for (size_t i = 0; i < nparams; i++) {
      result << nl_sol[i] << " ";
    }
    result << std::endl;
    
    result.close();

    std::cout << "      err[L2] = " << std::sqrt(err_L2) << "    #rmsd of full ts\\n"
              << "     err[wL2] = " << std::sqrt(err_wL2) << "   #weighted rmsd of full ts\\n"
              << "    err[Linf] = " << err_Linf << "   #highest error in full ts\\n"
              << "  err[L2,low] = " << std::sqrt(err_L2_lo) << "   #rmsd of low-energy ts \\n"
              << "err[Linf,low] = " << err_Linf_lo << "   #highest error in low-energy ts "
              << std::endl;

    return 0;
}

"""
ff.write(a)
ff.close()


# ## Makefile

# In[ ]:


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

FIT_OBJ = fit-utils.o coulomb.o electrostatics.o gammq.o io-xyz.o \\
kvstring.o mon1.o mon2.o ps.o wlsq.o stuff.o tang-toennies.o \\
training_set.o ttm4.o dispersion.o buckingham.o 

EVAL_OBJ = fit-utils.o coulomb.o electrostatics.o gammq.o io-xyz.o \\
kvstring.o mon1.o mon2.o ps.o wlsq.o stuff.o tang-toennies.o \\
training_set.o ttm4.o dispersion.o buckingham.o


all: libfit.a libeval.a fit-2b-ttm eval-2b-ttm

libfit.a: $(addprefix $(OBJDIR)/, $(FIT_OBJ))
\t$(AR) cru libfit.a $(addprefix $(OBJDIR)/, $(FIT_OBJ))

libeval.a: $(addprefix $(OBJDIR)/, $(EVAL_OBJ))
\t$(AR) cru libeval.a $(addprefix $(OBJDIR)/, $(EVAL_OBJ))

fit-2b-ttm: fit-2b-ttm.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) $< $(LIBS) -lfit -o $@

eval-2b-ttm: eval-2b-ttm.cpp
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) $< $(LIBS) -leval -o $@

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/.sentinel
\t$(CXX) $(CXXFLAGS) $(INCLUDE) -L$(LIBDIR) -c $< $(LIBS) -o $@

$(OBJDIR)/.sentinel:
\tmkdir -p $(OBJDIR)
\ttouch $@

clean:
\trm -rf $(addprefix $(OBJDIR)/, $(FIT_OBJ)) libfit*.a libeval*.a fit-2b-ttm eval-2b-ttm
"""
ff.write(a)
ff.close()


# ## Evaluation code

# In[ ]:


ff = open('eval-2b-ttm.cpp','w')
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
#include "mon2.h"
#include "training_set.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"
#include "buckingham.h"
#include "io-xyz.h"

#define GRADIENTS

static std::vector<double> elec_e;
static std::vector<double> disp_e;
static std::vector<double> buck_e;
static std::vector<double> aparams(""" + str(len(nlparam)) + """,0.0);
static std::vector<double> bparams(""" + str(len(nlparam)) + """,0.0);

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: eval ttm-params.txt dimer.xyz"
                  << std::endl;
        return 0;
    }
    std::cout << std::scientific << std::setprecision(9);
    std::vector<std::string> elements;
    std::vector<double> crd;

    try {
        ++argv;
        --argc;
        
        std::ifstream parameters;
        parameters.open(*argv);
        
        for (size_t i = 0; i < """ + str(len(nlparam)) + """; i++)
          parameters >> aparams[i];
          
        for (size_t i = 0; i < """ + str(len(nlparam)) + """; i++)
          parameters >> bparams[i];

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
    
    double xyz[""" + str(3*(nat1 + nat2)) + """];
    std::copy(crd.begin(), crd.end(), xyz);
    
    x::mon1 m1(xyz);
    x::mon2 m2(xyz + 3*m1.get_realsites());
    
    int system_nsites = m1.get_nsites() + m2.get_nsites();
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
    std::fill(system_is_w + m1.get_nsites(), system_is_w + system_nsites, m2.is_w);
    
    int * is_w_a = system_is_w;
    int * is_w_b = system_is_w + m1.get_nsites();
    
    std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
              system_sitecrds);
    std::copy(m2.get_sitecrds(), m2.get_sitecrds() + 3 * m2.get_nsites(),
              system_sitecrds + 3 * m1.get_nsites());

    std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
              system_charge);
    std::copy(m2.get_charges(), m2.get_charges() + m2.get_nsites(),
              system_charge + m1.get_nsites());

    std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
              system_pol);
    std::copy(m2.get_pol(), m2.get_pol() + m2.get_nsites(),
              system_pol + m1.get_nsites());

    std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
              system_polfac);
    std::copy(m2.get_polfacs(), m2.get_polfacs() + m2.get_nsites(),
              system_polfac + m1.get_nsites());

    excluded_set_type exclude12;
    excluded_set_type exclude13;
    excluded_set_type exclude14;

    excluded_set_type exclude12_a;
    excluded_set_type exclude13_a;
    excluded_set_type exclude14_a;
      
    excluded_set_type exclude12_b;
    excluded_set_type exclude13_b;
    excluded_set_type exclude14_b;

    for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
      exclude12.insert(*i);
      exclude12_a.insert(*i);
    }
    for (auto i = m2.get_begin_12(); i != m2.get_end_12(); i++) {
      std::pair<size_t,size_t> p =
                    std::make_pair(i->first + m1.get_nsites() ,
                                   i->second + m1.get_nsites());
      exclude12.insert(p);
      exclude12_b.insert(*i);
    }

    for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
      exclude13.insert(*i);
      exclude13_a.insert(*i);
    }
    for (auto i = m2.get_begin_13(); i != m2.get_end_13(); i++) {
      std::pair<size_t,size_t> p =
                    std::make_pair(i->first + m1.get_nsites() ,
                                   i->second + m1.get_nsites());
      exclude13.insert(p);
      exclude13_b.insert(*i);
    }

    for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
      exclude14.insert(*i);
      exclude14_a.insert(*i);
    }
    for (auto i = m2.get_begin_14(); i != m2.get_end_14(); i++) {
      std::pair<size_t,size_t> p =
                    std::make_pair(i->first + m1.get_nsites() ,
                                   i->second + m1.get_nsites());
      exclude14.insert(p);
      exclude14_b.insert(*i);
    }

    ttm::electrostatics m_electrostatics;

    ttm::smear_ttm4x smr; 
    smr.m_aDD_intra_12 = 0.3;
    smr.m_aDD_intra_13 = 0.3;
    smr.m_aDD_intra_14 = 0.055;

    double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
                                  system_sitecrds, exclude12, exclude13, exclude14, 
                                  system_is_w, smr, 0);
    double ener_a = m_electrostatics(m1.get_nsites(), m1.get_charges(), m1.get_polfacs(), m1.get_pol(),
                                  m1.get_sitecrds(), exclude12_a, exclude13_a, exclude14_a, 
                                  is_w_a, smr, 0);
    double ener_b = m_electrostatics(m2.get_nsites(), m2.get_charges(), m2.get_polfacs(), m2.get_pol(),
                                  m2.get_sitecrds(), exclude12_b, exclude13_b, exclude14_b, 
                                  is_w_b, smr, 0);

    // Take out electrostatic energy:
    elec_e.push_back(ener - ener_a - ener_b);
      
    // Now need to take out dispersion
    x2b_disp disp(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
    disp.set_parameters(bparams.data());
    ener = disp.get_dispersion();
    disp_e.push_back(ener);
    
    // Now need to take out buckingham
    x2b_buck buck(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
    buck.set_parameters(aparams.data(),bparams.data());
    ener = buck.get_buckingham();
    buck_e.push_back(ener);
    
    std::cout << "2B_nograd = " << buck_e[0] + elec_e[0] + disp_e[0] << std::endl;
    std::cout << "E_elec2b = " << elec_e[0] << std::endl;
    std::cout << "E_disp2b = " << disp_e[0] << std::endl;
    std::cout << "E_buck2b = " << buck_e[0] << std::endl;

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

