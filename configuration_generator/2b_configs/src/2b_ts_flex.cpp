/* Original code from Greg Medders. Adapted by Marc Riera in 20171005*/

#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "mt19937.h"
#include "random-rotation.h"
#include "io-xyz.h"
#include "io-utils.h"
#include "read-xyz.h"

namespace {
  const double step_size = 0.05;//Angstrom
  const double max_dist = 9.0;//Angstrom

  static std::vector<tools::xyz_frame> set1;
  static std::vector<tools::xyz_frame> set2;
}

////////////////////////////////////////////////////////////////////////////////

inline double dist(const double* r1, const double* r2)
{
    const double dx = r1[0] - r2[0];
    const double dy = r1[1] - r2[1];
    const double dz = r1[2] - r2[2];

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

////////////////////////////////////////////////////////////////////////////////

inline double getMass(std::string atom) {
  if (atom == "O") return 16.0;
  else if (atom == "H") return 1.0;
  else return 1.0;
}

////////////////////////////////////////////////////////////////////////////////

void com_to_origin( double* crd, std::vector<std::string> at_names) {
  // Center the molecule at the center of mass
  // Finc COM coordinate:
  std::vector<double> COM(3,0.0);
  double tot_mass = 0.0;
  for (size_t i = 0; i < at_names.size(); i++) {
    double mass = getMass(at_names[i]); 
    for (size_t j = 0; j < 3; j++) {
      COM[j] += crd[3*i + j]*mass;
    }
    tot_mass += mass;
  }

  for (size_t j = 0; j < 3; j++) {
    COM[j] = COM[j]/tot_mass;
  }

  // Center molecule
  for (size_t i = 0; i < at_names.size(); i++)
    for (size_t j = 0; j < 3; j++)
      crd[3*i + j] = crd[3*i + j] - COM[j];
}

////////////////////////////////////////////////////////////////////////////////

void rotate_atom(const double* rot, double* xyz) {
  double old[3] = {xyz[0], xyz[1], xyz[2]};

  for (size_t i = 0; i < 3; ++i) {
    xyz[i] = 0.0;
    for (size_t j = 0; j < 3; ++j)
      xyz[i] += rot[i*3 + j]*old[j];
  }
}

////////////////////////////////////////////////////////////////////////////////

void rotate_monomer(kit::mt19937& prg, double* xyz, size_t n) {
  double x[3];
  for (size_t k = 0; k < 3; ++k)
    x[k] = prg.random_double();

  double rot[9];
  random_rotation(x, rot);

  for (size_t k = 0; k < n; ++k)
    rotate_atom(rot, xyz + 3*k);
}

////////////////////////////////////////////////////////////////////////////////

void get_confs(const double R,
  std::vector<tools::xyz_frame> mon1, 
  std::vector<tools::xyz_frame> mon2, 
  size_t ns, kit::mt19937& prg, 
  const double d1, std::ofstream &ofs, int seed) {
  srand(seed);

  size_t count = 0;

  for (size_t mm = 0; mm < ns; mm++) { 
    size_t i1 = rand()%mon1.size();
    size_t i2 = rand()%mon2.size();

    size_t nat1 = mon1[i1].natm;
    size_t nat2 = mon2[i2].natm;

    std::vector<double> crd (3*nat1 + 3*nat2, 0.0);
    std::vector<std::string> ats (nat1 + nat2, "");

    double m1[3*nat1];
    double m2[3*nat2];

    double dim1[3*nat1 + 3*nat2];

    std::copy(mon1[i1].xyz, mon1[i1].xyz + 3*nat1, m1);
    std::copy(mon2[i2].xyz, mon2[i2].xyz + 3*nat2, m2);

    com_to_origin(m1, mon1[i1].at_name);
    com_to_origin(m2, mon2[i2].at_name);
    
    std::copy(mon1[i1].at_name.begin(), 
              mon1[i1].at_name.end(), ats.begin());
    std::copy(mon2[i2].at_name.begin(), 
              mon2[i2].at_name.end(), ats.begin() + nat1);

    std::copy(m1, m1 + 3*nat1, dim1);
    std::copy(m2, m2 + 3*nat2, dim1 + 3*nat1);

    rotate_monomer(prg, dim1, nat1);
    rotate_monomer(prg, dim1 + 3*nat1, nat2);

    // Mon1 stays at origin.
    // Mon2 is moved on the x axis by R +- 0.5 ang
    double s0 = prg.random_double() - 0.5;
    for (size_t j = 0; j < nat2; j++) {
      dim1[3*nat1 + 3*j] += R + s0;
    }

    bool good = true;
    for (size_t k = 0; k < nat1; k++) {
      for (size_t j = 0; j < nat2; j++) {
        if (dist(dim1 + 3*k, dim1 + 3*j + 3*nat1) < d1) {
          good = false;
          k = nat1;
          j = nat2;
        }
      }
    }
    
    if (good) {
      std::stringstream ss("0.0");
      std::copy(dim1, dim1 + 3*nat1 + 3*nat2, crd.begin());
      kit::io::save_xyz(ofs, ss.str(), ats, crd);
      count = 0;
    } else {
      mm--;
      count++;
    }
    
    if (count > 10*ns) break;
  }
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

  if (argc < 6) {
    std::cout << "Usage: " << argv[0] << std::endl
              << " <confs1.xyz> <confs2.xyz> " << std::endl
              << " <n_samples_x_distance> <min_d 12> <seed>"
              << std::endl;
    std::cout << "min_d ij is the minimum distance acceptable between any pair\n"
              << " of atoms between molecule i and j\n\n";

    return EXIT_FAILURE;
  }

  // 1. Load Configurations
  const double dR = step_size;
  size_t ns = 0;
  double d1 = 0.0;
  int seed = 0;

  try {
    // Checking each one of the input files in the command line
    tools::read_xyz(argv[1], set1);
    tools::read_xyz(argv[2], set2);
    ns = atoi(argv[3]);
    d1 = atof(argv[4]);
    seed = atoi(argv[5]);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }

  kit::mt19937 prg(seed);

  // 2. Do trapezoidal rule for numerical integration of radial coord
  //   Don't start at R = 0 because first (and final) points are zero
  std::ofstream ofs("configs.xyz", std::ios_base::app);
  for (double R = 0.6; R < max_dist; R += dR) {
    get_confs(R, set1, set2, ns, prg, d1, ofs, seed/2);

    std::cout << R << " done" << std::endl;
  }
  ofs.close();
  
  return 0;
}  
