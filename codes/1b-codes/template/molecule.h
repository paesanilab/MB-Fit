//This contains the base structure "molecule". 
//The derived structures will the specific molecules.

#ifndef MOLECULE_H
#define MOLECULE_H

#include <iterator>
#include <set>
#include <utility>
#include <cstddef>

typedef std::set<std::pair<size_t, size_t> > excluded_set_type;

class molecule {

protected : 
    int nsites; //# of sites
    int realsites;
    double* sitecrds; // site coordinates
    double* atmpolar; //atomic polarizabilities
    double* charge;
    double* polfac; //damping factor
    //specify the excluded pairs here
    excluded_set_type excluded12;
    excluded_set_type excluded13;
    excluded_set_type excluded14;

public :

    double* memory;
    int is_w; 

    // sets site positions based on atom coordinates
    virtual double* set_sitecrds(double* xyz) = 0;
    // sets geometry-dependent charges
    virtual double* set_charges(double* xyz) = 0;
    // sets polfac
    virtual double* set_polfacs(double* atmpolar) = 0;
    // sets polarizabilities
    virtual double* set_pol() = 0;
    // allocates memory
    virtual void allocate() = 0;

    virtual double* get_sitecrds() = 0;
    virtual int get_nsites() = 0;
    virtual int get_realsites() = 0;
    virtual double* get_charges() = 0;
    virtual double* get_polfacs() = 0;
    virtual double* get_pol() = 0;

    virtual excluded_set_type::iterator get_begin_12() = 0;
    virtual excluded_set_type::iterator get_begin_13() = 0;
    virtual excluded_set_type::iterator get_begin_14() = 0;
    virtual excluded_set_type::iterator get_end_12() = 0;
    virtual excluded_set_type::iterator get_end_13() = 0;
    virtual excluded_set_type::iterator get_end_14() = 0;


};

#endif // MOLECULE_H
