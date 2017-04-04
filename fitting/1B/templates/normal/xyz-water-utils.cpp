#include <cassert>
#include <algorithm>

#include "xyz-water-utils.h"

namespace {

void swap_atoms(size_t i, size_t j,
                std::vector<std::string>& elements, std::vector<double>& xyz)
{
    if (i != j) {
        std::swap(elements[i], elements[j]);
        for (size_t k = 0; k < 3; ++k)
            std::swap(xyz[3*i + k], xyz[3*j + k]);
    }
}

bool is_oxygen(const std::string& name)
{
    const size_t a = name.find_first_not_of(' ');
    assert(a != std::string::npos);

    return (name[a] == 'O' || name[a] == 'o' || name[a] == '8');
}

bool is_hydrogen(const std::string& name)
{
    const size_t a = name.find_first_not_of(' ');
    assert(a != std::string::npos);

    return (name[a] == 'H' || name[a] == 'h' || name[a] == '1');
}

} // namespace

namespace kit {

bool is_water(std::vector<std::string>& elements,
              std::vector<double>&      xyz,
              bool reorder, const double& rOH)
{
    const size_t natom = elements.size();
    if (natom%3 != 0 || natom == 0)
        return false;

    if (reorder) {
        // place Oxygens at the top of the xyz and Hydrogens at the bottom
        size_t iO = 0;

        for (size_t n = 0; n < natom; ++n)
            if (is_oxygen(elements[n]))
                swap_atoms(n, iO++, elements, xyz);

        if (iO != natom/3)
            return false;

        for (size_t n = iO; n < natom; ++n)
            if (!is_hydrogen(elements[n]))
                return false;

        // check geometry and rearrange the atoms
        // as O1 .. ON H1a H1b ... HNa HNb

        const double RcSQ = rOH*rOH;
        const size_t nw = natom/3;

        for (size_t i = 0; i < nw; ++i) {
            size_t i_conn = 0;
            for (size_t j = nw; j < natom; ++j) {
                double Rsq(0);
                for (size_t k = 0; k < 3; ++k) {
                    const double delta = xyz[3*i + k] - xyz[3*j + k];
                    Rsq += delta*delta;
                }

                // for the i-th Oxygen, the Hydrogens start at nw + 2*i
                if (Rsq < RcSQ) {
                    const size_t idx = nw + 2*i + i_conn++;

                    if (i_conn <= 2)
                        swap_atoms(j, idx, elements, xyz);
                }
            }

            if (i_conn != 2)
                return false;
        }

        // rearrange the atoms as Oa Ha1 Ha2 Ob Hb1 Hb2 ...

        for (size_t n = nw - 1; n != 0; --n) {
            size_t  first = n;
            size_t middle = n + 1;
            size_t   last = 3*n + 1;

            size_t next = middle;
            while (first != next) {
                swap_atoms(first++, next++, elements, xyz);
                if (next == last)
                    next = middle;
                else if (first == middle)
                    middle = next;
            }
        }

        // swap some hydrogens so that the output is the same
        // as input for the case when the latter is O H H O H H ...

        for (size_t n = 1 + 3*(nw%2); n < 3*nw; n += 6)
            swap_atoms(n, n + 1, elements, xyz);

        return true;

    } else {
        for (size_t n = 0; n < natom; n += 3)
            if (!is_oxygen(elements[n])
             || !is_hydrogen(elements[n + 1])
             || !is_hydrogen(elements[n + 2]))
                return false;
        return true;
    }
}

void make_water_elements(size_t nw, std::vector<std::string>& elements)
{
    elements.resize(0);

    for (size_t n = 0; n < nw; ++n) {
        elements.push_back("O");
        elements.push_back("H");
        elements.push_back("H");
    }
}

} // namespace kit
