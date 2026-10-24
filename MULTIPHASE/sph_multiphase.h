#ifndef LMP_SPH_MULTIPHASE_H
#define LMP_SPH_MULTIPHASE_H

#include "atom.h"
#include "error.h"

#include <cmath>

namespace LAMMPS_NS {
namespace SPHMultiphase {

  static constexpr double norm3d = 0.0716197243913529;
  static constexpr double norm2d = 0.04195297663091802;

  inline double kernel_poly(double s)
  {
    if (s < 1.0) return pow(3 - s, 5) - 6 * pow(2 - s, 5) + 15 * pow(1 - s, 5);
    if (s < 2.0) return pow(3 - s, 5) - 6 * pow(2 - s, 5);
    if (s < 3.0) return pow(3 - s, 5);
    return 0.0;
  }

  inline double dkernel_poly(double s)
  {
    if (s < 1.0) return -50 * pow(s, 4) + 120 * pow(s, 3) - 120 * s;
    if (s < 2.0) return 25 * pow(s, 4) - 180 * pow(s, 3) + 450 * pow(s, 2) - 420 * s + 75;
    if (s < 3.0) return -5 * pow(s, 4) + 60 * pow(s, 3) - 270 * pow(s, 2) + 540 * s - 405;
    return 0.0;
  }

  inline double w(int dim, double r)
  {
    return (dim == 3 ? norm3d : norm2d) * kernel_poly(3.0 * r);
  }

  inline double dw(int dim, double r)
  {
    return 3.0 * (dim == 3 ? norm3d : norm2d) * dkernel_poly(3.0 * r);
  }

  inline double energy2t(double e, double cv) { return e / cv; }

  inline double t2energy(double t, double cv) { return cv * t; }

  inline double pressure(double B, double rho0, double gamma, double rbackground, double rho)
  {
    return B * (pow(rho / rho0, gamma) - rbackground);
  }

  inline double **colorgradient(Atom *atom, Error *error)
  {
    int flag, cols;
    int index = atom->find_custom("cg", flag, cols);
    if (index < 0 || flag != 1 || cols != 3)
      error->all(FLERR, "SPH multiphase requires fix property/atom d2_cg 3 ghost yes");
    return atom->darray[index];
  }

  inline void require_rmass(Atom *atom, Error *error)
  {
    if (!atom->rmass_flag || !atom->rho_flag || !atom->esph_flag)
      error->all(FLERR,
                 "SPH multiphase requires atom_style sph and fix property/atom rmass ghost yes");
  }

}
}

#endif
