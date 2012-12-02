/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include "math.h"
#include "sph_kernel_quintic.h"

using namespace LAMMPS_NS;

double sph_kernel_quintic_2d(double r) {
  double norm2d = 0.041952976630918;
  double s = 3.0*r;
  if (s<1.0) {
    return norm2d*(pow(3 - s, 5) - 6*pow(2 - s, 5) + 15*pow(1 - s, 5));
  } else if (s<2.0) {
    return norm2d*(pow(3 - s, 5) - 6*pow(2 - s, 5));
  } else {
    return norm2d*pow(3 - s, 5);
  }
  return 0.0;
};

double sph_kernel_quintic_3d(double r) {
  double norm3d = 0.0026636810559313;
  double s = 3.0*r;
  if (s<1.0) {
    return norm3d*(pow(3 - s, 5) - 6*pow(2 - s, 5) + 15*pow(1 - s, 5));
  } else if (s<2.0) {
    return norm3d*(pow(3 - s, 5) - 6*pow(2 - s, 5));
  } else {
    return norm3d*pow(3 - s, 5);
  }
  return 0.0;
}
