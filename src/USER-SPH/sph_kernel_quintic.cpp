/* ----------------------------------------------------------------------
 LAMMPS-Large-scale Atomic/Molecular Massively Parallel Simulator
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

double LAMMPS_NS::sph_kernel_quintic3d(double r) {
  const double norm3d = 0.0716197243913529;
  const double s = 3.0*r;
  if (s<1.0) {
    return norm3d*(pow(3-s,5)-6*pow(2-s,5) + 15*pow(1-s,5));
  } else if (s<2.0) {
    return norm3d*(pow(3-s,5)-6*pow(2-s,5));
  } else if (s<3.0) {
    return norm3d*pow(3-s,5);
  }
  return 0.0;
}

double LAMMPS_NS::sph_kernel_quintic2d(double r) {
  const double norm2d = 0.04195297663091802;
  const double s = 3.0*r;
  if (s<1.0) {
    return norm2d*(pow(3-s,5)-6*pow(2-s,5) + 15*pow(1-s,5));
  } else if (s<2.0) {
    return norm2d*(pow(3-s,5)-6*pow(2-s,5));
  } else if (s<3.0) {
    return norm2d*pow(3-s,5);
  }
  return 0.0;
}

double LAMMPS_NS::sph_dw_quintic3d(double r) {
  const double norm3d = 3.0*0.0716197243913529;
  const double s = 3.0*r;
  double wfd;
  if (s<1) {
    wfd = -50*pow(s,4)+120*pow(s,3)-120*s;
  } else if (s<2) {
    wfd = 25*pow(s,4)-180*pow(s,3)+450*pow(s,2)-420*s+75;
  } else if (s<3.0) {
    wfd = -5*pow(s,4)+60*pow(s,3)-270*pow(s,2)+540*s-405;
  } else {
    wfd = 0.0;
  }
  return norm3d*wfd;
}

double LAMMPS_NS::sph_dw_quintic2d(double r) {
  const double norm2d = 3.0*0.04195297663091802;
  const double s = 3.0*r;
  double wfd;
  if (s<1) {
    wfd = -50*pow(s,4)+120*pow(s,3)-120*s;
  } else if (s<2) {
    wfd = 25*pow(s,4)-180*pow(s,3)+450*pow(s,2)-420*s+75;
  } else if (s<3.0) {
    wfd = -5*pow(s,4)+60*pow(s,3)-270*pow(s,2)+540*s-405;
  } else {
    wfd = 0.0;
  }
  return norm2d*wfd;
}

