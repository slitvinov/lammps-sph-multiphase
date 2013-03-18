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

#ifdef PAIR_CLASS

PairStyle(sph/taitwater/multiphase,PairSPHTaitwaterMultiphase)

#else

#ifndef LMP_PAIR_TAITWATER_MULTIPHASE_H
#define LMP_PAIR_TAITWATER_MULTIPHASE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHTaitwaterMultiphase : public Pair {
 public:
  PairSPHTaitwaterMultiphase(class LAMMPS *);
  virtual ~PairSPHTaitwaterMultiphase();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double *rho0, *soundspeed, *B, *rbackground, *gamma;
  double **cut,**viscosity;
  int first;

  void allocate();
};

  // equation of state (helper function)
  double sph_pressure(double B, double rho0, double gamma, double rbackground, double rho);
}


#endif
#endif
