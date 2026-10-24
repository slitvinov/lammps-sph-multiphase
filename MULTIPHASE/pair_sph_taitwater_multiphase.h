#ifdef PAIR_CLASS
PairStyle(sph/taitwater/multiphase,PairSPHTaitwaterMultiphase);
#else

#ifndef LMP_PAIR_SPH_TAITWATER_MULTIPHASE_H
#define LMP_PAIR_SPH_TAITWATER_MULTIPHASE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHTaitwaterMultiphase : public Pair {
 public:
  PairSPHTaitwaterMultiphase(class LAMMPS *);
  ~PairSPHTaitwaterMultiphase() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

 protected:
  double *rho0, *soundspeed, *B, *rbackground, *gamma;
  double **cut, **viscosity;
  void allocate();
};

}

#endif
#endif
