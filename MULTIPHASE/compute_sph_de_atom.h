#ifdef COMPUTE_CLASS
ComputeStyle(sph/de/atom,ComputeSPHDeAtom);
#else

#ifndef LMP_COMPUTE_SPH_DE_ATOM_H
#define LMP_COMPUTE_SPH_DE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSPHDeAtom : public Compute {
 public:
  ComputeSPHDeAtom(class LAMMPS *, int, char **);
  ~ComputeSPHDeAtom() override;
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax;
  double *devector;
};

}

#endif
#endif
