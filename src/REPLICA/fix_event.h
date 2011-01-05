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

#ifndef LMP_FIX_EVENT_H
#define LMP_FIX_EVENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEvent : public Fix {
 public:

  FixEvent(class LAMMPS *, int, char **);
  virtual ~FixEvent()=0;    // Use destructor to make base class virtual
  int setmask();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  virtual void write_restart(FILE *);
  virtual void restart(char *);

  // methods specific to FixEvent

  void store_event();         // store quenched atoms
  void restore_event();       // restore quenched atoms
  void store_state();         // store hot atoms
  void restore_state();       // restore hot atoms

 private:
  double **xevent;       // atom coords at last event
  double **xold;         // atom coords for reset/restore
  double **vold;         // atom vels for reset/restore
  int *imageold;         // image flags for reset/restore
};

}

#endif
