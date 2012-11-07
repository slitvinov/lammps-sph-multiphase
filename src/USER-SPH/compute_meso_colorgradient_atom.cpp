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

#include "string.h"
#include "compute_meso_colorgradient_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMesoColorGradientAtom::ComputeMesoColorGradientAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute meso_colorgradient/atom command");
  if (atom->colorgradient_flag != 1) 
    error->all(FLERR,"compute meso_colorgradient/atom command requires atom_style with colorgradient (e.g. meso)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  colorgradientVector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeMesoColorGradientAtom::~ComputeMesoColorGradientAtom()
{
  memory->sfree(colorgradientVector);
}

/* ---------------------------------------------------------------------- */

void ComputeMesoColorGradientAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"colorgradientVector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute colorgradientVector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeMesoColorGradientAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow colorgradientVector array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(colorgradientVector);
    nmax = atom->nmax;
    colorgradientVector = (double *) memory->smalloc(nmax*sizeof(double),"atom:colorgradientVector");
    vector_atom = colorgradientVector;
  }

  // compute kinetic energy for each atom in group

  double **colorgradient = atom->colorgradient;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	if (domain->dimension == 3) {
	  colorgradientVector[i] = sqrt(colorgradient[i][0]*colorgradient[i][0] +
					colorgradient[i][1]*colorgradient[i][1] +
					colorgradient[i][2]*colorgradient[i][2]);
	} else {
	  colorgradientVector[i] = sqrt(colorgradient[i][0]*colorgradient[i][0] +
	  				colorgradient[i][1]*colorgradient[i][1]);
	}
      }
      else {
      	colorgradientVector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMesoColorGradientAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
