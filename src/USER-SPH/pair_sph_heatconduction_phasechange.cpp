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
#include "assert.h"
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "pair_sph_heatconduction_phasechange.h"
#include "sph_kernel_quintic.h"
#include "sph_energy_equation.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "neigh_list.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHHeatConductionPhaseChange::PairSPHHeatConductionPhaseChange(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
}

/* ---------------------------------------------------------------------- */

PairSPHHeatConductionPhaseChange::~PairSPHHeatConductionPhaseChange() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(alpha);
    memory->destroy(tc);
    memory->destroy(fixflag);
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHHeatConductionPhaseChange::compute(int eflag, int vflag) {
  int ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih;
  double rsq, wfd, D;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *e = atom->e;
  double *de = atom->de;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  double *cv = atom->cv;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms and do heat diffusion

  for (ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    itype = type[i];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = rmass[i];

    for (jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        // kernel function
        if (domain->dimension == 3) {
	  wfd = sph_dw_quintic3d(sqrt(rsq)*ih);
          wfd = wfd * ih * ih * ih * ih / sqrt(rsq);
        } else {
	  wfd = sph_dw_quintic2d(sqrt(rsq)*ih);
          wfd = wfd * ih * ih * ih / sqrt(rsq);
        }

        jmass = rmass[j];
	assert(jmass>0);
        D = alpha[itype][jtype]; // diffusion coefficient

	double Ti = sph_energy2t(e[i], cv[i]);
	double Tj = sph_energy2t(e[j], cv[j]);
	
	if ( (fixflag[itype][jtype]==itype) && (Ti<Tj))  {
	  Ti = tc[itype][jtype];
	}
	if ( (fixflag[itype][jtype]==jtype ) && (Tj<Ti)){
	  Tj = tc[itype][jtype];
	}
	assert(rho[i]>0);
	assert(rho[j]>0);
        double deltaE = 2.0*D*(Ti - Tj)*wfd/(rho[i]*rho[j]);
        de[i] += deltaE*jmass;
        if (newton_pair || j < nlocal) {
          de[j] -= deltaE*imass;
        }

      }
    }
  }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHHeatConductionPhaseChange::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");
  memory->create(tc, n + 1, n + 1, "pair:tc");
  memory->create(fixflag, n + 1, n + 1, "pair:fixflag");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHHeatConductionPhaseChange::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/heatconduction/phasechange");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHHeatConductionPhaseChange::coeff(int narg, char **arg) {
  if ( (narg != 4) && (narg != 6) )
    error->all(FLERR,"Incorrect number of args for pair_style sph/heatconduction/phasechange coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(arg[0], atom->ntypes, ilo, ihi);
  force->bounds(arg[1], atom->ntypes, jlo, jhi);

  double alpha_one = force->numeric(FLERR, arg[2]);
  double cut_one   = force->numeric(FLERR, arg[3]);

  // default value
  int fixflag_one = 0;
  double tc_one = 0;
  if (narg==6) {
    // process "fix temperature" options
    if (strcmp(arg[4],"NULL") == 0) {
      fixflag_one = 2; // second type has fixed temperature
      tc_one = force->numeric(FLERR, arg[5]);
    } else if (strcmp(arg[5],"NULL") == 0) {
      fixflag_one = 1; // first type has fixed temperature
      tc_one = force->numeric(FLERR, arg[4]);
    } else {
      error->all(FLERR,"Incorrect args for pair coefficients");
    }
  }
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      alpha[i][j] = alpha_one;
      setflag[i][j] = 1;
      if (fixflag_one == 1) {
	fixflag[i][j] = i;
	tc[i][j] = tc_one;
      } else if (fixflag_one==2) {
	fixflag[i][j] = j;
	tc[i][j] = tc_one;
      }
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHHeatConductionPhaseChange::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/heatconduction/phasechange coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  alpha[j][i] = alpha[i][j];
  tc[j][i] = tc[i][j];
  fixflag[j][i] = fixflag[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHHeatConductionPhaseChange::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}
