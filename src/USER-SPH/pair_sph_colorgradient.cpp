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
#include "stdlib.h"
#include "pair_sph_colorgradient.h"
#include "sph_kernel_quintic.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"
#include <iostream>
#include <assert.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHColorGradient::PairSPHColorGradient(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair

  comm_forward = 3;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHColorGradient::~PairSPHColorGradient() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(alpha);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairSPHColorGradient::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHColorGradient::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, h, ih, ihsq;
  int *jlist;
  
  const int ndim = domain->dimension;
  double eij[ndim];
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *rho = atom->rho;
  double **colorgradient = atom->colorgradient;
  int *type = atom->type;
  double *rmass = atom->rmass;

  // check consistency of pair coefficients

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 0.0) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH particle types %d and %d interact, but not all of their single particle properties are set.\n",
                  i, j);
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density
  // we use a full neighborlist here

  if (nstep != 0) {
    if ((update->ntimestep % nstep) == 0) {

      // initialize color gradient with zeros
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        colorgradient[i][0] = 0.0;
        colorgradient[i][1] = 0.0;
	colorgradient[i][2] = 0.0;
      } // ii loop

      // add density at each atom via kernel function overlap
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];
	double Vi = rmass[i]/rho[i]; 
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;

          jtype = type[j];
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

	  if (rsq < cutsq[itype][jtype]) {
	    double r = sqrt(rsq);
	    if (ndim==2) {
	      eij[0]= delx/r; 
	      eij[1]= dely/r;
	    } else {
	      eij[0]= delx/r;
	      eij[1]= dely/r;
	      eij[2]= delz/r;
	    }

	    h = cut[itype][jtype];
	    ih = 1.0 / h;

	    double wfd;
	    if (domain->dimension == 3) {
	      // Quintic spline
	      wfd = sph_dw_quintic3d(r*ih);
	      wfd = wfd * ih * ih * ih * ih;
	    } else {
	      wfd = sph_dw_quintic2d(r*ih);
	      wfd = wfd * ih * ih * ih;
	    }
	    double Vj = rmass[j]/rho[j];
	    double Vj2 = Vj*Vj;
	    double dphi = -wfd*alpha[itype][jtype]*dphi*Vj2/Vi;
	    
	    colorgradient[i][0] += dphi*eij[0];
	    colorgradient[i][1] += dphi*eij[1];
	    if (ndim==3) {
	      colorgradient[i][2] += dphi*eij[2];
	    }
          } // rsq < cutsq
	  
        } // jj loop
      } // ii loop
    }
  }

  // communicate densities
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHColorGradient::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHColorGradient::settings(int narg, char **arg) {
  if (narg != 1)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/colorgradient");
  nstep = force->inumeric(FLERR, arg[0]);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHColorGradient::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,"Incorrect number of args for sph/colorgradient coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(arg[0], atom->ntypes, ilo, ihi);
  force->bounds(arg[1], atom->ntypes, jlo, jhi);

  double cut_one = force->numeric(FLERR, arg[2]);
  double alpha_one = force->numeric(FLERR, arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      alpha[i][j] = alpha_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHColorGradient::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/colorgradient coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  alpha[j][i] = alpha[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHColorGradient::single(int i, int j, int itype, int jtype, double rsq,
    double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHColorGradient::pack_comm(int n, int *list, double *buf, int pbc_flag,
    int *pbc) {
  int i, j, m;
  double **colorgradient = atom->colorgradient;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = colorgradient[j][0];
    buf[m++] = colorgradient[j][1];
    buf[m++] = colorgradient[j][2];
  }
  return 3;
}

/* ---------------------------------------------------------------------- */

void PairSPHColorGradient::unpack_comm(int n, int first, double *buf) {
  int i, m, last;
  double **colorgradient = atom->colorgradient;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    colorgradient[i][0] = buf[m++];
    colorgradient[i][1] = buf[m++];
    colorgradient[i][2] = buf[m++];
}

