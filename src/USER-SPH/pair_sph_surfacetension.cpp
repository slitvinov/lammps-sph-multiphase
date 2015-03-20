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
#include "pair_sph_surfacetension.h"
#include "sph_kernel_quintic.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "neigh_list.h"
#include "domain.h"
#include <iostream>

using namespace LAMMPS_NS;

#define EPSILON 1.0e-12

/* ---------------------------------------------------------------------- */

PairSPHSurfaceTension::PairSPHSurfaceTension(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
}

/* ---------------------------------------------------------------------- */

PairSPHSurfaceTension::~PairSPHSurfaceTension() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHSurfaceTension::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih;
  double rsq, wfd;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  double **surface_stress = atom->surface_stress;
  const int ndim = domain->dimension;
  double eij[ndim];
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms and do surface tension
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = rmass[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = rmass[j];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        if (ndim == 3) {
	  wfd = sph_dw_quintic3d(sqrt(rsq)*ih);
          wfd = wfd * ih * ih * ih * ih;
        } else {
	  wfd = sph_dw_quintic2d(sqrt(rsq)*ih);
          wfd = wfd * ih * ih * ih;
        }

	eij[0] = delx/sqrt(rsq); 
	eij[1] = dely/sqrt(rsq);    
	eij[2] = delz/sqrt(rsq);

	const double sigmai = rho[i]/imass;
	const double sigmaj = rho[j]/jmass;
	double fx, fy, fz;
	double* Pj = surface_stress[j];
	double* Pi = surface_stress[i];
	if (ndim == 3) {
	  fx=((eij[0]*Pj[0]+eij[1]*Pj[3]+eij[2]*Pj[4])*pow(sigmai,2)+(eij[0]*Pi[0]+eij[1]*Pi[3]+eij[2]*Pi[4])*pow(sigmaj,2))/pow(sigmai,2)/pow(sigmaj,2);
	  fy=((eij[1]*Pj[1]+eij[0]*Pj[3]+eij[2]*Pj[5])*pow(sigmai,2)+(eij[1]*Pi[1]+eij[0]*Pi[3]+eij[2]*Pi[5])*pow(sigmaj,2))/pow(sigmai,2)/pow(sigmaj,2);
	  fz=((eij[2]*Pj[2]+eij[0]*Pj[4]+eij[1]*Pj[5])*pow(sigmai,2)+(eij[2]*Pi[2]+eij[0]*Pi[4]+eij[1]*Pi[5])*pow(sigmaj,2))/pow(sigmai,2)/pow(sigmaj,2);
	} else {
	  fx=((eij[0]*Pj[0]+eij[1]*Pj[3])*pow(sigmai,2)+(eij[0]*Pi[0]+eij[1]*Pi[3])*pow(sigmaj,2))/pow(sigmai,2)/pow(sigmaj,2);
	  fy=((eij[1]*Pj[1]+eij[0]*Pj[3])*pow(sigmai,2)+(eij[1]*Pi[1]+eij[0]*Pi[3])*pow(sigmaj,2))/pow(sigmai,2)/pow(sigmaj,2);
	  fz = 0.0;
	}

	f[i][0] += fx*wfd;
	f[i][1] += fy*wfd;
	f[i][2] += fz*wfd;	
        if (newton_pair || j < nlocal) {
	  f[j][0] -= fx*wfd;
	  f[j][1] -= fy*wfd;
	  f[j][2] -= fz*wfd;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairSPHSurfaceTension::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairSPHSurfaceTension::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
	       "Illegal number of setting arguments for pair_style sph/surfacetension");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSPHSurfaceTension::coeff(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for pair_style sph/surfacetension coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(arg[0], atom->ntypes, ilo, ihi);
  force->bounds(arg[1], atom->ntypes, jlo, jhi);

  double cut_one   = force->numeric(FLERR, arg[2]);
 
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
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

double PairSPHSurfaceTension::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/surfacetension coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHSurfaceTension::single(int i, int j, int itype, int jtype,
				     double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;
  return 0.0;
}

