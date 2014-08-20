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
  double imass, jmass, h, ih, ihsq;
  double rsq, wfd;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  double **cg = atom->colorgradient;
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

    double abscgi;
    if (ndim == 3) {
      abscgi = sqrt(cg[i][0]*cg[i][0] +
		    cg[i][1]*cg[i][1] +
		    cg[i][2]*cg[i][2]);
    } else {
      abscgi = sqrt(cg[i][0]*cg[i][0] +
		    cg[i][1]*cg[i][1]);
    }
    // TODO: FixMe
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
	if (ndim==3) {
	  eij[2] = delz/sqrt(rsq);
	}

	double SurfaceForcei[ndim];
	double SurfaceForcej[ndim];
	SurfaceForcei[0]=0; SurfaceForcej[0]=0;
	SurfaceForcei[1]=0; SurfaceForcej[1]=0;
	if (ndim==3) {
	  SurfaceForcei[2]=0; SurfaceForcej[2]=0;
	}

	if (ndim==2) {
	  /// TODO: can be moved outside of the jj loop
	  double abscgj = sqrt(cg[j][0]*cg[j][0] + cg[j][1]*cg[j][1]);
	  if (abscgi > EPSILON) {
	    SurfaceForcei[0] = (eij[0]*((cg[i][1]*cg[i][1]+cg[i][0]*cg[i][0])/2-cg[i][0]*cg[i][0])-cg[i][0]*eij[1]*cg[i][1])/abscgi;
	    SurfaceForcei[1] = (eij[1]*((cg[i][1]*cg[i][1]+cg[i][0]*cg[i][0])/2-cg[i][1]*cg[i][1])-eij[0]*cg[i][0]*cg[i][1])/abscgi;
	  }

	  if (abscgj > EPSILON) {
	    SurfaceForcej[0] = (eij[0]*((cg[j][1]*cg[j][1]+cg[j][0]*cg[j][0])/2-cg[j][0]*cg[j][0])-cg[j][0]*eij[1]*cg[j][1])/abscgj;
	    SurfaceForcej[1] = (eij[1]*((cg[j][1]*cg[j][1]+cg[j][0]*cg[j][0])/2-cg[j][1]*cg[j][1])-eij[0]*cg[j][0]*cg[j][1])/abscgj;
	  }
	} else {
	  if (abscgi > EPSILON) {
	    SurfaceForcei[0] = (eij[0]*((cg[i][2]*cg[i][2]+cg[i][1]*cg[i][1]+cg[i][0]*cg[i][0])/3-cg[i][0]*cg[i][0])
				-cg[i][0]*eij[2]*cg[i][2]-cg[i][0]*eij[1]*cg[i][1])/ abscgi;
	    SurfaceForcei[1] = (eij[1]*((cg[i][2]*cg[i][2]+cg[i][1]*cg[i][1]+cg[i][0]*cg[i][0])/3-cg[i][1]*cg[i][1])
				-cg[i][1]*eij[2]*cg[i][2]-eij[0]*cg[i][0]*cg[i][1] ) / abscgi;
	    SurfaceForcei[2] = ( eij[2]*((cg[i][2]*cg[i][2]+cg[i][1]*cg[i][1]+cg[i][0]*cg[i][0])/3-cg[i][2]*cg[i][2])
				 -eij[1]*cg[i][1]*cg[i][2]-eij[0]*cg[i][0]*cg[i][2] ) /abscgi;
	  }
	  double abscgj = sqrt(cg[j][0]*cg[j][0] + cg[j][1]*cg[j][1] + cg[j][2]*cg[j][2]);
	  if (abscgj > EPSILON) {
	    SurfaceForcej[0] = (eij[0]*((cg[j][2]*cg[j][2]+cg[j][1]*cg[j][1]+cg[j][0]*cg[j][0])/3-cg[j][0]*cg[j][0])
				-cg[j][0]*eij[2]*cg[j][2]-cg[j][0]*eij[1]*cg[j][1])/ abscgj;
	    SurfaceForcej[1] = (eij[1]*((cg[j][2]*cg[j][2]+cg[j][1]*cg[j][1]+cg[j][0]*cg[j][0])/3-cg[j][1]*cg[j][1])
				-cg[j][1]*eij[2]*cg[j][2]-eij[0]*cg[j][0]*cg[j][1] ) / abscgj;
	    SurfaceForcej[2] = ( eij[2]*((cg[j][2]*cg[j][2]+cg[j][1]*cg[j][1]+cg[j][0]*cg[j][0])/3-cg[j][2]*cg[j][2])
				 -eij[1]*cg[j][1]*cg[j][2]-eij[0]*cg[j][0]*cg[j][2] ) /abscgj;
	  }
	}

	const double Vi = imass / rho[i];
	const double Vj = jmass / rho[j];
	//	const double rij = sqrt(rsq);	    

	f[i][0] += (SurfaceForcei[0]*Vi*Vi + SurfaceForcej[0]*Vj*Vj)*wfd;
	f[i][1] += (SurfaceForcei[1]*Vi*Vi + SurfaceForcej[1]*Vj*Vj)*wfd;
	if (ndim==3) {
	  f[i][2] += (SurfaceForcei[2]*Vi*Vi + SurfaceForcej[2]*Vj*Vj)*wfd;
	}

        if (newton_pair || j < nlocal) {
	  f[j][0] -= (SurfaceForcei[0]*Vi*Vi + SurfaceForcej[0]*Vj*Vj)*wfd;
	  f[j][1] -= (SurfaceForcei[1]*Vi*Vi + SurfaceForcej[1]*Vj*Vj)*wfd;
	  if (ndim==3) {
	    f[j][2] -= (SurfaceForcei[2]*Vi*Vi + SurfaceForcej[2]*Vj*Vj)*wfd;
	  }
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

/* calculate phase stress base on phase gradient */
void get_phase_stress(double* v, double* del_phi) {
  double interm0 = 1.0/ ( sqrt(v[0]*v[0] + v[1]*v[1]) + EPSILON );
  double interm1 = 0.5 * (v[0]*v[0] - v[1]*v[1]);
  double interm2 = v[0]*v[1];
  del_phi[0] = interm1*interm0;
  del_phi[1] = interm2*interm0;
}
