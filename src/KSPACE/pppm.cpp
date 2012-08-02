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

/* ----------------------------------------------------------------------
   Contributing authors: Roy Pollock (LLNL), Paul Crozier (SNL)
     per-atom energy/virial & group/group energy/force added by Stan Moore (BYU)
     analytic diff (2 FFT) option added by Rolf Isele-Holder (Aachen University)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pppm.h"
#include "math_const.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXORDER 7
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPM::PPPM(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg < 1) error->all(FLERR,"Illegal kspace_style pppm command");

  group_group_enable = 1;

  accuracy_relative = atof(arg[0]);

  nfactors = 3;
  factors = new int[nfactors];
  factors[0] = 2;
  factors[1] = 3;
  factors[2] = 5;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  density_brick = vdx_brick = vdy_brick = vdz_brick = NULL;
  density_fft = NULL;
  u_brick = NULL;
  v0_brick = v1_brick = v2_brick = v3_brick = v4_brick = v5_brick = NULL;
  greensfn = NULL;
  work1 = work2 = NULL;
  vg = NULL;
  fkx = fky = fkz = NULL;
  buf1 = buf2 = buf3 = buf4 = NULL;

  density_A_brick = density_B_brick = NULL;
  density_A_fft = density_B_fft = NULL;

  gf_b = NULL;
  rho1d = rho_coeff = drho1d = drho_coeff = NULL;

  fft1 = fft2 = NULL;
  remap = NULL;

  nmax = 0;
  part2grid = NULL;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPM::~PPPM()
{
  delete [] factors;
  deallocate();
  deallocate_peratom();
  deallocate_groups();
  memory->destroy(part2grid);
  memory->destroy(acons);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPM::init()
{
  if (me == 0) {
    if (screen) fprintf(screen,"PPPM initialization ...\n");
    if (logfile) fprintf(logfile,"PPPM initialization ...\n");
  }

  // error check

  if (domain->triclinic)
    error->all(FLERR,"Cannot (yet) use PPPM with triclinic box");
  if (domain->dimension == 2) error->all(FLERR,
                                         "Cannot use PPPM with 2d simulation");

  if (!atom->q_flag) error->all(FLERR,"Kspace style requires atom attribute q");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use nonperiodic boundaries with PPPM");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPM");
  }

  if (order < 2 || order > MAXORDER) {
    char str[128];
    sprintf(str,"PPPM order cannot be < 2 or > than %d",MAXORDER);
    error->all(FLERR,str);
  }

  // free all arrays previously allocated

  deallocate();
  deallocate_peratom();
  peratom_allocate_flag = 0;
  deallocate_groups();
  group_allocate_flag = 0;

  // extract short-range Coulombic cutoff from pair style

  scale = 1.0;

  if (force->pair == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  int itmp=0;
  double *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  if (p_cutoff == NULL)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // if kspace is TIP4P, extract TIP4P params from pair style
  // bond/angle are not yet init(), so insure equilibrium request is valid

  qdist = 0.0;

  if ( (strcmp(force->kspace_style,"pppm/tip4p") == 0) ||
       (strcmp(force->kspace_style,"pppm/tip4p/proxy") == 0) ) {
    if (force->pair == NULL)
      error->all(FLERR,"KSpace style is incompatible with Pair style");
    double *p_qdist = (double *) force->pair->extract("qdist",itmp);
    int *p_typeO = (int *) force->pair->extract("typeO",itmp);
    int *p_typeH = (int *) force->pair->extract("typeH",itmp);
    int *p_typeA = (int *) force->pair->extract("typeA",itmp);
    int *p_typeB = (int *) force->pair->extract("typeB",itmp);
    if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
      error->all(FLERR,"KSpace style is incompatible with Pair style");
    qdist = *p_qdist;
    typeO = *p_typeO;
    typeH = *p_typeH;
    int typeA = *p_typeA;
    int typeB = *p_typeB;

    if (force->angle == NULL || force->bond == NULL)
      error->all(FLERR,"Bond and angle potentials must be defined for TIP4P");
    if (typeA < 1 || typeA > atom->nangletypes ||
        force->angle->setflag[typeA] == 0)
      error->all(FLERR,"Bad TIP4P angle type for PPPM/TIP4P");
    if (typeB < 1 || typeB > atom->nbondtypes ||
        force->bond->setflag[typeB] == 0)
      error->all(FLERR,"Bad TIP4P bond type for PPPM/TIP4P");
    double theta = force->angle->equilibrium_angle(typeA);
    double blen = force->bond->equilibrium_distance(typeB);
    alpha = qdist / (cos(0.5*theta) * blen);
  }

  // if we have a /proxy pppm version check if the pair style is compatible

  if ((strcmp(force->kspace_style,"pppm/proxy") == 0) ||
      (strcmp(force->kspace_style,"pppm/tip4p/proxy") == 0) ) {
    if (force->pair == NULL)
      error->all(FLERR,"KSpace style is incompatible with Pair style");
    if (strstr(force->pair_style,"pppm/") == NULL )
      error->all(FLERR,"KSpace style is incompatible with Pair style");
  }

  // compute qsum & qsqsum and warn if not charge-neutral

  qsum = qsqsum = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    qsum += atom->q[i];
    qsqsum += atom->q[i]*atom->q[i];
  }

  double tmp;
  MPI_Allreduce(&qsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum = tmp;
  MPI_Allreduce(&qsqsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsqsum = tmp;
  q2 = qsqsum * force->qqrd2e / force->dielectric;

  if (qsqsum == 0.0)
    error->all(FLERR,"Cannot use kspace solver on system with no charge");
  if (fabs(qsum) > SMALL && me == 0) {
    char str[128];
    sprintf(str,"System is not charge neutral, net charge = %g",qsum);
    error->warning(FLERR,str);
  }

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // define acons coefficients for estimation of kspace errors
  // see JCP 109, pg 7698 for derivation of coefficients
  // higher order coefficients may be computed if needed

  memory->destroy(acons);
  memory->create(acons,8,7,"pppm:acons");
  acons[1][0] = 2.0 / 3.0;
  acons[2][0] = 1.0 / 50.0;
  acons[2][1] = 5.0 / 294.0;
  acons[3][0] = 1.0 / 588.0;
  acons[3][1] = 7.0 / 1440.0;
  acons[3][2] = 21.0 / 3872.0;
  acons[4][0] = 1.0 / 4320.0;
  acons[4][1] = 3.0 / 1936.0;
  acons[4][2] = 7601.0 / 2271360.0;
  acons[4][3] = 143.0 / 28800.0;
  acons[5][0] = 1.0 / 23232.0;
  acons[5][1] = 7601.0 / 13628160.0;
  acons[5][2] = 143.0 / 69120.0;
  acons[5][3] = 517231.0 / 106536960.0;
  acons[5][4] = 106640677.0 / 11737571328.0;
  acons[6][0] = 691.0 / 68140800.0;
  acons[6][1] = 13.0 / 57600.0;
  acons[6][2] = 47021.0 / 35512320.0;
  acons[6][3] = 9694607.0 / 2095994880.0;
  acons[6][4] = 733191589.0 / 59609088000.0;
  acons[6][5] = 326190917.0 / 11700633600.0;
  acons[7][0] = 1.0 / 345600.0;
  acons[7][1] = 3617.0 / 35512320.0;
  acons[7][2] = 745739.0 / 838397952.0;
  acons[7][3] = 56399353.0 / 12773376000.0;
  acons[7][4] = 25091609.0 / 1560084480.0;
  acons[7][5] = 1755948832039.0 / 36229939200000.0;
  acons[7][6] = 4887769399.0 / 37838389248.0;
  
  // setup FFT grid resolution and g_ewald
  // normally one iteration thru while loop is all that is required
  // if grid stencil extends beyond neighbor proc, reduce order and try again

  int iteration = 0;
  triclinic = domain->triclinic;

  while (order > 1) {
    if (iteration && me == 0)
      error->warning(FLERR,"Reducing PPPM order b/c stencil extends "
                     "beyond neighbor processor");
    iteration++;

    set_grid();

    if (nx_pppm >= OFFSET || ny_pppm >= OFFSET || nz_pppm >= OFFSET)
      error->all(FLERR,"PPPM grid is too large");

    set_fft_parameters();

    // test that ghost overlap is not bigger than my sub-domain

    int flag = 0;
    if (nxlo_ghost > nxhi_in-nxlo_in+1) flag = 1;
    if (nxhi_ghost > nxhi_in-nxlo_in+1) flag = 1;
    if (nylo_ghost > nyhi_in-nylo_in+1) flag = 1;
    if (nyhi_ghost > nyhi_in-nylo_in+1) flag = 1;
    if (nzlo_ghost > nzhi_in-nzlo_in+1) flag = 1;
    if (nzhi_ghost > nzhi_in-nzlo_in+1) flag = 1;

    int flag_all;
    MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);

    if (flag_all == 0) break;
    order--;
  }

  if (order == 0) error->all(FLERR,"PPPM order has been reduced to 0");

  // adjust g_ewald

  if (!gewaldflag) adjust_gewald();

  // calculate the final accuracy

  double estimated_accuracy = final_accuracy();
  
  // print stats

  int ngrid_max,nfft_both_max,nbuf_max;
  MPI_Allreduce(&ngrid,&ngrid_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nfft_both,&nfft_both_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nbuf,&nbuf_max,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {

  #ifdef FFT_SINGLE
    const char fft_prec[] = "single";
  #else
    const char fft_prec[] = "double";
  #endif

    if (screen) {
      fprintf(screen,"  G vector (1/distance)= %g\n",g_ewald);
      fprintf(screen,"  grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
      fprintf(screen,"  stencil order = %d\n",order);
      fprintf(screen,"  estimated absolute RMS force accuracy = %g\n",
              estimated_accuracy);
      fprintf(screen,"  estimated relative force accuracy = %g\n",
              estimated_accuracy/two_charge_force);
      fprintf(screen,"  using %s precision FFTs\n",fft_prec);
      fprintf(screen,"  brick FFT buffer size/proc = %d %d %d\n",
                        ngrid_max,nfft_both_max,nbuf_max);
    }
    if (logfile) {
      fprintf(logfile,"  G vector (1/distance) = %g\n",g_ewald);
      fprintf(logfile,"  grid = %d %d %d\n",nx_pppm,ny_pppm,nz_pppm);
      fprintf(logfile,"  stencil order = %d\n",order);
      fprintf(logfile,"  estimated absolute RMS force accuracy = %g\n",
              estimated_accuracy);
      fprintf(logfile,"  estimated relative force accuracy = %g\n",
              estimated_accuracy/two_charge_force);
      fprintf(logfile,"  using %s precision FFTs\n",fft_prec);
      fprintf(logfile,"  brick FFT buffer size/proc = %d %d %d\n",
                         ngrid_max,nfft_both_max,nbuf_max);
    }
  }

  // allocate K-space dependent memory
  // don't invoke allocate_peratom() here, wait to see if needed

  allocate();

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  compute_rho_coeff();
}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void PPPM::setup()
{
  int i,j,k,l,m,n;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  delxinv = nx_pppm/xprd;
  delyinv = ny_pppm/yprd;
  delzinv = nz_pppm/zprd_slab;

  delvolinv = delxinv*delyinv*delzinv;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  // fkx,fky,fkz for my FFT grid pts

  double per;

  for (i = nxlo_fft; i <= nxhi_fft; i++) {
    per = i - nx_pppm*(2*i/nx_pppm);
    fkx[i] = unitkx*per;
  }

  for (i = nylo_fft; i <= nyhi_fft; i++) {
    per = i - ny_pppm*(2*i/ny_pppm);
    fky[i] = unitky*per;
  }

  for (i = nzlo_fft; i <= nzhi_fft; i++) {
    per = i - nz_pppm*(2*i/nz_pppm);
    fkz[i] = unitkz*per;
  }

  // virial coefficients

  double sqk,vterm;

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++) {
    for (j = nylo_fft; j <= nyhi_fft; j++) {
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        sqk = fkx[i]*fkx[i] + fky[j]*fky[j] + fkz[k]*fkz[k];
        if (sqk == 0.0) {
          vg[n][0] = 0.0;
          vg[n][1] = 0.0;
          vg[n][2] = 0.0;
          vg[n][3] = 0.0;
          vg[n][4] = 0.0;
          vg[n][5] = 0.0;
        } else {
          vterm = -2.0 * (1.0/sqk + 0.25/(g_ewald*g_ewald));
          vg[n][0] = 1.0 + vterm*fkx[i]*fkx[i];
          vg[n][1] = 1.0 + vterm*fky[j]*fky[j];
          vg[n][2] = 1.0 + vterm*fkz[k]*fkz[k];
          vg[n][3] = vterm*fkx[i]*fky[j];
          vg[n][4] = vterm*fkx[i]*fkz[k];
          vg[n][5] = vterm*fky[j]*fkz[k];
        }
        n++;
      }
    }
  }
  
  if (differentiation_flag == 1) {
    compute_gf_ad();
    compute_sf_coeff();
  } else {
    compute_gf_ik();
  }
}

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPM::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    peratom_allocate_flag = 1;
  }

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  poisson();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  fillbrick();

  // extra per-atom energy/virial communication

  if (evflag_atom) fillbrick_peratom();

  // calculate the force on my particles

  fieldforce();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in volume-dependent term

  const double qscale = force->qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qscale;
  }

  // sum global virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    double *q = atom->q;
    int nlocal = atom->nlocal;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum /
          (g_ewald*g_ewald*volume);
        eatom[i] *= qscale;
      }
    }

    if (vflag_atom) {
      for (i = 0; i < nlocal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPM::allocate()
{
  memory->create3d_offset(density_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_brick");

  memory->create(density_fft,nfft_both,"pppm:density_fft");
  memory->create(greensfn,nfft_both,"pppm:greensfn");
  memory->create(work1,2*nfft_both,"pppm:work1");
  memory->create(work2,2*nfft_both,"pppm:work2");
  memory->create(vg,nfft_both,6,"pppm:vg");

  memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"pppm:fkx");
  memory->create1d_offset(fky,nylo_fft,nyhi_fft,"pppm:fky");
  memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"pppm:fkz");

  memory->create(buf1,nbuf,"pppm:buf1");
  memory->create(buf2,nbuf,"pppm:buf2");

  if (differentiation_flag == 1) {
    memory->create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:u_brick");
  } else {
    memory->create3d_offset(vdx_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdx_brick");
    memory->create3d_offset(vdy_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdy_brick");
    memory->create3d_offset(vdz_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdz_brick");
  }
  
  // summation coeffs

  memory->create(gf_b,order,"pppm:gf_b");
  memory->create2d_offset(rho1d,3,-order/2,order/2,"pppm:rho1d");
  memory->create2d_offset(drho1d,3,-order/2,order/2,"pppm:drho1d");
  memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"pppm:rho_coeff");
  memory->create2d_offset(drho_coeff,order,(1-order)/2,order/2,"pppm:drho_coeff");

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decompostion
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;

  fft1 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   0,0,&tmp);

  fft2 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   0,0,&tmp);

  remap = new Remap(lmp,world,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                    1,0,0,FFT_PRECISION);
}

/* ----------------------------------------------------------------------
   allocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPM::allocate_peratom()
{
  if (differentiation_flag != 1)
    memory->create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:u_brick");

  memory->create3d_offset(v0_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v0_brick");
  memory->create3d_offset(v1_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v1_brick");
  memory->create3d_offset(v2_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v2_brick");
  memory->create3d_offset(v3_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v3_brick");
  memory->create3d_offset(v4_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v4_brick");
  memory->create3d_offset(v5_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v5_brick");

  memory->create(buf3,nbuf_peratom,"pppm:buf3");
  memory->create(buf4,nbuf_peratom,"pppm:buf4");
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPM::deallocate()
{
  memory->destroy3d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);
  if (differentiation_flag == 1) {
    memory->destroy3d_offset(u_brick,nzlo_out,nylo_out,nxlo_out);
  } else {
    memory->destroy3d_offset(vdx_brick,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdy_brick,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdz_brick,nzlo_out,nylo_out,nxlo_out);
  }

  memory->destroy(density_fft);
  memory->destroy(greensfn);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(vg);

  memory->destroy1d_offset(fkx,nxlo_fft);
  memory->destroy1d_offset(fky,nylo_fft);
  memory->destroy1d_offset(fkz,nzlo_fft);

  memory->destroy(buf1);
  memory->destroy(buf2);

  memory->destroy(gf_b);
  memory->destroy2d_offset(rho1d,-order/2);
  memory->destroy2d_offset(drho1d,-order/2);
  memory->destroy2d_offset(rho_coeff,(1-order)/2);
  memory->destroy2d_offset(drho_coeff,(1-order)/2);

  delete fft1;
  delete fft2;
  delete remap;
}

/* ----------------------------------------------------------------------
   deallocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPM::deallocate_peratom()
{
  memory->destroy3d_offset(v0_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v1_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v2_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v3_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v4_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v5_brick,nzlo_out,nylo_out,nxlo_out);

  if (differentiation_flag != 1)
    memory->destroy3d_offset(u_brick,nzlo_out,nylo_out,nxlo_out);
  
  memory->destroy(buf3);
  memory->destroy(buf4);
}
  
/* ----------------------------------------------------------------------
   set size of FFT grid (nx,ny,nz_pppm)
------------------------------------------------------------------------- */

void PPPM::set_grid()
{
  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab PPPM
  // 3d PPPM just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;

  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  double h,h_x,h_y,h_z;
  bigint natoms = atom->natoms;

  if (!gewaldflag) {
    if (accuracy <= 0.0)
      error->all(FLERR,"KSpace accuracy must be > 0");
    g_ewald = accuracy*sqrt(natoms*cutoff*xprd*yprd*zprd) / (2.0*q2);
    if (g_ewald >= 1.0) g_ewald = (1.35 - 0.15*log(accuracy))/cutoff;
    else g_ewald = sqrt(-log(g_ewald)) / cutoff;
  }

  // set optimal nx_pppm,ny_pppm,nz_pppm based on order and accuracy
  // nz_pppm uses extended zprd_slab instead of zprd
  // reduce it until accuracy target is met

  if (!gridflag) {
  
    if (differentiation_flag == 1) {

      h = h_x = h_y = h_z = 4.0/g_ewald;
      int count = 0;
      while (1) {
      
        // set grid dimension
        nx_pppm = static_cast<int> (xprd/h_x);
        ny_pppm = static_cast<int> (yprd/h_y);
        nz_pppm = static_cast<int> (zprd_slab/h_z);
      
        if (nx_pppm <= 1) nx_pppm = 2;
        if (ny_pppm <= 1) ny_pppm = 2;
        if (nz_pppm <= 1) nz_pppm = 2;
      
        //set local grid dimension
        int npey_fft,npez_fft;
        if (nz_pppm >= nprocs) {
          npey_fft = 1;
          npez_fft = nprocs;
        } else procs2grid2d(nprocs,ny_pppm,nz_pppm,&npey_fft,&npez_fft);
      
        int me_y = me % npey_fft;
        int me_z = me / npey_fft;
      
        nxlo_fft = 0;
        nxhi_fft = nx_pppm - 1;
        nylo_fft = me_y*ny_pppm/npey_fft;
        nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1;
        nzlo_fft = me_z*nz_pppm/npez_fft;
        nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1;
      
        double df_kspace = compute_df_kspace();
      
        count++;
      
        // break loop if the accuracy has been reached or too many loops have been performed
        if (df_kspace <= accuracy) break;
        if (count > 500) error->all(FLERR, "Could not compute grid size!");
        h *= 0.95;
        h_x = h_y = h_z = h;
      }

    } else {

      double err;
      h_x = h_y = h_z = 1.0/g_ewald;

      nx_pppm = static_cast<int> (xprd/h_x) + 1;
      ny_pppm = static_cast<int> (yprd/h_y) + 1;
      nz_pppm = static_cast<int> (zprd_slab/h_z) + 1;

      err = rms(h_x,xprd,natoms,q2,acons);
      while (err > accuracy) {
        err = rms(h_x,xprd,natoms,q2,acons);
        nx_pppm++;
        h_x = xprd/nx_pppm;
      }

      err = rms(h_y,yprd,natoms,q2,acons);
      while (err > accuracy) {
        err = rms(h_y,yprd,natoms,q2,acons);
        ny_pppm++;
        h_y = yprd/ny_pppm;
      }

      err = rms(h_z,zprd_slab,natoms,q2,acons);
      while (err > accuracy) {
        err = rms(h_z,zprd_slab,natoms,q2,acons);
        nz_pppm++;
        h_z = zprd_slab/nz_pppm;
      }
    }
  }

  // boost grid size until it is factorable

  while (!factorable(nx_pppm)) nx_pppm++;
  while (!factorable(ny_pppm)) ny_pppm++;
  while (!factorable(nz_pppm)) nz_pppm++;
}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int PPPM::factorable(int n)
{
  int i;

  while (n > 1) {
    for (i = 0; i < nfactors; i++) {
      if (n % factors[i] == 0) {
        n /= factors[i];
        break;
      }
    }
    if (i == nfactors) return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   compute estimated kspace force error
------------------------------------------------------------------------- */

double PPPM::compute_df_kspace()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;
  double df_kspace = 0.0;
  if (differentiation_flag == 1) {
    double qopt = compute_qopt();
    df_kspace = sqrt(qopt/natoms)*q2/(xprd*yprd*zprd_slab);
  } else {
    double lprx = rms(xprd/nx_pppm,xprd,natoms,q2,acons);
    double lpry = rms(yprd/ny_pppm,yprd,natoms,q2,acons);
    double lprz = rms(zprd_slab/nz_pppm,zprd_slab,natoms,q2,acons);
    df_kspace = sqrt(lprx*lprx + lpry*lpry + lprz*lprz) / sqrt(3.0);
  }
  return df_kspace;
}

/* ----------------------------------------------------------------------
   compute qopt
------------------------------------------------------------------------- */

double PPPM::compute_qopt()
{
  double qopt = 0.0;
  int i,j,k,l,m,n;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;
  double delvolcell = nx_pppm*ny_pppm*nz_pppm/volume;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double u2, sqk;
  double sum1,sum2,sum3,sum4,dot2;
  double numerator,denominator;

  int nbx = 2;
  int nby = 2;
  int nbz = 2;
  double form = 1.0;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);

        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) +
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
          numerator = form*12.5663706;

          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*pow(qx/g_ewald,2.0));
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*pow(qy/g_ewald,2.0));
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*pow(qz/g_ewald,2.0));
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order);

                dot2 = qx*qx+qy*qy+qz*qz;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += sx*sy*sz*sx*sy*sz/dot2*4.0*4.0*MY_PI*MY_PI;
                sum2 += sx*sy*sz * u2*4.0*MY_PI;
                sum3 += u2;
                sum4 += dot2*u2;
              }
            }
          }
          sum2 *= sum2;
          qopt += sum1 - sum2/(sum3*sum4);
        }
      }
    }
  }
  double qopt_all;
  MPI_Allreduce(&qopt,&qopt_all,1,MPI_DOUBLE,MPI_SUM,world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   estimate kspace force error for ik method
------------------------------------------------------------------------- */

double PPPM::rms(double h, double prd, bigint natoms,
                 double q2_local, double **acons_local)
{
  double sum = 0.0;
  for (int m = 0; m < order; m++)
    sum += acons_local[order][m] * pow(h*g_ewald,2.0*m);
  double value = q2_local * pow(h*g_ewald,(double)order) *
    sqrt(g_ewald*prd*sqrt(2.0*MY_PI)*sum/natoms) / (prd*prd);

  return value;
}

/* ----------------------------------------------------------------------
   compute difference in real-space and KSpace RMS accuracy
------------------------------------------------------------------------- */

double PPPM::diffpr(double h_x, double h_y, double h_z, double q2_local,
                    double **acons_local)
{
  double lprx,lpry,lprz,kspace_prec,real_prec;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  bigint natoms = atom->natoms;

  lprx = rms(h_x,xprd,natoms,q2_local,acons_local);
  lpry = rms(h_y,yprd,natoms,q2_local,acons_local);
  lprz = rms(h_z,zprd*slab_volfactor,natoms,q2_local,acons_local);
  kspace_prec = sqrt(lprx*lprx + lpry*lpry + lprz*lprz) / sqrt(3.0);
  real_prec = 2.0*q2_local * exp(-g_ewald*g_ewald*cutoff*cutoff) /
   sqrt(static_cast<double>(natoms)*cutoff*xprd*yprd*zprd);
  double value = kspace_prec - real_prec;
  return value;
}

/* ----------------------------------------------------------------------
   adjust the g_ewald parameter to near its optimal value
   using a Newton-Raphson solver
------------------------------------------------------------------------- */
void PPPM::adjust_gewald()
{
  double dx;
  
  for (int i = 0; i < LARGE; i++) {
    dx = newton_raphson_f() / derivf();
    g_ewald -= dx; 
    if (fabs(newton_raphson_f()) < SMALL) return;
  }
  
  char str[128];
  sprintf(str, "Could not compute g_ewald");
  error->all(FLERR, str);
}

/* ----------------------------------------------------------------------
 Calculate f(x) using Newton�Raphson solver
 ------------------------------------------------------------------------- */

double PPPM::newton_raphson_f()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;

  double df_rspace = 2.0*q2*exp(-g_ewald*g_ewald*cutoff*cutoff) /
       sqrt(natoms*cutoff*xprd*yprd*zprd);

  double df_kspace = compute_df_kspace();

  return df_rspace - df_kspace;
}

/* ----------------------------------------------------------------------
 Calculate numerical derivative f'(x) using forward difference
 [f(x + h) - f(x)] / h
 ------------------------------------------------------------------------- */

double PPPM::derivf()
{
  double h = 0.000001;  //Derivative step-size
  double df,f1,f2,g_ewald_old;

  f1 = newton_raphson_f();
  g_ewald_old = g_ewald;
  g_ewald += h;
  f2 = newton_raphson_f();
  g_ewald = g_ewald_old;
  df = (f2 - f1)/h;

  return df;
}

/* ----------------------------------------------------------------------
   Calculate the final estimate of the accuracy
------------------------------------------------------------------------- */

double PPPM::final_accuracy()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;
  
  double df_kspace = compute_df_kspace();
  double q2_over_sqrt = q2 / sqrt(natoms*cutoff*xprd*yprd*zprd_slab);
  double df_rspace = 2.0 * q2_over_sqrt * exp(-g_ewald*g_ewald*cutoff*cutoff);
  double df_table = estimate_table_accuracy(q2_over_sqrt,df_rspace);
  double estimated_accuracy = sqrt(df_kspace*df_kspace + df_rspace*df_rspace + 
   df_table*df_table);

  return estimated_accuracy;
}

/* ----------------------------------------------------------------------
   set the FFT parameters
------------------------------------------------------------------------- */

void PPPM::set_fft_parameters()
{
  // global indices of PPPM grid range from 0 to N-1
  // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that I own without ghost cells
  // for slab PPPM, assign z grid as if it were not extended
  
  nxlo_in = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx_pppm);
  nxhi_in = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx_pppm) - 1;

  nylo_in = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny_pppm);
  nyhi_in = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny_pppm) - 1;

  nzlo_in = static_cast<int>
      (comm->zsplit[comm->myloc[2]] * nz_pppm/slab_volfactor);
  nzhi_in = static_cast<int>
      (comm->zsplit[comm->myloc[2]+1] * nz_pppm/slab_volfactor) - 1;      

  // nlower,nupper = stencil size for mapping particles to PPPM grid

  nlower = -(order-1)/2;
  nupper = order/2;

  // shift values for particle <-> grid mapping
  // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

  if (order % 2) shift = OFFSET + 0.5;
  else shift = OFFSET;
  if (order % 2) shiftone = 0.0;
  else shiftone = 0.5;
  
  // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that my particles can contribute charge to
  // effectively nlo_in,nhi_in + ghost cells
  // nlo,nhi = global coords of grid pt to "lower left" of smallest/largest
  //           position a particle in my box can be at
  // dist[3] = particle position bound = subbox + skin/2.0 + qdist
  //   qdist = offset due to TIP4P fictitious charge
  //   convert to triclinic if necessary
  // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping
  // for slab PPPM, assign z grid as if it were not extended

  double *prd,*sublo,*subhi;

  if (triclinic == 0) {
    prd = domain->prd;
    boxlo = domain->boxlo;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }
  
  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double dist[3];
  double cuthalf = 0.5*neighbor->skin + qdist;
  if (triclinic == 0) dist[0] = dist[1] = dist[2] = cuthalf;
  else {
    dist[0] = cuthalf/domain->prd[0];
    dist[1] = cuthalf/domain->prd[1];
    dist[2] = cuthalf/domain->prd[2];
  }
  
  int nlo,nhi;

  nlo = static_cast<int> ((sublo[0]-dist[0]-boxlo[0]) *
                            nx_pppm/xprd + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[0]+dist[0]-boxlo[0]) *
                            nx_pppm/xprd + shift) - OFFSET;
  nxlo_out = nlo + nlower;
  nxhi_out = nhi + nupper;

  nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) *
                            ny_pppm/yprd + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) *
                            ny_pppm/yprd + shift) - OFFSET;
  nylo_out = nlo + nlower;
  nyhi_out = nhi + nupper;

  nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) *
                            nz_pppm/zprd_slab + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) *
                            nz_pppm/zprd_slab + shift) - OFFSET;
  nzlo_out = nlo + nlower;
  nzhi_out = nhi + nupper;

  // for slab PPPM, change the grid boundary for processors at +z end
  //   to include the empty volume between periodically repeating slabs
  // for slab PPPM, want charge data communicated from -z proc to +z proc,
  //   but not vice versa, also want field data communicated from +z proc to
  //   -z proc, but not vice versa
  // this is accomplished by nzhi_in = nzhi_out on +z end (no ghost cells)

  if (slabflag && (comm->myloc[2] == comm->procgrid[2]-1)) {
    nzhi_in = nz_pppm - 1;
    nzhi_out = nz_pppm - 1;
  }

  // nlo_ghost,nhi_ghost = # of planes I will recv from 6 directions
  //   that overlay domain I own
  // proc in that direction tells me via sendrecv()
  // if no neighbor proc, value is from self since I have ghosts regardless
  
  int nplanes;
  MPI_Status status;

  nplanes = nxlo_in - nxlo_out;
  if (comm->procneigh[0][0] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[0][0],0,
                   &nxhi_ghost,1,MPI_INT,comm->procneigh[0][1],0,
                   world,&status);
  else nxhi_ghost = nplanes;

  nplanes = nxhi_out - nxhi_in;
  if (comm->procneigh[0][1] != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[0][1],0,
                   &nxlo_ghost,1,MPI_INT,comm->procneigh[0][0],
                   0,world,&status);
  else nxlo_ghost = nplanes;

  nplanes = nylo_in - nylo_out;
  if (comm->procneigh[1][0] != me)
    MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[1][0],0,
                 &nyhi_ghost,1,MPI_INT,comm->procneigh[1][1],0,
                 world,&status);
  else nyhi_ghost = nplanes;

  nplanes = nyhi_out - nyhi_in;
  if (comm->procneigh[1][1] != me)
    MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[1][1],0,
                 &nylo_ghost,1,MPI_INT,comm->procneigh[1][0],0,
                 world,&status);
  else nylo_ghost = nplanes;

  nplanes = nzlo_in - nzlo_out;
  if (comm->procneigh[2][0] != me)
    MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[2][0],0,
                 &nzhi_ghost,1,MPI_INT,comm->procneigh[2][1],0,
                 world,&status);
  else nzhi_ghost = nplanes;

  nplanes = nzhi_out - nzhi_in;
  if (comm->procneigh[2][1] != me)
    MPI_Sendrecv(&nplanes,1,MPI_INT,comm->procneigh[2][1],0,
                 &nzlo_ghost,1,MPI_INT,comm->procneigh[2][0],0,
                 world,&status);
  else nzlo_ghost = nplanes;
  
  // decomposition of FFT mesh
  // global indices range from 0 to N-1
  // proc owns entire x-dimension, clump of columns in y,z dimensions
  // npey_fft,npez_fft = # of procs in y,z dims
  // if nprocs is small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own

  int npey_fft,npez_fft;
  if (nz_pppm >= nprocs) {
    npey_fft = 1;
    npez_fft = nprocs;
  } else procs2grid2d(nprocs,ny_pppm,nz_pppm,&npey_fft,&npez_fft);

  int me_y = me % npey_fft;
  int me_z = me / npey_fft;

  nxlo_fft = 0;
  nxhi_fft = nx_pppm - 1;
  nylo_fft = me_y*ny_pppm/npey_fft;
  nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1;
  nzlo_fft = me_z*nz_pppm/npez_fft;
  nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1;

  // PPPM grid for this proc, including ghosts
 
  ngrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);

  // FFT arrays on this proc, without ghosts
  // nfft = FFT points in FFT decomposition on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  // nfft_both = greater of 2 values

  nfft = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) *
    (nzhi_fft-nzlo_fft+1);
  int nfft_brick = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) *
    (nzhi_in-nzlo_in+1);
  nfft_both = MAX(nfft,nfft_brick);

  // buffer space for use in brick2fft and fillbrick
  // idel = max # of ghost planes to send or recv in +/- dir of each dim
  // nx,ny,nz = owned planes (including ghosts) in each dim
  // nxx,nyy,nzz = max # of grid cells to send in each dim
  // nbuf = max in any dim

  int idelx,idely,idelz,nx,ny,nz,nxx,nyy,nzz;

  idelx = MAX(nxlo_ghost,nxhi_ghost);
  idelx = MAX(idelx,nxhi_out-nxhi_in);
  idelx = MAX(idelx,nxlo_in-nxlo_out);

  idely = MAX(nylo_ghost,nyhi_ghost);
  idely = MAX(idely,nyhi_out-nyhi_in);
  idely = MAX(idely,nylo_in-nylo_out);

  idelz = MAX(nzlo_ghost,nzhi_ghost);
  idelz = MAX(idelz,nzhi_out-nzhi_in);
  idelz = MAX(idelz,nzlo_in-nzlo_out);

  nx = nxhi_out - nxlo_out + 1;
  ny = nyhi_out - nylo_out + 1;
  nz = nzhi_out - nzlo_out + 1;

  nxx = idelx * ny * nz;
  nyy = idely * nx * nz;
  nzz = idelz * nx * ny;

  nbuf = MAX(nxx,nyy);
  nbuf = MAX(nbuf,nzz);

  nbuf_peratom = 6*nbuf;
  if (differentiation_flag != 1) {
    nbuf *= 3;
    nbuf_peratom++;
  }
}

/* ----------------------------------------------------------------------
   pre-compute Green's function denominator expansion coeffs, Gamma(2n)
------------------------------------------------------------------------- */

void PPPM::compute_gf_denom()
{
  int k,l,m;

  for (l = 1; l < order; l++) gf_b[l] = 0.0;
  gf_b[0] = 1.0;

  for (m = 1; m < order; m++) {
    for (l = m; l > 0; l--)
      gf_b[l] = 4.0 * (gf_b[l]*(l-m)*(l-m-0.5)-gf_b[l-1]*(l-m-1)*(l-m-1));
    gf_b[0] = 4.0 * (gf_b[0]*(l-m)*(l-m-0.5));
  }

  int ifact = 1;
  for (k = 1; k < 2*order; k++) ifact *= k;
  double gaminv = 1.0/ifact;
  for (l = 0; l < order; l++) gf_b[l] *= gaminv;
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
------------------------------------------------------------------------- */

void PPPM::compute_gf_ik()
{
  int nx,ny,nz,kper,lper,mper;
  double snx,sny,snz,snx2,sny2,snz2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,dot1,dot2;
  double numerator,denominator;
  double sqk;

  int k,l,m,n;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;
  
  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);
  
  int nbx = static_cast<int> ((g_ewald*xprd/(MY_PI*nx_pppm)) *
                              pow(-log(EPS_HOC),0.25));
  int nby = static_cast<int> ((g_ewald*yprd/(MY_PI*ny_pppm)) *
                              pow(-log(EPS_HOC),0.25));
  int nbz = static_cast<int> ((g_ewald*zprd_slab/(MY_PI*nz_pppm)) *
                              pow(-log(EPS_HOC),0.25));

  double form = 1.0;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    snz = sin(0.5*unitkz*mper*zprd_slab/nz_pppm);
    snz2 = snz*snz;

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      sny = sin(0.5*unitky*lper*yprd/ny_pppm);
      sny2 = sny*sny;

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        snx = sin(0.5*unitkx*kper*xprd/nx_pppm);
        snx2 = snx*snx;

        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) +
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
          numerator = form*12.5663706/sqk;
          denominator = gf_denom(snx2,sny2,snz2);
          sum1 = 0.0;
          const double dorder = static_cast<double>(order);
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*pow(qx/g_ewald,2.0));
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm;
            if (argx != 0.0) wx = pow(sin(argx)/argx,dorder);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*pow(qy/g_ewald,2.0));
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm;
              if (argy != 0.0) wy = pow(sin(argy)/argy,dorder);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*pow(qz/g_ewald,2.0));
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm;
                if (argz != 0.0) wz = pow(sin(argz)/argz,dorder);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx+qy*qy+qz*qz;
                sum1 += (dot1/dot2) * sx*sy*sz * pow(wx*wy*wz,2.0);
              }
            }
          }
          greensfn[n++] = numerator*sum1/denominator;
        } else greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute optimized Green's function for energy calculation
------------------------------------------------------------------------- */

void PPPM::compute_gf_ad()
{
  int i,j,k,l,m,n;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double snx,sny,snz,snx2,sny2,snz2;
  double sqk, u2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double numerator,denominator;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    qz = unitkz*mper;
    snz = sin(0.5*qz*zprd_slab/nz_pppm);
    snz2 = snz*snz;
    sz = exp(-0.25*pow(qz/g_ewald,2.0));
    wz = 1.0;
    argz = 0.5*qz*zprd_slab/nz_pppm;
    if (argz != 0.0) wz = pow(sin(argz)/argz,order);
    wz *= wz;

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      qy = unitky*lper;
      sny = sin(0.5*qy*yprd/ny_pppm);
      sny2 = sny*sny;
      sy = exp(-0.25*pow(qy/g_ewald,2.0));
      wy = 1.0;
      argy = 0.5*qy*yprd/ny_pppm;
      if (argy != 0.0) wy = pow(sin(argy)/argy,order);
      wy *= wy;

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        qx = unitkx*kper;
        snx = sin(0.5*qx*xprd/nx_pppm);
        snx2 = snx*snx;
        sx = exp(-0.25*pow(qx/g_ewald,2.0));
        wx = 1.0;
        argx = 0.5*qx*xprd/nx_pppm;
        if (argx != 0.0) wx = pow(sin(argx)/argx,order);
        wx *= wx;

        sqk = pow(qx,2.0) + pow(qy,2.0) + pow(qz,2.0);

        if (sqk != 0.0) {
          numerator = 4.0*MY_PI/sqk;
          denominator = gf_denom(snx2,sny2,snz2);
          greensfn[n++] = numerator*sx*sy*sz*wx*wy*wz/denominator;
        } else greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute self force coefficients for ad-differentiation scheme
------------------------------------------------------------------------- */

void PPPM::compute_sf_coeff()
{

  int i,j,k,l,m,n;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double argx,argy,argz;
  double wx0[5],wy0[5],wz0[5],wx1[5],wy1[5],wz1[5],wx2[5],wy2[5],wz2[5];
  double qx0,qy0,qz0,qx1,qy1,qz1,qx2,qy2,qz2;
  double u0,u1,u2,u3,u4,u5,u6;
  double sum1,sum2,sum3,sum4,sum5,sum6;

  int nb = 2;

  double form = 1.0;
  for (n = 0; n < 6; n++) sf_coeff[n] = 0.0;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);

        sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = 0.0;
        for (i = -nb; i <= nb; i++) {

          qx0 = unitkx*(kper+nx_pppm*i);
          qx1 = unitkx*(kper+nx_pppm*(i+1));
          qx2 = unitkx*(kper+nx_pppm*(i+2));
          wx0[i+2] = 1.0;
          wx1[i+2] = 1.0;
          wx2[i+2] = 1.0;
          argx = 0.5*qx0*xprd/nx_pppm;
          if (argx != 0.0) wx0[i+2] = pow(sin(argx)/argx,order);
          argx = 0.5*qx1*xprd/nx_pppm;
          if (argx != 0.0) wx1[i+2] = pow(sin(argx)/argx,order);
          argx = 0.5*qx2*xprd/nx_pppm;
          if (argx != 0.0) wx2[i+2] = pow(sin(argx)/argx,order);

          qy0 = unitky*(lper+ny_pppm*i);
          qy1 = unitky*(lper+ny_pppm*(i+1));
          qy2 = unitky*(lper+ny_pppm*(i+2));
          wy0[i+2] = 1.0;
          wy1[i+2] = 1.0;
          wy2[i+2] = 1.0;
          argy = 0.5*qy0*yprd/ny_pppm;
          if (argy != 0.0) wy0[i+2] = pow(sin(argy)/argy,order);
          argy = 0.5*qy1*yprd/ny_pppm;
          if (argy != 0.0) wy1[i+2] = pow(sin(argy)/argy,order);
          argy = 0.5*qy2*yprd/ny_pppm;
          if (argy != 0.0) wy2[i+2] = pow(sin(argy)/argy,order);

          qz0 = unitkz*(mper+nz_pppm*i);
          qz1 = unitkz*(mper+nz_pppm*(i+1));
          qz2 = unitkz*(mper+nz_pppm*(i+2));
          wz0[i+2] = 1.0;
          wz1[i+2] = 1.0;
          wz2[i+2] = 1.0;
          argz = 0.5*qz0*zprd_slab/nz_pppm;
          if (argz != 0.0) wz0[i+2] = pow(sin(argz)/argz,order);
          argz = 0.5*qz1*zprd_slab/nz_pppm;
          if (argz != 0.0) wz1[i+2] = pow(sin(argz)/argz,order);
           argz = 0.5*qz2*zprd_slab/nz_pppm;
          if (argz != 0.0) wz2[i+2] = pow(sin(argz)/argz,order);
        }

        for (nx = 0; nx <= 4; nx++) {
          for (ny = 0; ny <= 4; ny++) {
            for (nz = 0; nz <= 4; nz++) {
              u0 = wx0[nx]*wy0[ny]*wz0[nz];
              u1 = wx1[nx]*wy0[ny]*wz0[nz];
              u2 = wx2[nx]*wy0[ny]*wz0[nz];
              u3 = wx0[nx]*wy1[ny]*wz0[nz];
              u4 = wx0[nx]*wy2[ny]*wz0[nz];
              u5 = wx0[nx]*wy0[ny]*wz1[nz];
              u6 = wx0[nx]*wy0[ny]*wz2[nz];

              sum1 += u0*u1;
              sum2 += u0*u2;
              sum3 += u0*u3;
              sum4 += u0*u4;
              sum5 += u0*u5;
              sum6 += u0*u6;
            }
          }
        }

        // multiplication with Green's function

        sf_coeff[0] += sum1 * greensfn[n];
        sf_coeff[1] += sum2 * greensfn[n];
        sf_coeff[2] += sum3 * greensfn[n];
        sf_coeff[3] += sum4 * greensfn[n];
        sf_coeff[4] += sum5 * greensfn[n];
        sf_coeff[5] += sum6 * greensfn[n++];
      }
    }
  }

  // perform multiplication with prefactors

  double prex, prey, prez;
  prex = prey = prez = MY_PI/volume;
  prex *= nx_pppm/xprd;
  prey *= ny_pppm/yprd;
  prez *= nz_pppm/zprd_slab;
  sf_coeff[0] *= prex;
  sf_coeff[1] *= prex*2;
  sf_coeff[2] *= prey;
  sf_coeff[3] *= prey*2;
  sf_coeff[4] *= prez;
  sf_coeff[5] *= prez*2;

  // communicate values with other procs

  double tmp[6];
  MPI_Allreduce(sf_coeff,tmp,6,MPI_DOUBLE,MPI_SUM,world);
  for (n = 0; n < 6; n++) sf_coeff[n] = tmp[n];
}

/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

void PPPM::brick2fft()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my ghosts for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in+1; ix <= nxhi_out; ix++)
        buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix < nxlo_in+nxlo_ghost; ix++)
        density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix < nxlo_in; ix++)
        buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in-nxhi_ghost+1; ix <= nxhi_in; ix++)
        density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        density_brick[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        buf1[n++] = density_brick[iz][iy][ix];

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        density_brick[iz][iy][ix] += buf2[n++];

  // remap from 3d brick decomposition to FFT decomposition
  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        density_fft[n++] = density_brick[iz][iy][ix];

  remap->perform(density_fft,density_fft,work1);
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with field values
------------------------------------------------------------------------- */

void PPPM::fillbrick()
{
  if (differentiation_flag == 1) fillbrick_ad();
  else fillbrick_ik();
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with field values for ik
------------------------------------------------------------------------- */

void PPPM::fillbrick_ik()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my real cells for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf1[n++] = vdx_brick[iz][iy][ix];
        buf1[n++] = vdy_brick[iz][iy][ix];
        buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        vdx_brick[iz][iy][ix] = buf2[n++];
        vdy_brick[iz][iy][ix] = buf2[n++];
        vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf1[n++] = vdx_brick[iz][iy][ix];
        buf1[n++] = vdy_brick[iz][iy][ix];
        buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        vdx_brick[iz][iy][ix] = buf2[n++];
        vdy_brick[iz][iy][ix] = buf2[n++];
        vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf1[n++] = vdx_brick[iz][iy][ix];
        buf1[n++] = vdy_brick[iz][iy][ix];
        buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        vdx_brick[iz][iy][ix] = buf2[n++];
        vdy_brick[iz][iy][ix] = buf2[n++];
        vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf1[n++] = vdx_brick[iz][iy][ix];
        buf1[n++] = vdy_brick[iz][iy][ix];
        buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        vdx_brick[iz][iy][ix] = buf2[n++];
        vdy_brick[iz][iy][ix] = buf2[n++];
        vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in-nxhi_ghost+1; ix <= nxhi_in; ix++) {
        buf1[n++] = vdx_brick[iz][iy][ix];
        buf1[n++] = vdy_brick[iz][iy][ix];
        buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix < nxlo_in; ix++) {
        vdx_brick[iz][iy][ix] = buf2[n++];
        vdy_brick[iz][iy][ix] = buf2[n++];
        vdz_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix < nxlo_in+nxlo_ghost; ix++) {
        buf1[n++] = vdx_brick[iz][iy][ix];
        buf1[n++] = vdy_brick[iz][iy][ix];
        buf1[n++] = vdz_brick[iz][iy][ix];
      }

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in+1; ix <= nxhi_out; ix++) {
        vdx_brick[iz][iy][ix] = buf2[n++];
        vdy_brick[iz][iy][ix] = buf2[n++];
        vdz_brick[iz][iy][ix] = buf2[n++];
      }
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with field values for ad
------------------------------------------------------------------------- */

void PPPM::fillbrick_ad()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my real cells for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf1[n++] = u_brick[iz][iy][ix];
      }

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        u_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf1[n++] = u_brick[iz][iy][ix];
      }

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        u_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf1[n++] = u_brick[iz][iy][ix];
      }

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        u_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf1[n++] = u_brick[iz][iy][ix];
      }

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        u_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in-nxhi_ghost+1; ix <= nxhi_in; ix++) {
        buf1[n++] = u_brick[iz][iy][ix];
      }

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix < nxlo_in; ix++) {
        u_brick[iz][iy][ix] = buf2[n++];
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix < nxlo_in+nxlo_ghost; ix++) {
        buf1[n++] = u_brick[iz][iy][ix];
      }

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in+1; ix <= nxhi_out; ix++) {
        u_brick[iz][iy][ix] = buf2[n++];
      }
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with per-atom field values
------------------------------------------------------------------------- */

void PPPM::fillbrick_peratom()
{
  if (differentiation_flag == 1) fillbrick_peratom_ad();
  else fillbrick_peratom_ik();
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with per-atom field values for ik
------------------------------------------------------------------------- */

void PPPM::fillbrick_peratom_ik()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my real cells for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (eflag_atom) buf3[n++] = u_brick[iz][iy][ix];
        if (vflag_atom) {
          buf3[n++] = v0_brick[iz][iy][ix];
          buf3[n++] = v1_brick[iz][iy][ix];
          buf3[n++] = v2_brick[iz][iy][ix];
          buf3[n++] = v3_brick[iz][iy][ix];
          buf3[n++] = v4_brick[iz][iy][ix];
          buf3[n++] = v5_brick[iz][iy][ix];
        }
      }

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (eflag_atom) u_brick[iz][iy][ix] = buf4[n++];
        if (vflag_atom) {
          v0_brick[iz][iy][ix] = buf4[n++];
          v1_brick[iz][iy][ix] = buf4[n++];
          v2_brick[iz][iy][ix] = buf4[n++];
          v3_brick[iz][iy][ix] = buf4[n++];
          v4_brick[iz][iy][ix] = buf4[n++];
          v5_brick[iz][iy][ix] = buf4[n++];
        }
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (eflag_atom) buf3[n++] = u_brick[iz][iy][ix];
        if (vflag_atom) {
          buf3[n++] = v0_brick[iz][iy][ix];
          buf3[n++] = v1_brick[iz][iy][ix];
          buf3[n++] = v2_brick[iz][iy][ix];
          buf3[n++] = v3_brick[iz][iy][ix];
          buf3[n++] = v4_brick[iz][iy][ix];
          buf3[n++] = v5_brick[iz][iy][ix];
        }
      }

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (eflag_atom) u_brick[iz][iy][ix] = buf4[n++];
        if (vflag_atom) {
          v0_brick[iz][iy][ix] = buf4[n++];
          v1_brick[iz][iy][ix] = buf4[n++];
          v2_brick[iz][iy][ix] = buf4[n++];
          v3_brick[iz][iy][ix] = buf4[n++];
          v4_brick[iz][iy][ix] = buf4[n++];
          v5_brick[iz][iy][ix] = buf4[n++];
        }
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (eflag_atom) buf3[n++] = u_brick[iz][iy][ix];
        if (vflag_atom) {
          buf3[n++] = v0_brick[iz][iy][ix];
          buf3[n++] = v1_brick[iz][iy][ix];
          buf3[n++] = v2_brick[iz][iy][ix];
          buf3[n++] = v3_brick[iz][iy][ix];
          buf3[n++] = v4_brick[iz][iy][ix];
          buf3[n++] = v5_brick[iz][iy][ix];
        }
      }

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (eflag_atom) u_brick[iz][iy][ix] = buf4[n++];
        if (vflag_atom) {
          v0_brick[iz][iy][ix] = buf4[n++];
          v1_brick[iz][iy][ix] = buf4[n++];
          v2_brick[iz][iy][ix] = buf4[n++];
          v3_brick[iz][iy][ix] = buf4[n++];
          v4_brick[iz][iy][ix] = buf4[n++];
          v5_brick[iz][iy][ix] = buf4[n++];
        }
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (eflag_atom) buf3[n++] = u_brick[iz][iy][ix];
        if (vflag_atom) {
          buf3[n++] = v0_brick[iz][iy][ix];
          buf3[n++] = v1_brick[iz][iy][ix];
          buf3[n++] = v2_brick[iz][iy][ix];
          buf3[n++] = v3_brick[iz][iy][ix];
          buf3[n++] = v4_brick[iz][iy][ix];
          buf3[n++] = v5_brick[iz][iy][ix];
        }
      }

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (eflag_atom) u_brick[iz][iy][ix] = buf4[n++];
        if (vflag_atom) {
          v0_brick[iz][iy][ix] = buf4[n++];
          v1_brick[iz][iy][ix] = buf4[n++];
          v2_brick[iz][iy][ix] = buf4[n++];
          v3_brick[iz][iy][ix] = buf4[n++];
          v4_brick[iz][iy][ix] = buf4[n++];
          v5_brick[iz][iy][ix] = buf4[n++];
        }
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in-nxhi_ghost+1; ix <= nxhi_in; ix++) {
        if (eflag_atom) buf3[n++] = u_brick[iz][iy][ix];
        if (vflag_atom) {
          buf3[n++] = v0_brick[iz][iy][ix];
          buf3[n++] = v1_brick[iz][iy][ix];
          buf3[n++] = v2_brick[iz][iy][ix];
          buf3[n++] = v3_brick[iz][iy][ix];
          buf3[n++] = v4_brick[iz][iy][ix];
          buf3[n++] = v5_brick[iz][iy][ix];
        }
      }

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix < nxlo_in; ix++) {
        if (eflag_atom) u_brick[iz][iy][ix] = buf4[n++];
        if (vflag_atom) {
          v0_brick[iz][iy][ix] = buf4[n++];
          v1_brick[iz][iy][ix] = buf4[n++];
          v2_brick[iz][iy][ix] = buf4[n++];
          v3_brick[iz][iy][ix] = buf4[n++];
          v4_brick[iz][iy][ix] = buf4[n++];
          v5_brick[iz][iy][ix] = buf4[n++];
        }
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix < nxlo_in+nxlo_ghost; ix++) {
        if (eflag_atom) buf3[n++] = u_brick[iz][iy][ix];
        if (vflag_atom) {
          buf3[n++] = v0_brick[iz][iy][ix];
          buf3[n++] = v1_brick[iz][iy][ix];
          buf3[n++] = v2_brick[iz][iy][ix];
          buf3[n++] = v3_brick[iz][iy][ix];
          buf3[n++] = v4_brick[iz][iy][ix];
          buf3[n++] = v5_brick[iz][iy][ix];
        }
      }

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in+1; ix <= nxhi_out; ix++) {
        if (eflag_atom) u_brick[iz][iy][ix] = buf4[n++];
        if (vflag_atom) {
          v0_brick[iz][iy][ix] = buf4[n++];
          v1_brick[iz][iy][ix] = buf4[n++];
          v2_brick[iz][iy][ix] = buf4[n++];
          v3_brick[iz][iy][ix] = buf4[n++];
          v4_brick[iz][iy][ix] = buf4[n++];
          v5_brick[iz][iy][ix] = buf4[n++];
        }
      }
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with per-atom field values for ad
------------------------------------------------------------------------- */

void PPPM::fillbrick_peratom_ad()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my real cells for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf3[n++] = v0_brick[iz][iy][ix];
        buf3[n++] = v1_brick[iz][iy][ix];
        buf3[n++] = v2_brick[iz][iy][ix];
        buf3[n++] = v3_brick[iz][iy][ix];
        buf3[n++] = v4_brick[iz][iy][ix];
        buf3[n++] = v5_brick[iz][iy][ix];
      }

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        v0_brick[iz][iy][ix] = buf4[n++];
        v1_brick[iz][iy][ix] = buf4[n++];
        v2_brick[iz][iy][ix] = buf4[n++];
        v3_brick[iz][iy][ix] = buf4[n++];
        v4_brick[iz][iy][ix] = buf4[n++];
        v5_brick[iz][iy][ix] = buf4[n++];
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf3[n++] = v0_brick[iz][iy][ix];
        buf3[n++] = v1_brick[iz][iy][ix];
        buf3[n++] = v2_brick[iz][iy][ix];
        buf3[n++] = v3_brick[iz][iy][ix];
        buf3[n++] = v4_brick[iz][iy][ix];
        buf3[n++] = v5_brick[iz][iy][ix];
      }

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        v0_brick[iz][iy][ix] = buf4[n++];
        v1_brick[iz][iy][ix] = buf4[n++];
        v2_brick[iz][iy][ix] = buf4[n++];
        v3_brick[iz][iy][ix] = buf4[n++];
        v4_brick[iz][iy][ix] = buf4[n++];
        v5_brick[iz][iy][ix] = buf4[n++];
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf3[n++] = v0_brick[iz][iy][ix];
        buf3[n++] = v1_brick[iz][iy][ix];
        buf3[n++] = v2_brick[iz][iy][ix];
        buf3[n++] = v3_brick[iz][iy][ix];
        buf3[n++] = v4_brick[iz][iy][ix];
        buf3[n++] = v5_brick[iz][iy][ix];
      }

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        v0_brick[iz][iy][ix] = buf4[n++];
        v1_brick[iz][iy][ix] = buf4[n++];
        v2_brick[iz][iy][ix] = buf4[n++];
        v3_brick[iz][iy][ix] = buf4[n++];
        v4_brick[iz][iy][ix] = buf4[n++];
        v5_brick[iz][iy][ix] = buf4[n++];
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf3[n++] = v0_brick[iz][iy][ix];
        buf3[n++] = v1_brick[iz][iy][ix];
        buf3[n++] = v2_brick[iz][iy][ix];
        buf3[n++] = v3_brick[iz][iy][ix];
        buf3[n++] = v4_brick[iz][iy][ix];
        buf3[n++] = v5_brick[iz][iy][ix];
      }

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        v0_brick[iz][iy][ix] = buf4[n++];
        v1_brick[iz][iy][ix] = buf4[n++];
        v2_brick[iz][iy][ix] = buf4[n++];
        v3_brick[iz][iy][ix] = buf4[n++];
        v4_brick[iz][iy][ix] = buf4[n++];
        v5_brick[iz][iy][ix] = buf4[n++];
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in-nxhi_ghost+1; ix <= nxhi_in; ix++) {
        buf3[n++] = v0_brick[iz][iy][ix];
        buf3[n++] = v1_brick[iz][iy][ix];
        buf3[n++] = v2_brick[iz][iy][ix];
        buf3[n++] = v3_brick[iz][iy][ix];
        buf3[n++] = v4_brick[iz][iy][ix];
        buf3[n++] = v5_brick[iz][iy][ix];
      }

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix < nxlo_in; ix++) {
        v0_brick[iz][iy][ix] = buf4[n++];
        v1_brick[iz][iy][ix] = buf4[n++];
        v2_brick[iz][iy][ix] = buf4[n++];
        v3_brick[iz][iy][ix] = buf4[n++];
        v4_brick[iz][iy][ix] = buf4[n++];
        v5_brick[iz][iy][ix] = buf4[n++];
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix < nxlo_in+nxlo_ghost; ix++) {
        buf3[n++] = v0_brick[iz][iy][ix];
        buf3[n++] = v1_brick[iz][iy][ix];
        buf3[n++] = v2_brick[iz][iy][ix];
        buf3[n++] = v3_brick[iz][iy][ix];
        buf3[n++] = v4_brick[iz][iy][ix];
        buf3[n++] = v5_brick[iz][iy][ix];
      }

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf4[i] = buf3[i];
  else {
    MPI_Irecv(buf4,nbuf_peratom,MPI_FFT_SCALAR,
              comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf3,n,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in+1; ix <= nxhi_out; ix++) {
        v0_brick[iz][iy][ix] = buf4[n++];
        v1_brick[iz][iy][ix] = buf4[n++];
        v2_brick[iz][iy][ix] = buf4[n++];
        v3_brick[iz][iy][ix] = buf4[n++];
        v4_brick[iz][iy][ix] = buf4[n++];
        v5_brick[iz][iy][ix] = buf4[n++];
      }
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

void PPPM::particle_map()
{
  int nx,ny,nz;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift) - OFFSET;
    ny = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift) - OFFSET;
    nz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
        nz+nlower < nzlo_out || nz+nupper > nzhi_out)
      flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute PPPM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPM::make_rho()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  memset(&(density_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          density_brick[mz][my][mx] += x0*rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver
------------------------------------------------------------------------- */

void PPPM::poisson()
{
  if (differentiation_flag == 1) poisson_ad();
  else poisson_ik();
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ik
------------------------------------------------------------------------- */

void PPPM::poisson_ik()
{
  int i,j,k,n;
  double eng;

  // transform charge density (r -> k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] = density_fft[i];
    work1[n++] = ZEROF;
  }

  fft1->compute(work1,work1,1);

  // global energy and virial contribution

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft; i++) {
        eng = s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        for (j = 0; j < 6; j++) virial[j] += eng*vg[i][j];
        if (eflag_global) energy += eng;
        n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft; i++) {
        energy +=
          s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        n += 2;
      }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] *= scaleinv * greensfn[i];
    work1[n++] *= scaleinv * greensfn[i];
  }

  // extra FFTs for per-atom energy/virial

  if (evflag_atom) poisson_peratom();

  // compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = fkx[i]*work1[n+1];
        work2[n+1] = -fkx[i]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdx_brick[k][j][i] = work2[n];
        n += 2;
      }

  // y direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = fky[j]*work1[n+1];
        work2[n+1] = -fky[j]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdy_brick[k][j][i] = work2[n];
        n += 2;
      }

  // z direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = fkz[k]*work1[n+1];
        work2[n+1] = -fkz[k]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdz_brick[k][j][i] = work2[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ad
------------------------------------------------------------------------- */

void PPPM::poisson_ad()
{
  int i,j,k,n;
  double eng;

  // transform charge density (r -> k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] = density_fft[i];
    work1[n++] = ZEROF;
  }

  fft1->compute(work1,work1,1);

  // global energy and virial contribution

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft; i++) {
        eng = s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        for (j = 0; j < 6; j++) virial[j] += eng*vg[i][j];
        if (eflag_global) energy += eng;
        n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft; i++) {
        energy +=
          s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        n += 2;
      }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] *= scaleinv * greensfn[i];
    work1[n++] *= scaleinv * greensfn[i];
  }

  // extra FFTs for per-atom energy/virial

  if (vflag_atom) poisson_peratom();

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n];
    work2[n+1] = work1[n+1];
    n += 2;
  }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        u_brick[k][j][i] = work2[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for per-atom energy/virial
------------------------------------------------------------------------- */

void PPPM::poisson_peratom()
{
  int i,j,k,n;

  // energy

  if (eflag_atom && differentiation_flag != 1) {
    n = 0;
    for (i = 0; i < nfft; i++) {
      work2[n] = work1[n];
      work2[n+1] = work1[n+1];
      n += 2;
    }

    fft2->compute(work2,work2,-1);

    n = 0;
    for (k = nzlo_in; k <= nzhi_in; k++)
      for (j = nylo_in; j <= nyhi_in; j++)
        for (i = nxlo_in; i <= nxhi_in; i++) {
          u_brick[k][j][i] = work2[n];
          n += 2;
        }
  }

  // 6 components of virial in v0 thru v5

  if (!vflag_atom) return;

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][0];
    work2[n+1] = work1[n+1]*vg[i][0];
    n += 2;
  }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v0_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][1];
    work2[n+1] = work1[n+1]*vg[i][1];
    n += 2;
  }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v1_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][2];
    work2[n+1] = work1[n+1]*vg[i][2];
    n += 2;
  }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v2_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][3];
    work2[n+1] = work1[n+1]*vg[i][3];
    n += 2;
  }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v3_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][4];
    work2[n+1] = work1[n+1]*vg[i][4];
    n += 2;
  }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v4_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][5];
    work2[n+1] = work1[n+1]*vg[i][5];
    n += 2;
  }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v5_brick[k][j][i] = work2[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
------------------------------------------------------------------------- */

void PPPM::fieldforce()
{
  if (differentiation_flag == 1) fieldforce_ad();
  else fieldforce_ik();
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPM::fieldforce_ik()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          ekx -= x0*vdx_brick[mz][my][mx];
          eky -= x0*vdy_brick[mz][my][mx];
          ekz -= x0*vdz_brick[mz][my][mx];
        }
      }
    }

    // convert E-field to force

    const double qfactor = force->qqrd2e * scale * q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    if (slabflag != 2) f[i][2] += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void PPPM::fieldforce_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0,dx0,dy0,dz0;
  FFT_SCALAR ekx,eky,ekz;
  double s1,s2,s3;
  double sf = 0.0;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double hx_inv = nx_pppm/xprd;
  double hy_inv = ny_pppm/yprd;
  double hz_inv = nz_pppm/zprd;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);
    compute_drho1d(dx,dy,dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          ekx += drho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          eky += rho1d[0][l]*drho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekz += rho1d[0][l]*rho1d[1][m]*drho1d[2][n]*u_brick[mz][my][mx];
        }
      }
    }
    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;
    // convert E-field to force and substract self forces
    const double qfactor = force->qqrd2e * scale;

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;
    sf = sf_coeff[0]*sin(2*MY_PI*s1);
    sf += sf_coeff[1]*sin(4*MY_PI*s1);
    sf *= 2*q[i]*q[i];
    f[i][0] += qfactor*(ekx*q[i] - sf);

    sf = sf_coeff[2]*sin(2*MY_PI*s2);
    sf += sf_coeff[3]*sin(4*MY_PI*s2);
    sf *= 2*q[i]*q[i];
    f[i][1] += qfactor*(eky*q[i] - sf);


    sf = sf_coeff[4]*sin(2*MY_PI*s3);
    sf += sf_coeff[5]*sin(4*MY_PI*s3);
    sf *= 2*q[i]*q[i];
    if (slabflag != 2) f[i][2] += qfactor*(ekz*q[i] - sf);
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

void PPPM::fieldforce_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u,v0,v1,v2,v3,v4,v5;

  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    u = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          if (eflag_atom) u += x0*u_brick[mz][my][mx];
          if (vflag_atom) {
            v0 += x0*v0_brick[mz][my][mx];
            v1 += x0*v1_brick[mz][my][mx];
            v2 += x0*v2_brick[mz][my][mx];
            v3 += x0*v3_brick[mz][my][mx];
            v4 += x0*v4_brick[mz][my][mx];
            v5 += x0*v5_brick[mz][my][mx];
          }
        }
      }
    }

    if (eflag_atom) eatom[i] += q[i]*u;
    if (vflag_atom) {
      vatom[i][0] += q[i]*v0;
      vatom[i][1] += q[i]*v1;
      vatom[i][2] += q[i]*v2;
      vatom[i][3] += q[i]*v3;
      vatom[i][4] += q[i]*v4;
      vatom[i][5] += q[i]*v5;
    }
  }
}

/* ----------------------------------------------------------------------
   map nprocs to NX by NY grid as PX by PY procs - return optimal px,py
------------------------------------------------------------------------- */

void PPPM::procs2grid2d(int nprocs, int nx, int ny, int *px, int *py)
{
  // loop thru all possible factorizations of nprocs
  // surf = surface area of largest proc sub-domain
  // innermost if test minimizes surface area and surface/volume ratio

  int bestsurf = 2 * (nx + ny);
  int bestboxx = 0;
  int bestboxy = 0;

  int boxx,boxy,surf,ipx,ipy;

  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      ipy = nprocs/ipx;
      boxx = nx/ipx;
      if (nx % ipx) boxx++;
      boxy = ny/ipy;
      if (ny % ipy) boxy++;
      surf = boxx + boxy;
      if (surf < bestsurf ||
          (surf == bestsurf && boxx*boxy > bestboxx*bestboxy)) {
        bestsurf = surf;
        bestboxx = boxx;
        bestboxy = boxy;
        *px = ipx;
        *py = ipy;
      }
    }
    ipx++;
  }
}

/* ----------------------------------------------------------------------
   charge assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */

void PPPM::compute_rho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                         const FFT_SCALAR &dz)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = order-1; l >= 0; l--) {
      r1 = rho_coeff[l][k] + r1*dx;
      r2 = rho_coeff[l][k] + r2*dy;
      r3 = rho_coeff[l][k] + r3*dz;
    }
    rho1d[0][k] = r1;
    rho1d[1][k] = r2;
    rho1d[2][k] = r3;
  }
}

/* ----------------------------------------------------------------------
   charge assignment into drho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */

void PPPM::compute_drho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                          const FFT_SCALAR &dz)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = order-2; l >= 0; l--) {
      r1 = drho_coeff[l][k] + r1*dx;
      r2 = drho_coeff[l][k] + r2*dy;
      r3 = drho_coeff[l][k] + r3*dz;
    }
    drho1d[0][k] = r1;
    drho1d[1][k] = r2;
    drho1d[2][k] = r3;
  }
}

/* ----------------------------------------------------------------------
   generate coeffients for the weight function of order n

              (n-1)
  Wn(x) =     Sum    wn(k,x) , Sum is over every other integer
           k=-(n-1)
  For k=-(n-1),-(n-1)+2, ....., (n-1)-2,n-1
      k is odd integers if n is even and even integers if n is odd
              ---
             | n-1
             | Sum a(l,j)*(x-k/2)**l   if abs(x-k/2) < 1/2
  wn(k,x) = <  l=0
             |
             |  0                       otherwise
              ---
  a coeffients are packed into the array rho_coeff to eliminate zeros
  rho_coeff(l,((k+mod(n+1,2))/2) = a(l,k)
------------------------------------------------------------------------- */

void PPPM::compute_rho_coeff()
{
  int j,k,l,m;
  FFT_SCALAR s;

  FFT_SCALAR **a;
  memory->create2d_offset(a,order,-order,order,"pppm:a");

  for (k = -order; k <= order; k++)
    for (l = 0; l < order; l++)
      a[l][k] = 0.0;

  a[0][0] = 1.0;
  for (j = 1; j < order; j++) {
    for (k = -j; k <= j; k += 2) {
      s = 0.0;
      for (l = 0; l < j; l++) {
        a[l+1][k] = (a[l][k+1]-a[l][k-1]) / (l+1);
#ifdef FFT_SINGLE
        s += powf(0.5,(float) l+1) *
          (a[l][k-1] + powf(-1.0,(float) l) * a[l][k+1]) / (l+1);
#else
        s += pow(0.5,(double) l+1) *
          (a[l][k-1] + pow(-1.0,(double) l) * a[l][k+1]) / (l+1);
#endif
      }
      a[0][k] = s;
    }
  }

  m = (1-order)/2;
  for (k = -(order-1); k < order; k += 2) {
    for (l = 0; l < order; l++)
      rho_coeff[l][m] = a[l][k];
    for (l = 1; l < order; l++)
      drho_coeff[l-1][m] = l*a[l][k];
    m++;
  }

  memory->destroy2d_offset(a,-order);
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane.
------------------------------------------------------------------------- */

void PPPM::slabcorr()
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // compute corrections

  const double e_slabcorr = 2.0*MY_PI*dipole_all*dipole_all/volume;
  const double qscale = force->qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = 2.0*MY_PI*dipole_all/volume;
    for (int i = 0; i < nlocal; i++) eatom[i] += qscale * q[i]*x[i][2]*efact;
  }

  // add on force corrections

  double ffact = -4.0*MY_PI*dipole_all/volume;
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) f[i][2] += qscale * q[i]*ffact;
}

/* ----------------------------------------------------------------------
   perform and time the FFTs required for N timesteps
------------------------------------------------------------------------- */

void PPPM::timing(int n, double &time3d, double &time1d)
{
  double time1,time2;

  for (int i = 0; i < 2*nfft_both; i++) work1[i] = ZEROF;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  for (int i = 0; i < n; i++) {
    fft1->compute(work1,work1,1);
    fft2->compute(work1,work1,-1);
    if (differentiation_flag != 1) {
      fft2->compute(work1,work1,-1);
      fft2->compute(work1,work1,-1);
    }
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time3d = time2 - time1;

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  for (int i = 0; i < n; i++) {
    fft1->timing1d(work1,nfft_both,1);
    fft2->timing1d(work1,nfft_both,-1);
    if (differentiation_flag != 1) {
      fft2->timing1d(work1,nfft_both,-1);
      fft2->timing1d(work1,nfft_both,-1);
    }
  }

  MPI_Barrier(world);
  time2 = MPI_Wtime();
  time1d = time2 - time1;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double PPPM::memory_usage()
{
  double bytes = nmax*3 * sizeof(double);
  int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);
  if (differentiation_flag == 1) {
    bytes += 2 * nbrick * sizeof(FFT_SCALAR);
  } else {
    bytes += 4 * nbrick * sizeof(FFT_SCALAR);
  }
  bytes += 6 * nfft_both * sizeof(double);
  bytes += nfft_both * sizeof(double);
  bytes += nfft_both*5 * sizeof(FFT_SCALAR);
  bytes += 2 * nbuf * sizeof(FFT_SCALAR);

  if (peratom_allocate_flag) {
    bytes += 6 * nbrick * sizeof(FFT_SCALAR);
    bytes += 2 * nbuf_peratom * sizeof(FFT_SCALAR);
  }

  if (group_allocate_flag) {
    bytes += 2 * nbrick * sizeof(FFT_SCALAR);
    bytes += 2 * nfft_both * sizeof(FFT_SCALAR);;
  }

  return bytes;
}

/* ----------------------------------------------------------------------
   group-group interactions
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   compute the PPPM total long-range force and energy for groups A and B
 ------------------------------------------------------------------------- */

void PPPM::compute_group_group(int groupbit_A, int groupbit_B, int BA_flag)
{
  if (slabflag)
    error->all(FLERR,"Cannot (yet) use K-space slab "
               "correction with compute group/group");

  int i,j;

  if (!group_allocate_flag) {
    allocate_groups();
    group_allocate_flag = 1;
  }

  e2group = 0; //energy
  f2group[0] = 0; //force in x-direction
  f2group[1] = 0; //force in y-direction
  f2group[2] = 0; //force in z-direction

  double *q = atom->q;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;


  // map my particle charge onto my local 3d density grid

  make_rho_groups(groupbit_A,groupbit_B,BA_flag);

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  // temporarily store and switch pointers so we can
  //  use brick2fft() for groups A and B (without
  //  writing an additional function)

  FFT_SCALAR ***density_brick_real = density_brick;
  FFT_SCALAR *density_fft_real = density_fft;

  // group A

  density_brick = density_A_brick;
  density_fft = density_A_fft;

  brick2fft();

  // group B

  density_brick = density_B_brick;
  density_fft = density_B_fft;

  brick2fft();

  // switch back pointers

  density_brick = density_brick_real;
  density_fft = density_fft_real;

  // compute potential gradient on my FFT grid and
  //   portion of group-group energy/force on this proc's FFT grid

  poisson_groups(BA_flag);

  const double qscale = force->qqrd2e * scale;

  // total group A <--> group B energy
  // self and boundary correction terms are in compute_group_group.cpp

  double e2group_all;
  MPI_Allreduce(&e2group,&e2group_all,1,MPI_DOUBLE,MPI_SUM,world);
  e2group = e2group_all;

  e2group *= qscale*0.5*volume;

  // total group A <--> group B force

  double f2group_all[3];
  MPI_Allreduce(f2group,f2group_all,3,MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < 3; i++) f2group[i] = qscale*volume*f2group_all[i];
}

/* ----------------------------------------------------------------------
 allocate group-group memory that depends on # of K-vectors and order
 ------------------------------------------------------------------------- */

void PPPM::allocate_groups()
{
  memory->create3d_offset(density_A_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_A_brick");
  memory->create3d_offset(density_B_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_B_brick");
  memory->create(density_A_fft,nfft_both,"pppm:density_A_fft");
  memory->create(density_B_fft,nfft_both,"pppm:density_B_fft");
}

/* ----------------------------------------------------------------------
 deallocate group-group memory that depends on # of K-vectors and order
 ------------------------------------------------------------------------- */

void PPPM::deallocate_groups()
{
  memory->destroy3d_offset(density_A_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(density_B_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy(density_A_fft);
  memory->destroy(density_B_fft);
}

/* ----------------------------------------------------------------------
 create discretized "density" on section of global grid due to my particles
 density(x,y,z) = charge "density" at grid points of my 3d brick
 (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
 in global grid for group-group interactions
 ------------------------------------------------------------------------- */

void PPPM::make_rho_groups(int groupbit_A, int groupbit_B, int BA_flag)
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density arrays

  memset(&(density_A_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  memset(&(density_B_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  for (int i = 0; i < nlocal; i++) {

    if ((mask[i] & groupbit_A) && (mask[i] & groupbit_B))
      if (BA_flag) continue;

    if ((mask[i] & groupbit_A) || (mask[i] & groupbit_B)) {

      nx = part2grid[i][0];
      ny = part2grid[i][1];
      nz = part2grid[i][2];
      dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
      dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
      dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

      compute_rho1d(dx,dy,dz);

      z0 = delvolinv * q[i];
      for (n = nlower; n <= nupper; n++) {
        mz = n+nz;
        y0 = z0*rho1d[2][n];
        for (m = nlower; m <= nupper; m++) {
          my = m+ny;
          x0 = y0*rho1d[1][m];
          for (l = nlower; l <= nupper; l++) {
            mx = l+nx;

            // group A

            if (mask[i] & groupbit_A)
              density_A_brick[mz][my][mx] += x0*rho1d[0][l];

            // group B

            if (mask[i] & groupbit_B)
              density_B_brick[mz][my][mx] += x0*rho1d[0][l];
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for group-group interactions
 ------------------------------------------------------------------------- */

void PPPM::poisson_groups(int BA_flag)
{
  int i,j,k,n;
  double eng;

  // reuse memory (already declared)

  FFT_SCALAR *work_A = work1;
  FFT_SCALAR *work_B = work2;

  // transform charge density (r -> k)

  // group A

  n = 0;
  for (i = 0; i < nfft; i++) {
    work_A[n++] = density_A_fft[i];
    work_A[n++] = ZEROF;
  }

  fft1->compute(work_A,work_A,1);

  // group B

  n = 0;
  for (i = 0; i < nfft; i++) {
    work_B[n++] = density_B_fft[i];
    work_B[n++] = ZEROF;
  }

  fft1->compute(work_B,work_B,1);

  // group-group energy and force contribution,
  //  keep everything in reciprocal space so
  //  no inverse FFTs needed

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
  double s2 = scaleinv*scaleinv;

  // energy

  n = 0;
  for (i = 0; i < nfft; i++) {
    e2group += s2 * greensfn[i] *
      (work_A[n]*work_B[n] + work_A[n+1]*work_B[n+1]);
    n += 2;
  }

  if (BA_flag) return;


  // multiply by Green's function and s2
  //  (only for work_A so it is not squared below)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work_A[n++] *= s2 * greensfn[i];
    work_A[n++] *= s2 * greensfn[i];
  }

  double partial_group;

  // force, x direction

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        partial_group = work_A[n+1]*work_B[n] - work_A[n]*work_B[n+1];
        f2group[0] += fkx[i] * partial_group;
        n += 2;
      }

  // force, y direction

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        partial_group = work_A[n+1]*work_B[n] - work_A[n]*work_B[n+1];
        f2group[1] += fky[j] * partial_group;
        n += 2;
      }

  // force, z direction

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        partial_group = work_A[n+1]*work_B[n] - work_A[n]*work_B[n+1];
        f2group[2] += fkz[k] * partial_group;
        n += 2;
      }
}
