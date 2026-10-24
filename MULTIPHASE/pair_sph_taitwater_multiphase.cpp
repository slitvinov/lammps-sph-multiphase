#include "pair_sph_taitwater_multiphase.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "sph_multiphase.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace SPHMultiphase;

PairSPHTaitwaterMultiphase::PairSPHTaitwaterMultiphase(LAMMPS *lmp) :
    Pair(lmp), rho0(nullptr), soundspeed(nullptr), B(nullptr), rbackground(nullptr),
    gamma(nullptr), cut(nullptr), viscosity(nullptr)
{
  restartinfo = 0;
}

PairSPHTaitwaterMultiphase::~PairSPHTaitwaterMultiphase()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(gamma);
    memory->destroy(rbackground);
    memory->destroy(viscosity);
  }
}

void PairSPHTaitwaterMultiphase::init_style()
{
  require_rmass(atom, error);
  Pair::init_style();
}

void PairSPHTaitwaterMultiphase::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const int dim = domain->dimension;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int itype = type[i];
    double pi = pressure(B[itype], rho0[itype], gamma[itype], rbackground[itype], rho[i]);
    double Vi = rmass[i] / rho[i];
    double Vi2 = Vi * Vi;
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      int jtype = type[j];
      double delx = x[i][0] - x[j][0];
      double dely = x[i][1] - x[j][1];
      double delz = x[i][2] - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;
      if (rsq >= cutsq[itype][jtype]) continue;
      double r = sqrt(rsq);
      double ih = 1.0 / cut[itype][jtype];
      double wfd = dw(dim, r * ih) * (dim == 3 ? ih * ih * ih * ih : ih * ih * ih) / r;
      double Vj = rmass[j] / rho[j];
      double Vj2 = Vj * Vj;
      double pj = pressure(B[jtype], rho0[jtype], gamma[jtype], rbackground[jtype], rho[j]);
      double pij = (rho[j] * pi + rho[i] * pj) / (rho[i] + rho[j]);
      double velx = v[i][0] - v[j][0];
      double vely = v[i][1] - v[j][1];
      double velz = v[i][2] - v[j][2];
      double fvisc = (Vi2 + Vj2) * viscosity[itype][jtype] * wfd;
      double fpair = -(Vi2 + Vj2) * pij * wfd;
      f[i][0] += delx * fpair + velx * fvisc;
      f[i][1] += dely * fpair + vely * fvisc;
      f[i][2] += delz * fpair + velz * fvisc;
      if (newton_pair || j < nlocal) {
        f[j][0] -= delx * fpair + velx * fvisc;
        f[j][1] -= dely * fpair + vely * fvisc;
        f[j][2] -= delz * fpair + velz * fvisc;
      }
      if (evflag) ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

void PairSPHTaitwaterMultiphase::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(rho0, n + 1, "pair:rho0");
  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(gamma, n + 1, "pair:gamma");
  memory->create(rbackground, n + 1, "pair:rbackground");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");
}

void PairSPHTaitwaterMultiphase::settings(int narg, char **)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style sph/taitwater/multiphase command");
}

void PairSPHTaitwaterMultiphase::coeff(int narg, char **arg)
{
  if (narg != 8)
    error->all(FLERR, "Incorrect args for pair_style sph/taitwater/multiphase coefficients");
  if (!allocated) allocate();
  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  double rho0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double soundspeed_one = utils::numeric(FLERR, arg[3], false, lmp);
  double viscosity_one = utils::numeric(FLERR, arg[4], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[5], false, lmp);
  double cut_one = utils::numeric(FLERR, arg[6], false, lmp);
  double rbackground_one = utils::numeric(FLERR, arg[7], false, lmp);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / gamma_one;
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    gamma[i] = gamma_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    rbackground[i] = rbackground_one;
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

double PairSPHTaitwaterMultiphase::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, "All pair sph/taitwater/multiphase coeffs are not set");
  cut[j][i] = cut[i][j];
  viscosity[j][i] = viscosity[i][j];
  return cut[i][j];
}
