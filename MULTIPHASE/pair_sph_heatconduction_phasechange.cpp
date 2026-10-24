#include "pair_sph_heatconduction_phasechange.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "sph_multiphase.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace SPHMultiphase;

PairSPHHeatConductionPhaseChange::PairSPHHeatConductionPhaseChange(LAMMPS *lmp) : Pair(lmp), cut(nullptr), alpha(nullptr)
{
  restartinfo = 0;
}

PairSPHHeatConductionPhaseChange::~PairSPHHeatConductionPhaseChange()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(alpha);
    memory->destroy(tc);
    memory->destroy(fixflag);
  }
}

void PairSPHHeatConductionPhaseChange::init_style()
{
  require_rmass(atom, error);
  if (!atom->cv_flag) error->all(FLERR, "Pair sph/heatconduction/phasechange requires atom attribute cv");
  Pair::init_style();
}

void PairSPHHeatConductionPhaseChange::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  double **x = atom->x;
  double *e = atom->esph;
  double *de = atom->desph;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  double *cv = atom->cv;
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
      double Ti = energy2t(e[i], cv[i]);
      double Tj = energy2t(e[j], cv[j]);
      if (fixflag[itype][jtype] == itype && Ti < Tj) Ti = tc[itype][jtype];
      if (fixflag[itype][jtype] == jtype && Tj < Ti) Tj = tc[itype][jtype];
      double deltaE = 2.0 * alpha[itype][jtype] * (Ti - Tj) * wfd / (rho[i] * rho[j]);
      de[i] += deltaE * rmass[j];
      if (newton_pair || j < nlocal) de[j] -= deltaE * rmass[i];
    }
  }
}

void PairSPHHeatConductionPhaseChange::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");
  memory->create(tc, n + 1, n + 1, "pair:tc");
  memory->create(fixflag, n + 1, n + 1, "pair:fixflag");
  for (int i = 0; i <= n; i++)
    for (int j = 0; j <= n; j++) {
      tc[i][j] = 0.0;
      fixflag[i][j] = 0;
    }
}

void PairSPHHeatConductionPhaseChange::settings(int narg, char **)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style sph/heatconduction/phasechange command");
}

void PairSPHHeatConductionPhaseChange::coeff(int narg, char **arg)
{
  if (narg != 4 && narg != 6)
    error->all(FLERR, "Incorrect args for pair_style sph/heatconduction/phasechange coefficients");
  if (!allocated) allocate();
  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  double alpha_one = utils::numeric(FLERR, arg[2], false, lmp);
  double cut_one = utils::numeric(FLERR, arg[3], false, lmp);
  int fixflag_one = 0;
  double tc_one = 0.0;
  if (narg == 6) {
    if (strcmp(arg[4], "NULL") == 0) {
      fixflag_one = 2;
      tc_one = utils::numeric(FLERR, arg[5], false, lmp);
    } else if (strcmp(arg[5], "NULL") == 0) {
      fixflag_one = 1;
      tc_one = utils::numeric(FLERR, arg[4], false, lmp);
    } else
      error->all(FLERR, "Incorrect args for pair coefficients");
  }
  int count = 0;
  for (int i = ilo; i <= ihi; i++)
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      alpha[i][j] = alpha_one;
      fixflag[i][j] = fixflag_one == 1 ? i : (fixflag_one == 2 ? j : 0);
      tc[i][j] = tc_one;
      setflag[i][j] = 1;
      count++;
    }
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

double PairSPHHeatConductionPhaseChange::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, "All pair sph/heatconduction/phasechange coeffs are not set");
  cut[j][i] = cut[i][j];
  alpha[j][i] = alpha[i][j];
  tc[j][i] = tc[i][j];
  fixflag[j][i] = fixflag[i][j];
  return cut[i][j];
}
