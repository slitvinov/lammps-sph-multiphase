#include "fix_phase_change.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_park.h"
#include "region.h"
#include "sph_multiphase.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace SPHMultiphase;

static constexpr double CG_SMALL = 1.0e-20;

FixPhaseChange::FixPhaseChange(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), idregion(nullptr), dmass(nullptr), list(nullptr), random(nullptr)
{
  if (narg < 14) error->all(FLERR, "Illegal fix phase_change command");
  comm_reverse = 1;
  comm_forward = 1;
  restart_global = 1;
  time_depend = 1;

  int m = 3;
  Tc = utils::numeric(FLERR, arg[m++], false, lmp);
  Tt = utils::numeric(FLERR, arg[m++], false, lmp);
  Hwv = utils::numeric(FLERR, arg[m++], false, lmp);
  dr = utils::numeric(FLERR, arg[m++], false, lmp);
  to_mass = utils::numeric(FLERR, arg[m++], false, lmp);
  cutoff = utils::numeric(FLERR, arg[m++], false, lmp);
  from_type = utils::inumeric(FLERR, arg[m++], false, lmp);
  to_type = utils::inumeric(FLERR, arg[m++], false, lmp);
  nfreq = utils::inumeric(FLERR, arg[m++], false, lmp);
  seed = utils::inumeric(FLERR, arg[m++], false, lmp);
  if (seed <= 0) error->all(FLERR, "Illegal value for seed");
  if (strcmp(arg[m], "ENERGY") == 0) {
    energy_chance_flag = true;
    if (m + 1 >= narg) error->all(FLERR, "Illegal fix phase_change command");
    phase_change_rate = utils::numeric(FLERR, arg[m + 1], false, lmp);
    m += 2;
  } else {
    energy_chance_flag = false;
    change_chance = utils::numeric(FLERR, arg[m++], false, lmp);
    if (change_chance < 0) error->all(FLERR, "Illegal value for change_chance");
  }

  maxattempt = 10;
  convert_flag = 0;
  convert_cv = 0.0;
  while (m < narg) {
    if (strcmp(arg[m], "region") == 0) {
      if (m + 2 > narg) error->all(FLERR, "Illegal fix phase_change command");
      delete[] idregion;
      idregion = utils::strdup(arg[m + 1]);
      m += 2;
    } else if (strcmp(arg[m], "attempt") == 0) {
      if (m + 2 > narg) error->all(FLERR, "Illegal fix phase_change command");
      maxattempt = utils::inumeric(FLERR, arg[m + 1], false, lmp);
      m += 2;
    } else if (strcmp(arg[m], "convert") == 0) {
      if (m + 2 > narg) error->all(FLERR, "Illegal fix phase_change command");
      convert_flag = 1;
      convert_cv = utils::numeric(FLERR, arg[m + 1], false, lmp);
      m += 2;
    } else if (strcmp(arg[m], "units") == 0) {
      if (m + 2 > narg || strcmp(arg[m + 1], "box") != 0)
        error->all(FLERR, "Fix phase_change supports only units box");
      m += 2;
    } else
      error->all(FLERR, "Illegal fix phase_change command");
  }

  if (!idregion) error->all(FLERR, "Must specify a region in fix phase_change");
  auto region = domain->get_region_by_id(idregion);
  if (!region) error->all(FLERR, "Region ID for fix phase_change does not exist");
  if (region->bboxflag == 0)
    error->all(FLERR, "Fix phase_change region does not support a bounding box");
  if (region->dynamic_check()) error->all(FLERR, "Fix phase_change region cannot be dynamic");
  if (region->extent_xlo < domain->boxlo[0] || region->extent_xhi > domain->boxhi[0] ||
      region->extent_ylo < domain->boxlo[1] || region->extent_yhi > domain->boxhi[1] ||
      region->extent_zlo < domain->boxlo[2] || region->extent_zhi > domain->boxhi[2])
    error->all(FLERR, "Phase change region extends outside simulation box");

  random = new RanPark(lmp, seed + comm->me);
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
}

FixPhaseChange::~FixPhaseChange()
{
  delete random;
  delete[] idregion;
}

int FixPhaseChange::setmask()
{
  return PRE_EXCHANGE;
}

void FixPhaseChange::init()
{
  require_rmass(atom, error);
  colorgradient(atom, error);
  if (!domain->get_region_by_id(idregion))
    error->all(FLERR, "Region ID for fix phase_change does not exist");
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

void FixPhaseChange::init_list(int, NeighList *ptr)
{
  list = ptr;
}

void FixPhaseChange::pre_exchange()
{
  if (next_reneighbor != update->ntimestep) return;
  comm->forward_comm(this);

  double *sublo = domain->triclinic ? domain->sublo_lamda : domain->sublo;
  double *subhi = domain->triclinic ? domain->subhi_lamda : domain->subhi;
  const int dim = domain->dimension;

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double **x = atom->x;
  double **v = atom->v;
  double **vest = atom->vest;
  double **cg = colorgradient(atom, error);
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  double *cv = atom->cv;
  double *e = atom->esph;
  int *type = atom->type;
  dmass = atom->drho;
  for (int i = 0; i < nall; i++) dmass[i] = 0.0;

  std::vector<double> pending;
  for (int i = 0; i < nlocal; i++) {
    double Ti = energy2t(e[i], cv[i]);
    bool isphasechange;
    if (Ti < Tc || type[i] != to_type)
      isphasechange = false;
    else if (energy_chance_flag) {
      double threshold = (e[i] - t2energy(Tc, cv[i])) / Hwv * update->dt * phase_change_rate;
      isphasechange = random->uniform() < threshold && isfromphasearound(i);
    } else
      isphasechange = random->uniform() < change_chance && Ti > Tt && isfromphasearound(i);
    if (!isphasechange) continue;

    double coord[3];
    bool ok = false;
    double delta = dr;
    for (int n = 0; n < maxattempt && !ok; n++, delta *= 0.75) {
      create_newpos(x[i], cg[i], delta, coord);
      ok = owns(coord, sublo, subhi);
    }
    delta = dr;
    for (int n = 0; n < maxattempt && !ok; n++, delta *= 0.75) {
      create_newpos_simple(x[i], delta, coord);
      ok = owns(coord, sublo, subhi);
    }
    if (!ok) continue;

    int *jlist = list->firstneigh[i];
    int jnum = list->numneigh[i];
    std::vector<double> wj(jnum, 0.0);
    double wtotal = 0.0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      if (type[j] != from_type || rmass[j] - dmass[j] <= to_mass) continue;
      double delx = x[i][0] - x[j][0];
      double dely = x[i][1] - x[j][1];
      double delz = x[i][2] - x[j][2];
      wj[jj] = w(dim, sqrt(delx * delx + dely * dely + delz * delz) / cutoff);
      wtotal += wj[jj];
    }
    if (wtotal <= 0.0) continue;

    double dmom[3] = {0.0, 0.0, 0.0};
    double dmomest[3] = {0.0, 0.0, 0.0};
    double denergy = 0.0;
    for (int jj = 0; jj < jnum; jj++) {
      if (wj[jj] == 0.0) continue;
      int j = jlist[jj] & NEIGHMASK;
      double dm = to_mass * wj[jj] / wtotal;
      dmass[j] += dm;
      denergy += e[j] * dm;
      for (int d = 0; d < 3; d++) {
        dmom[d] += v[j][d] * dm;
        dmomest[d] += vest[j][d] * dm;
      }
    }

    double enew = (rmass[i] * e[i] + denergy - to_mass * Hwv) / (rmass[i] + to_mass);
    e[i] = enew;
    double one[12] = {coord[0],          coord[1],          coord[2],          dmom[0] / to_mass,
                      dmom[1] / to_mass, dmom[2] / to_mass, dmomest[0] / to_mass, dmomest[1] / to_mass,
                      dmomest[2] / to_mass, rho[i],         cv[i],             enew};
    pending.insert(pending.end(), one, one + 12);
  }

  comm->reverse_comm(this);
  for (int i = 0; i < nlocal; i++) {
    rmass[i] -= dmass[i];
    if (rmass[i] <= 0) error->one(FLERR, "Fix phase_change made a negative mass");
    dmass[i] = 0.0;
    if (convert_flag && type[i] == from_type && rmass[i] < 2.0 * to_mass) {
      type[i] = to_type;
      atom->mask[i] |= groupbit;
      e[i] -= Hwv;
      cv[i] = convert_cv;
    }
  }

  int nins = 0;
  for (size_t k = 0; k < pending.size(); k += 12) {
    double *one = &pending[k];
    if (!insert_one_atom(one, sublo, subhi)) continue;
    nins++;
    int m = atom->nlocal - 1;
    atom->rmass[m] = to_mass;
    for (int d = 0; d < 3; d++) {
      atom->v[m][d] = one[3 + d];
      atom->vest[m][d] = one[6 + d];
    }
    atom->rho[m] = one[9];
    atom->cv[m] = one[10];
    atom->esph[m] = one[11];
  }

  next_reneighbor += nfreq;
  int ninsall;
  MPI_Allreduce(&nins, &ninsall, 1, MPI_INT, MPI_SUM, world);
  if (ninsall > 0) {
    atom->natoms += ninsall;
    if (atom->tag_enable) atom->tag_extend();
    atom->nghost = 0;
    if (atom->map_style != Atom::MAP_NONE) {
      atom->map_init();
      atom->map_set();
    }
  }
}

bool FixPhaseChange::owns(double *coord, double *sublo, double *subhi)
{
  double lamda[3];
  double *c = coord;
  if (domain->triclinic) {
    domain->x2lamda(coord, lamda);
    c = lamda;
  }
  if (c[0] >= sublo[0] && c[0] < subhi[0] && c[1] >= sublo[1] && c[1] < subhi[1] &&
      c[2] >= sublo[2] && c[2] < subhi[2])
    return true;
  if (domain->dimension == 3 && c[2] >= domain->boxhi[2] &&
      comm->myloc[2] == comm->procgrid[2] - 1 && c[0] >= sublo[0] && c[0] < subhi[0] &&
      c[1] >= sublo[1] && c[1] < subhi[1])
    return true;
  if (domain->dimension == 2 && c[1] >= domain->boxhi[1] &&
      comm->myloc[1] == comm->procgrid[1] - 1 && c[0] >= sublo[0] && c[0] < subhi[0])
    return true;
  return false;
}

bool FixPhaseChange::insert_one_atom(double *coord, double *sublo, double *subhi)
{
  if (!owns(coord, sublo, subhi)) return false;
  atom->avec->create_atom(to_type, coord);
  int m = atom->nlocal - 1;
  atom->mask[m] = 1 | groupbit;
  for (auto &fix : modify->get_fix_list())
    if (fix->create_attribute) fix->set_arrays(m);
  return true;
}

void FixPhaseChange::create_newpos_simple(double *xone, double delta, double *coord)
{
  coord[0] = xone[0] + (random->uniform() - 0.5) * delta;
  coord[1] = xone[1] + (random->uniform() - 0.5) * delta;
  coord[2] = xone[2] + (random->uniform() - 0.5) * delta;
}

void FixPhaseChange::create_newpos(double *xone, double *g, double delta, double *coord)
{
  double eij[3];
  if (domain->dimension == 3) {
    double a[3] = {1.0, 0.0, 0.0};
    if (fabs(g[0]) > fabs(g[1]) || fabs(g[0]) > fabs(g[2])) {
      a[0] = 0.0;
      if (fabs(g[1]) < fabs(g[2]))
        a[1] = 1.0;
      else
        a[2] = 1.0;
    }
    double b1[3] = {g[1] * a[2] - g[2] * a[1], g[2] * a[0] - g[0] * a[2],
                    g[0] * a[1] - g[1] * a[0]};
    double b2[3] = {g[1] * b1[2] - g[2] * b1[1], g[2] * b1[0] - g[0] * b1[2],
                    g[0] * b1[1] - g[1] * b1[0]};
    double b1abs = sqrt(b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2]);
    double b2abs = sqrt(b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]);
    if (b1abs < CG_SMALL || b2abs < CG_SMALL) {
      create_newpos_simple(xone, delta, coord);
      return;
    }
    double atmp = random->uniform() - 0.5;
    double btmp = random->uniform() - 0.5;
    for (int d = 0; d < 3; d++) eij[d] = atmp * b1[d] / b1abs + btmp * b2[d] / b2abs;
  } else {
    double s = random->uniform() > 0.5 ? 1.0 : -1.0;
    eij[0] = -s * g[1];
    eij[1] = s * g[0];
    eij[2] = 0.0;
  }
  double eijabs = sqrt(eij[0] * eij[0] + eij[1] * eij[1] + eij[2] * eij[2]);
  if (eijabs < CG_SMALL) {
    create_newpos_simple(xone, delta, coord);
    return;
  }
  for (int d = 0; d < 3; d++) coord[d] = xone[d] + eij[d] * delta / eijabs;
}

bool FixPhaseChange::isfromphasearound(int i)
{
  double **x = atom->x;
  int *type = atom->type;
  int *jlist = list->firstneigh[i];
  int jnum = list->numneigh[i];
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj] & NEIGHMASK;
    if (type[j] != from_type) continue;
    double delx = x[i][0] - x[j][0];
    double dely = x[i][1] - x[j][1];
    double delz = x[i][2] - x[j][2];
    if (delx * delx + dely * dely + delz * delz <= cutoff * cutoff) return true;
  }
  return false;
}

void FixPhaseChange::write_restart(FILE *fp)
{
  double list[3] = {(double) random->state(), (double) nfirst, (double) next_reneighbor};
  if (comm->me == 0) {
    int size = 3 * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), 3, fp);
  }
}

void FixPhaseChange::restart(char *buf)
{
  auto list = (double *) buf;
  seed = static_cast<int>(list[0]);
  nfirst = static_cast<int>(list[1]);
  next_reneighbor = static_cast<bigint>(list[2]);
  random->reset(seed + comm->me);
}

int FixPhaseChange::pack_forward_comm(int n, int *list, double *buf, int, int *)
{
  double *e = atom->esph;
  for (int i = 0; i < n; i++) buf[i] = e[list[i]];
  return n;
}

void FixPhaseChange::unpack_forward_comm(int n, int first, double *buf)
{
  double *e = atom->esph;
  for (int i = 0; i < n; i++) e[first + i] = buf[i];
}

int FixPhaseChange::pack_reverse_comm(int n, int first, double *buf)
{
  for (int i = 0; i < n; i++) buf[i] = dmass[first + i];
  return n;
}

void FixPhaseChange::unpack_reverse_comm(int n, int *list, double *buf)
{
  for (int i = 0; i < n; i++) dmass[list[i]] += buf[i];
}
