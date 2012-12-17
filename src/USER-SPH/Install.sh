# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p fix_phase_change.cpp ..
  cp -p pair_sph_colorgradient.cpp ..
  cp -p atom_vec_meso.cpp ..
  cp -p pair_sph_heatconduction.cpp ..
  cp -p pair_sph_idealgas.cpp ..
  cp -p pair_sph_lj.cpp ..
  cp -p pair_sph_rhosum.cpp ..
  cp -p pair_sph_rhosum_multiphase.cpp ..
  cp -p pair_sph_taitwater.cpp ..
  cp -p pair_sph_taitwater_morris.cpp ..
  cp -p pair_sph_surfacetension.cpp ..
  cp -p sph_kernel_quintic.cpp ..
  cp -p sph_energy_equation.cpp ..
  cp -p compute_meso_e_atom.cpp ..
  cp -p compute_meso_rho_atom.cpp ..
  cp -p compute_meso_colorgradient_atom.cpp ..
  cp -p compute_meso_t_atom.cpp ..
  cp -p fix_meso.cpp ..
  cp -p fix_setmeso.cpp ..
  cp -p fix_setmesode.cpp ..
  cp -p fix_meso_stationary.cpp ..

  cp -p fix_phase_change.h ..
  cp -p pair_sph_colorgradient.h ..
  cp -p atom_vec_meso.h ..
  cp -p pair_sph_heatconduction.h ..
  cp -p pair_sph_idealgas.h ..
  cp -p pair_sph_lj.h ..
  cp -p pair_sph_rhosum.h ..
  cp -p pair_sph_rhosum_multiphase.h ..
  cp -p pair_sph_taitwater.h ..
  cp -p pair_sph_taitwater_morris.h ..
  cp -p pair_sph_surfacetension.h ..
  cp -p sph_kernel_quintic.h ..
  cp -p sph_energy_equation.h ..
  cp -p compute_meso_e_atom.h ..
  cp -p compute_meso_rho_atom.h ..
  cp -p compute_meso_colorgradient_atom.h ..
  cp -p compute_meso_t_atom.h ..
  cp -p fix_setmeso.h ..
  cp -p fix_setmesode.h ..
  cp -p fix_meso.h ..
  cp -p fix_meso_stationary.h ..

elif (test $1 = 0) then
  rm -f ../fix_phase_change.cpp
  rm -f ../pair_sph_colorgradient.cpp
  rm -f ../atom_vec_meso.cpp
  rm -f ../pair_sph_heatconduction.cpp
  rm -f ../pair_sph_idealgas.cpp
  rm -f ../pair_sph_lj.cpp
  rm -f ../pair_sph_rhosum.cpp
  rm -f ../pair_sph_rhosum_multiphase.cpp
  rm -f ../pair_sph_taitwater.cpp
  rm -f ../pair_sph_taitwater_morris.cpp
  rm -f ../pair_sph_surfacetension.cpp
  rm -f ../compute_meso_e_atom.cpp
  rm -f ../compute_meso_rho_atom.cpp
  rm -f ../compute_meso_colorgradient_atom.cpp
  rm -f ../compute_meso_t_atom.cpp
  rm -f ../sph_kernel_quintic.cpp
  rm -f ../sph_energy_equation.cpp
  rm -f ../fix_setmeso.cpp
  rm -f ../fix_setmesode.cpp
  rm -f ../fix_meso.cpp
  rm -f ../fix_meso_stationary.cpp

  rm -f ../fix_phase_change.h
  rm -f ../pair_sph_colorgradient.h
  rm -f ../atom_vec_meso.h
  rm -f ../pair_sph_heatconduction.h
  rm -f ../pair_sph_idealgas.h
  rm -f ../pair_sph_lj.h
  rm -f ../pair_sph_rhosum.h
  rm -f ../pair_sph_rhosum_multiphase.h
  rm -f ../pair_sph_taitwater.h
  rm -f ../pair_sph_taitwater_morris.h
  rm -f ../pair_sph_surfacetension.h
  rm -f ../compute_meso_e_atom.h
  rm -f ../compute_meso_rho_atom.h
  rm -f ../compute_meso_colorgradient_atom.h
  rm -f ../compute_meso_t_atom.h
  rm -f ../sph_kernel_quintic.h
  rm -f ../sph_energy_equation.h
  rm -f ../fix_meso.h
  rm -f ../fix_setmeso.h
  rm -f ../fix_setmesode.h
  rm -f ../fix_meso_stationary.h

fi
