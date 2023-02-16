module m_qednamespace
  use m_globalnamespace
  implicit none

#ifdef RADIATION
  real :: emit_gamma_syn, emit_gamma_ic, cool_gamma_syn, cool_gamma_ic, rad_beta_rec
  real :: rad_dens_lim, rad_cool_lim
  integer :: rad_photon_sp
  integer :: rad_interval
#endif

#ifdef QED
  real :: QED_tau0
#endif

#ifdef BWPAIRPRODUCTION
  ! BW pair production parameters
  integer :: BW_interval, BW_electron_sp, BW_positron_sp
  integer :: BW_algorithm
#endif

#ifdef COMPTONSCATTERING
  ! Compton scattering parameters
  integer :: Compton_interval, Compton_algorithm
  logical :: Compton_el_recoil
  real(kind=8) :: Thomson_lim
  real :: Compton_nph_over_ne
  logical :: Compton_clone_sets
  logical :: Compton_cool_el
#endif

#ifdef PAIRANNIHILATION
  ! Pair annihilation parameters
  integer :: Annihilation_interval, Annihilation_photon_sp
  integer :: Annihilation_algorithm
  logical :: Annihilation_sporadic
#endif

end module m_qednamespace
