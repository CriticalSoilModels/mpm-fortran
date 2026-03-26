module mpm_update_config
   !! Particle velocity update scheme configuration.
   !!
   !! Governs how grid velocities are mapped back to particles in the G2P step
   !! (and for APIC, also affects P2G). See notes/update_schemes.md for a full
   !! comparison of schemes, dissipation properties, and implementation notes.
   use mpm_precision, only: wp
   implicit none
   private

   ! --- scheme constants ---
   integer, parameter, public :: SCHEME_FLIP_PIC = 1  !! alpha-FLIP/PIC blend (default)
   integer, parameter, public :: SCHEME_APIC     = 2  !! Affine PIC (Jiang et al. 2015)
   integer, parameter, public :: SCHEME_XPIC     = 3  !! eXtended PIC (Hammerquist & Nairn 2017)

   public :: mpm_update_config_t
   public :: is_valid_config

   type :: mpm_update_config_t
      !! Configuration for the particle velocity update scheme.
      !!
      !! Default is alpha-FLIP with alpha=0.99 — low dissipation with slow
      !! bleed of particle noise. Set alpha=1.0 for elastic wave problems
      !! (pure FLIP, no artificial damping). Set alpha=0.0 for pure PIC
      !! (maximum dissipation, maximum stability).
      integer  :: scheme = SCHEME_FLIP_PIC  !! active update scheme
      real(wp) :: alpha  = 0.99_wp          !! FLIP/PIC blend [0,1]; ignored for APIC/XPIC
   end type mpm_update_config_t

contains

   pure function is_valid_config(config) result(valid)
      !! Return .true. if the configuration is self-consistent.
      !! alpha must be in [0, 1] for FLIP_PIC; scheme must be a known constant.
      type(mpm_update_config_t), intent(in) :: config
      logical :: valid

      valid = .false.

      select case (config%scheme)
         case (SCHEME_FLIP_PIC)
            valid = (config%alpha >= 0.0_wp) .and. (config%alpha <= 1.0_wp)
         case (SCHEME_APIC, SCHEME_XPIC)
            valid = .true.
         case default
            valid = .false.
      end select

   end function is_valid_config

end module mpm_update_config
