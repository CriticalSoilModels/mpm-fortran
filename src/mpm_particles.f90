module mpm_particles
   !! Material point storage and phase configuration.
   !!
   !! Particles are stored in SOA (Structure of Arrays) layout: first index is
   !! the component, second is the particle number. This ensures contiguous
   !! memory access per component in GPU kernels (Fortran column-major).
   !!
   !! Fluid arrays (uw, ua, Sr) are only allocated when the phase configuration
   !! requests them. Single-phase problems pay no allocation cost for multiphase.
   use mpm_precision, only: wp, n_voigt, n_dims
   implicit none
   private

   ! --- phase labels ---
   integer, parameter, public :: PHASE_SOLID  = 1
   integer, parameter, public :: PHASE_LIQUID = 2
   integer, parameter, public :: PHASE_GAS    = 3

   public :: mpm_phase_config_t
   public :: mpm_particle_set_t
   public :: allocate_particle_set
   public :: deallocate_particle_set

   type :: mpm_phase_config_t
      !! Controls which phases are active in the simulation.
      !! Fluid arrays are only allocated when the corresponding flag is true.
      logical :: has_liquid  = .false.  !! allocate pore water fields (uw, Sr)
      logical :: has_gas     = .false.  !! allocate pore air field (ua)
      logical :: has_contact = .false.  !! enable contact detection [future]
   end type mpm_phase_config_t

   type :: mpm_particle_set_t
      !! Material point storage for one phase in SOA layout.
      !!
      !! All arrays have shape (component, n_mp). Stress uses Voigt ordering
      !! [11,22,33,12,13,23] consistent with critical-soil-models (3D),
      !! or the first nvoigt components for lower-dimensional problems.

      ! --- configuration (set at allocation, read-only thereafter) ---
      integer :: phase   = PHASE_SOLID  !! PHASE_SOLID, PHASE_LIQUID, PHASE_GAS
      integer :: body_id = 1            !! body identifier for contact detection [future]
      integer :: n_mp    = 0            !! number of material points
      integer :: ndim    = 0            !! spatial dimension
      integer :: nvoigt  = 0            !! Voigt stress/strain components
      integer :: nstatev = 0            !! constitutive state variable count

      ! --- kinematics (always allocated) ---
      real(wp), allocatable :: xp(:,:)   !! position    [m]      (ndim, n_mp)
      real(wp), allocatable :: vp(:,:)   !! velocity    [m/s]    (ndim, n_mp)
      real(wp), allocatable :: ap(:,:)   !! acceleration [m/s^2] (ndim, n_mp)
      real(wp), allocatable :: mp(:)     !! mass        [kg]     (n_mp)
      real(wp), allocatable :: vol(:)    !! volume      [m^d]    (n_mp)
      real(wp), allocatable :: rho(:)    !! density     [kg/m^d] (n_mp)

      ! --- solid mechanics (always allocated) ---
      real(wp), allocatable :: sig(:,:)    !! Cauchy stress  [kPa] (nvoigt, n_mp)
      real(wp), allocatable :: eps(:,:)    !! total strain   [-]   (nvoigt, n_mp)
      real(wp), allocatable :: eps_p(:,:)  !! plastic strain [-]   (nvoigt, n_mp)
      real(wp), allocatable :: statev(:,:) !! constitutive state   (nstatev, n_mp)

      ! --- fluid (allocated only when phase config requests them) ---
      real(wp), allocatable :: uw(:)  !! pore water pressure [kPa] (n_mp)
      real(wp), allocatable :: ua(:)  !! pore air pressure   [kPa] (n_mp)
      real(wp), allocatable :: Sr(:)  !! degree of saturation [-]  (n_mp)
   end type mpm_particle_set_t

contains

   subroutine allocate_particle_set(set, n_mp, problem_type, nstatev, config)
      !! Allocate all arrays in a particle set and initialise to zero.
      !!
      !! Fluid arrays are only allocated when the corresponding flag is set in
      !! config. If config is not provided, only solid fields are allocated.
      type(mpm_particle_set_t),          intent(out) :: set
      integer,                           intent(in)  :: n_mp         !! number of material points
      integer,                           intent(in)  :: problem_type !! PROB_1D, PROB_2D_*, PROB_3D
      integer,                           intent(in)  :: nstatev      !! constitutive state variables
      type(mpm_phase_config_t), optional, intent(in) :: config

      integer :: nd, nv

      nd = n_dims(problem_type)
      nv = n_voigt(problem_type)

      set%n_mp    = n_mp
      set%ndim    = nd
      set%nvoigt  = nv
      set%nstatev = nstatev

      ! kinematics
      allocate(set%xp(nd, n_mp),   source=0.0_wp)
      allocate(set%vp(nd, n_mp),   source=0.0_wp)
      allocate(set%ap(nd, n_mp),   source=0.0_wp)
      allocate(set%mp(n_mp),       source=0.0_wp)
      allocate(set%vol(n_mp),      source=0.0_wp)
      allocate(set%rho(n_mp),      source=0.0_wp)

      ! solid mechanics
      allocate(set%sig(nv, n_mp),      source=0.0_wp)
      allocate(set%eps(nv, n_mp),      source=0.0_wp)
      allocate(set%eps_p(nv, n_mp),    source=0.0_wp)
      allocate(set%statev(nstatev, n_mp), source=0.0_wp)

      ! fluid (conditional)
      if (present(config)) then
         if (config%has_liquid) then
            allocate(set%uw(n_mp), source=0.0_wp)
            allocate(set%Sr(n_mp), source=1.0_wp)   ! fully saturated default
         end if
         if (config%has_gas) then
            allocate(set%ua(n_mp), source=0.0_wp)
         end if
      end if

   end subroutine allocate_particle_set

   subroutine deallocate_particle_set(set)
      !! Deallocate all arrays in a particle set.
      type(mpm_particle_set_t), intent(inout) :: set

      if (allocated(set%xp))     deallocate(set%xp)
      if (allocated(set%vp))     deallocate(set%vp)
      if (allocated(set%ap))     deallocate(set%ap)
      if (allocated(set%mp))     deallocate(set%mp)
      if (allocated(set%vol))    deallocate(set%vol)
      if (allocated(set%rho))    deallocate(set%rho)
      if (allocated(set%sig))    deallocate(set%sig)
      if (allocated(set%eps))    deallocate(set%eps)
      if (allocated(set%eps_p))  deallocate(set%eps_p)
      if (allocated(set%statev)) deallocate(set%statev)
      if (allocated(set%uw))     deallocate(set%uw)
      if (allocated(set%ua))     deallocate(set%ua)
      if (allocated(set%Sr))     deallocate(set%Sr)

      set%n_mp    = 0
      set%ndim    = 0
      set%nvoigt  = 0
      set%nstatev = 0

   end subroutine deallocate_particle_set

end module mpm_particles
