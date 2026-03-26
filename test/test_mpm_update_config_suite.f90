module test_mpm_update_config_suite
   !! Tests for mpm_update_config: defaults, scheme constants, validation.
   use testdrive,         only: new_unittest, unittest_type, error_type, check
   use mpm_precision,     only: wp
   use mpm_update_config, only: mpm_update_config_t, &
                                 SCHEME_FLIP_PIC, SCHEME_APIC, SCHEME_XPIC, &
                                 is_valid_config
   implicit none
   private

   public :: collect_mpm_update_config_suite

contains

   subroutine collect_mpm_update_config_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("default_scheme_is_flip_pic",   test_default_scheme),   &
         new_unittest("default_alpha_is_0p99",        test_default_alpha),    &
         new_unittest("scheme_constants_are_distinct", test_schemes_distinct), &
         new_unittest("alpha_1p0_is_valid",           test_alpha_1p0_valid),  &
         new_unittest("alpha_0p0_is_valid",           test_alpha_0p0_valid),  &
         new_unittest("alpha_above_1_is_invalid",     test_alpha_too_high),   &
         new_unittest("alpha_below_0_is_invalid",     test_alpha_too_low)     &
      ]
   end subroutine collect_mpm_update_config_suite

   subroutine test_default_scheme(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_update_config_t) :: config
      call check(error, config%scheme == SCHEME_FLIP_PIC)
   end subroutine

   subroutine test_default_alpha(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_update_config_t) :: config
      call check(error, config%alpha == 0.99_wp)
   end subroutine

   subroutine test_schemes_distinct(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, SCHEME_FLIP_PIC /= SCHEME_APIC);  if (allocated(error)) return
      call check(error, SCHEME_FLIP_PIC /= SCHEME_XPIC);  if (allocated(error)) return
      call check(error, SCHEME_APIC     /= SCHEME_XPIC)
   end subroutine

   subroutine test_alpha_1p0_valid(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_update_config_t) :: config
      config%alpha = 1.0_wp
      call check(error, is_valid_config(config))
   end subroutine

   subroutine test_alpha_0p0_valid(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_update_config_t) :: config
      config%alpha = 0.0_wp
      call check(error, is_valid_config(config))
   end subroutine

   subroutine test_alpha_too_high(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_update_config_t) :: config
      config%alpha = 1.0001_wp
      call check(error, .not. is_valid_config(config))
   end subroutine

   subroutine test_alpha_too_low(error)
      type(error_type), allocatable, intent(out) :: error
      type(mpm_update_config_t) :: config
      config%alpha = -0.0001_wp
      call check(error, .not. is_valid_config(config))
   end subroutine

end module test_mpm_update_config_suite
