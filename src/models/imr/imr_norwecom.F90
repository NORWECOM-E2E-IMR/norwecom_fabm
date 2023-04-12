#include "fabm_driver.h"

module imr_norwecom
    !! Module for declaring the FABM implementation of NORWECOM and associated
    !! procedures.

    use fabm_types

    implicit none

    private

    ! Global parameters
    logical, parameter :: zooplankton = .true. !! Include or exclude zooplankton
    real(rk), parameter :: eps = 1.0e-6_rk
    real(rk), parameter :: day_sec = 86400.0_rk

    
    type, extends(type_base_model), public :: type_imr_norwecom
        !! NORWECOM FABM model type definition

        ! State variables
        type(type_state_variable_id) :: id_nit !! Nitrate
        type(type_state_variable_id) :: id_pho !! Phosphate
        type(type_state_variable_id) :: id_sil !! Silicate
        type(type_state_variable_id) :: id_sis !! Biogenic silica
        type(type_state_variable_id) :: id_det !! Nitrogen detritus
        type(type_state_variable_id) :: id_detp !! Phosphorus detritus
        type(type_state_variable_id) :: id_oxy !! Dissolved oxygen
        type(type_state_variable_id) :: id_dia !! Diatoms
        type(type_state_variable_id) :: id_fla !! Flagellates
        type(type_state_variable_id) :: id_mic !! Microzooplankton
        type(type_state_variable_id) :: id_mes !! Mesozooplankton

        ! Diagnostic variables
        type(type_diagnostic_variable_id) :: id_chla !! Chlorophyll a
        type(type_diagnostic_variable_id) :: id_gpp !! Gross primary production
        type(type_diagnostic_variable_id) :: id_npp !! Net primary production
        type(type_diagnostic_variable_id) :: id_gsp !! Gross secondary production
        type(type_diagnostic_variable_id) :: id_nsp !! Net secondary production
        type(type_diagnostic_variable_id) :: id_totsil !! Total silicate concentration
        type(type_diagnostic_variable_id) :: id_totpho !! Total phosphate concentration
        type(type_diagnostic_variable_id) :: id_radlim !! Light limitation
        type(type_diagnostic_variable_id) :: id_nitlim !! Nitrate limitation
        type(type_diagnostic_variable_id) :: id_pholim !! Phosphate limitation
        type(type_diagnostic_variable_id) :: id_sillim !! Silicate limitation

        ! Dependencies
        type(type_dependency_id) :: id_temp !! Temperature
        type(type_dependency_id) :: id_salt !! Salinity
        type(type_dependency_id) :: id_dens !! Density
        type(type_dependency_id) :: id_par !! Photoactive radition

        ! Parameters
        real(rk) :: cnit !! Nitrogen atomic weight
        real(rk) :: cpho !! Phosphorus atomic weight
        real(rk) :: csil !! Silicium atomic weight
        real(rk) :: a1 !! Diatoms pmax at 0 degC
        real(rk) :: a2 !! Diatoms pmax temperature dependence
        real(rk) :: a3 !! Flagellates pmax at 0 degC
        real(rk) :: a4 !! Flagellates pmax temperature dependence
        real(rk) :: a5 !! Phytoplankton respiration at 0 degC
        real(rk) :: a6 !! Phytoplankton respiration temperature dependence
        real(rk) :: dia_kr !! Diatoms half-saturation constant for light "uptake" 
        real(rk) :: dia_kn !! Diatoms half-saturation constant for nitrate uptake
        real(rk) :: dia_kp !! Diatoms half-saturation constant for phosphate uptake
        real(rk) :: dia_ks !! Diatoms half-saturation constant for silicate uptake
        real(rk) :: fla_kr !! Flagellates half-saturation constant for light "uptake"
        real(rk) :: fla_kn !! Flagellates half-saturation constant for nitrate uptake
        real(rk) :: fla_kp !! Flagellates half-saturation constant for phosphate uptake
        real(rk) :: tmean !! Reference temperature for substrate affinity
        real(rk) :: dia_ar !! Diatoms affinity for light
        real(rk) :: dia_an !! Diatoms affinity for nitrate
        real(rk) :: dia_ap !! Diatoms affinity for phosphate
        real(rk) :: dia_as !! Diatoms affinity for silicate
        real(rk) :: fla_ar !! Flagellates affinity for light
        real(rk) :: fla_an !! Flagellates affinity for nitrate
        real(rk) :: fla_ap !! Flagellates affinity for phosphate
        real(rk) :: cc1 !! Intercellular P/N ratio
        real(rk) :: cc2 !! Intercellular Si/N ratio
        real(rk) :: cc3 !! Phytoplankton mortality rate
        real(rk) :: cc4 !! Detritus decomposition rate
        real(rk) :: diamin !! Mimimum diatoms concentration
        real(rk) :: flamin !! Mimimum flagellates concentration
        real(rk) :: n2chla !! Cellular fraction of nitrate and chlorophyll a
        real(rk) :: srdia_min !! Minimum diatoms sinking rate
        real(rk) :: srdia_max !! Maximum diatoms sinking rate
        real(rk) :: srdet !! Detritus sinking rate
        real(rk) :: sib !! Si concentration with diatom max sinking rate
        real(rk) :: scc1 !! O/N consumption rate
        real(rk) :: scc2 !! Intercellular C/N ratio
        real(rk) :: scc4 !! Bioenegic silica decomposition rate
        real(rk) :: pi11 !! Mesozooplankton prey preference for diatoms
        real(rk) :: pi12 !! Mesozooplankton prey preference for microzooplankton
        real(rk) :: pi13 !! Mesozooplankton prey preference for detritus
        real(rk) :: pi21 !! Microzooplankton prey preference for flagellates
        real(rk) :: pi22 !! Microzooplankton prey preference for detritus
        real(rk) :: k3 !! Half-saturation constant for zooplankton ingestion
        real(rk) :: k6 !! Half-saturation constant for zooplankton loss
        real(rk) :: q10 !! Temperature dependence on zooplankton growth
        real(rk) :: mju2 !! Maximum loss rate of zooplankton
        real(rk) :: beta !! Assimilation efficiency for zooplankton
        real(rk) :: delta !! Fraction of zooplankton losses to detritus
        real(rk) :: eps !! Fraction of zooplankton losses to nitrate
        real(rk) :: mes_g !! Mesozooplankton maximum growth rate
        real(rk) :: mic_g !! Microzooplankton maximum growth rate
        real(rk) :: v !! Chlorophyll a extinction coefficient
        real(rk) :: pvel !! Air-water oxygen exchange factor
    contains
        procedure :: initialize
        procedure :: do_surface
        procedure :: do
        procedure :: get_vertical_movement
        procedure :: get_light_extinction
    end type

    interface
        module subroutine initialize(self, configunit)
            class(type_imr_norwecom), intent(inout), target :: self
            integer, intent(in) :: configunit
        end subroutine initialize
        module subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
            class(type_imr_norwecom), intent(in) :: self
            _DECLARE_ARGUMENTS_DO_SURFACE_
        end subroutine do_surface
        module subroutine do(self, _ARGUMENTS_DO_)
            class(type_imr_norwecom), intent(in) :: self
            _DECLARE_ARGUMENTS_DO_
        end subroutine do
        module subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
            class(type_imr_norwecom), intent(in) :: self
            _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
        end subroutine get_vertical_movement
        module subroutine get_light_extinction(self, _ARGUMENTS_GET_EXTINCTION_)
            class (type_imr_norwecom), intent(in) :: self
            _DECLARE_ARGUMENTS_GET_EXTINCTION_
        end subroutine get_light_extinction
    end interface

end module imr_norwecom
