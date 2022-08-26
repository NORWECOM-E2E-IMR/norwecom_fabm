#include "fabm_driver.h"

module imr_norwecom

    use fabm_types

    implicit none

    private

    logical, parameter :: zooplankton = .false. !! Include or exclude zooplankton

    !> NORWECOM FABM model type definition
    type, extends(type_base_model), public :: type_imr_norwecom
        ! Define state variables
        type(type_state_variable_id) :: id_nit !! Nitrate
        type(type_state_variable_id) :: id_pho !! Phosphate
        type(type_state_variable_id) :: id_sil !! Silicate
        type(type_state_variable_id) :: id_sis !! Biogenic silica
        type(type_state_variable_id) :: id_det !! Nitrogen detritus
        type(type_state_variable_id) :: id_detp !! Phosphorus detritus
        type(type_state_variable_id) :: id_oxy !! Dissolved oxygen
        type(type_state_variable_id) :: id_dia !! Diatoms
        type(type_state_variable_id) :: id_fla !! Flagellates

        ! Define diagnostic variables
        type(type_diagnostic_variable_id) :: id_chla !! Chlorophyll a
        type(type_diagnostic_variable_id) :: id_gpp !! Gross primary production
        type(type_diagnostic_variable_id) :: id_npp !! Net primary production
        type(type_diagnostic_variable_id) :: id_totsil !! Total silicate concentration

        ! Define dependencies
        type(type_dependency_id) :: id_temp !! Temperature
        type(type_dependency_id) :: id_salt !! Salinity
        type(type_dependency_id) :: id_dens !! Density
        type(type_dependency_id) :: id_par !! Photoactive radition

        ! Define parameters
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
        real(rk) :: sib !! Si concentration with diatom max sinking rate
        real(rk) :: scc1 !! O/N consumption rate
        real(rk) :: scc2 !! Intercellular C/N ratio
        real(rk) :: scc4 !! Bioenegic silica decomposition rate

    contains
        procedure :: initialize
        procedure :: do_surface
        procedure :: do
        procedure :: get_vertical_movement
    end type

contains

    !> Initialize the NORWECOM FABM model
    !!
    !! Register state variables, diagnostic variables and paramters and returns the 
    !! initialized model.
    !!
    subroutine initialize(self, configunit)
        class(type_imr_norwecom), intent(inout), target :: self !! NORWECOM FABM model class
        integer, intent(in) :: configunit

        real(rk) :: pmax, dia_ar, dia_an, dia_ap, dia_as, fla_ar, fla_an, fla_ap

        ! Initialize state variables
        call self%register_state_variable(self%id_nit, "nit", "mgN m-3", "Nitrate concentration", &
            minimum = 0.0_rk, initial_value = 168.0_rk)
        call self%register_state_variable(self%id_pho, "pho", "mgP m-3", "Phosphate concentration", &
            minimum = 0.0_rk, initial_value = 25.0_rk)
        call self%register_state_variable(self%id_sil, "sil", "mgSi m-3", "Silicate concentraton", &
            minimum = 0.0_rk, initial_value = 155.0_rk)
        call self%register_state_variable(self%id_sis, "sis", "mgSi m-3", "Biogenic silica concentration", &
            minimum = 0.0_rk, initial_value = 3.0_rk, vertical_movement = -3.47e-5_rk)
        call self%register_state_variable(self%id_det, "det", "mgN m-3", "Nitrogen detritus concentration", &
            minimum = 0.0_rk, initial_value = 0.1_rk, vertical_movement = -3.47e-5_rk)
        call self%register_state_variable(self%id_detp, "detp", "mgP m-3", "Phosphorus detritus concentration", &
            minimum = 0.0_rk, initial_value = 0.1_rk, vertical_movement = -3.45e-5_rk)
        call self%register_state_variable(self%id_oxy, "oxy", "mg l-1", "Dissolved oxygen concentration", &
            minimum = 0.0_rk, initial_value = 10.0_rk)
        call self%register_state_variable(self%id_dia, "dia", "mgN m-3", "Diatoms concentration", &
            minimum = 0.0001_rk, initial_value = 0.1_rk)
        call self%register_state_variable(self%id_fla, "fla", "mgN m-3", "Flagellates concentration", &
            minimum = 0.0001_rk, initial_value = 0.1_rk, vertical_movement = -2.89e-6_rk)

        !----- Initialize diagnostic variables -----!
        call self%register_diagnostic_variable(self%id_chla, "chla", "mgChla m-3", "Chlorophyll a concentration")
        call self%register_diagnostic_variable(self%id_gpp, "gpp", "mgC m-3", "Gross primary production")
        call self%register_diagnostic_variable(self%id_npp, "npp", "mgC m-3", "Net primary production")
        call self%register_diagnostic_variable(self%id_totsil, "totsis", "mgSi m-3", "Total silicate concentration")

        !----- Initialize dependencies -----!
        call self%register_dependency(self%id_temp, standard_variables%temperature)
        call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
        call self%register_dependency(self%id_dens, standard_variables%density)
        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

        !----- Initialize parameters -----!

        ! Atomic weights
        call self%get_parameter(self%cnit, "cnit", "mg mmol-1", "Nitrogen atomic weight", default = 14.01_rk)
        call self%get_parameter(self%cpho, "cpho", "mg mmol-1", "Phosphorus atomic weight", default = 30.97_rk)
        call self%get_parameter(self%csil, "csil", "mg mmol-1", "Silicium atomic weight", default = 28.09_rk)

        ! Phytoplankton growth
        call self%get_parameter(self%a1, "a1", "s-1", "Diatoms pmax at 0 degC", default = 1.53e-5_rk)
        call self%get_parameter(self%a2, "a2", "degC-1", "Diatoms pmax temperature dependence", default = 0.063_rk)
        call self%get_parameter(self%a3, "a3", "s-1", "Flagellates pmax at 0 degC", default = 1.02e-5_rk)
        call self%get_parameter(self%a4, "a4", "degC-1", "Flagellates pmax temperature dependence", default = 0.063_rk)
        call self%get_parameter(self%a5, "a5", "s-1", "Phytoplankton respiration loss at 0 degC", default = 8.05e-7_rk)
        call self%get_parameter(self%a6, "a6", "degC-1", "Phytoplankton respiration loss temperature dependence", default = 0.07_rk)
        call self%get_parameter(self%dia_kr, "dia_kr", "m2 uE-1", "Diatoms half-saturation constant for light uptake", default = 96.0_rk)
        call self%get_parameter(self%dia_kn, "dia_kn", "mmol m-3", "Diatoms half-saturation constant for nitrate uptake", default = 2.0_rk)
        call self%get_parameter(self%dia_kp, "dia_kp", "s-1 uM-1", "Diatoms half-saturation constant for phosphate uptake", default = 0.125_rk)
        call self%get_parameter(self%dia_ks, "dia_ks", "s-1 uM-1", "Diatoms half-saturation constant for silicate uptake", default = 1.4_rk)
        call self%get_parameter(self%fla_kr, "fla_kr", "m2 uE-1", "Flagellates half-saturation constant for light uptake", default = 209.0_rk)
        call self%get_parameter(self%fla_kn, "fla_kn", "s-1 uM-1", "Flagellates half-saturation constant for nitrate uptake", default = 1.5_rk)
        call self%get_parameter(self%fla_kp, "fla_kp", "s-1 uM-1", "Flagellates half-saturation constant for phosphate uptake", default = 0.094_rk)
        call self%get_parameter(self%tmean, "tmean", "degC", "Reference temperature for substrate affinity", default = 13.0_rk)
        call get_affinities(dia_ar, dia_an, dia_ap, dia_as, fla_ar, fla_an, fla_ap)
        call self%get_parameter(self%dia_ar, "dia_ar", "m2 uE-1", "Diatoms affinity for light", default = dia_ar)
        call self%get_parameter(self%dia_an, "dia_an", "s-1 (mgN m-3)-1", "Diatoms affinity for nitrate", default = dia_an)
        call self%get_parameter(self%dia_ap, "dia_ap", "s-1 (mgP m-3)-1", "Diatoms affinity for phosphate", default = dia_ap)
        call self%get_parameter(self%dia_as, "dia_as", "s-1 (mgSi m-3)-1", "Diatoms affinity for silicate", default = dia_as)
        call self%get_parameter(self%fla_ar, "fla_ar", "m2 uE-1", "Flagellates affinity for light", default = fla_ar)
        call self%get_parameter(self%fla_an, "fla_an", "s-1 (mgN m-3)-1", "Flagellates affinity for nitrate", default = fla_an)
        call self%get_parameter(self%fla_ap, "fla_ap", "s-1 (mgP m-3)-1", "Flagellates affinity for phosphate", default = fla_ap)
        print *, self%dia_ar, self%dia_an, self%dia_ap, self%dia_as, self%fla_ar, self%fla_an, self%fla_ap
        call self%get_parameter(self%cc1, "cc1", "mgP mgN-1", "Intercellular P/N ratio", default = 0.138_rk)
        call self%get_parameter(self%cc2, "cc2", "mgSi mgN-1", "Intercellular Si/N ratio", default = 1.75_rk)
        if (.not. zooplankton) then
            call self%get_parameter(self%cc3, "cc3", "s-1", "Phytoplankton mortality rate", default = 1.6e-6_rk)
        end if
        call self%get_parameter(self%cc4, "cc4", "s-1", "Detritus decomposition rate", default = 1.52e-7_rk)
        call self%get_parameter(self%diamin, "diamin", "mgN m-3", "Minimum diatoms concentration", default = 0.1_rk)
        call self%get_parameter(self%flamin, "flamin", "mgN m-3", "Minimum flagellates concentration", default = 0.1_rk)
        call self%get_parameter(self%n2chla, "n2chla", "mgN mgChla-1", "Cellular fraction of nitrate and chlorophyll a", default = 11.0_rk)
        call self%get_parameter(self%srdia_min, "srdia_min", "mgN m-3", "Minimum diatoms sinking rate", default = 3.47e-6_rk)
        call self%get_parameter(self%srdia_max, "srdia_max", "mgN m-3", "Maximum diatoms sinking rate", default = 3.47e-5_rk)
        call self%get_parameter(self%sib, "sib", "uM", "Si concentration with diatom max sinking rate", default = 1.0_rk)
        call self%get_parameter(self%scc1, "scc1", "mgO mgN-1", "O/N consumption rate", default = 19.71_rk)
        call self%get_parameter(self%scc2, "scc2", "mgC mgN-1", "Intercellular C/N ratio", default = 5.68_rk)
        call self%get_parameter(self%scc4, "scc4", "s-1", "Biogenic silica decomposition rate", default = 6.41e-8_rk)

        ! Initiate aggregate variables
        call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_nit)
        call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_det)
        call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_dia)
        call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_fla)

    contains
        
        !> Calculate substrate affinities for diatoms and flagellates
        subroutine get_affinities(dia_ar, dia_an, dia_ap, dia_as, fla_ar, fla_an, fla_ap)
            real(rk), intent(inout) :: dia_ar, dia_an, dia_ap, dia_as !! Affinities for diatoms
            real(rk), intent(inout) :: fla_ar, fla_an, fla_ap !! Affinities for flagellates

            ! Diatoms
            pmax = self%a1 * exp(self%tmean * self%a2)
            dia_ar = pmax / self%dia_kr
            dia_an = pmax / (self%cnit * self%dia_kn)
            dia_ap = pmax / (self%cpho * self%dia_kp)
            dia_as = pmax / (self%csil * self%dia_ks)

            ! Flagellates
            pmax = self%a3 * exp(self%tmean * self%a4)
            fla_ar = pmax / self%fla_kr
            fla_an = pmax / (self%cnit * self%fla_kn)
            fla_ap = pmax / (self%cpho * self%fla_kp)
        end subroutine get_affinities

    end subroutine initialize

    !> Runs surface-water exchange processes
    !!
    !! Currently only the oxygen air-sea flux
    subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
        class(type_imr_norwecom), intent(in) :: self !! NORWECOM FABM model class
        _DECLARE_ARGUMENTS_DO_SURFACE_

        ! Temporary variables
        real(rk) :: temp, salt, dens, oxy
        real(rk) :: tempk, pvel, osat, doxy
        pvel = 5.0_rk

        _SURFACE_LOOP_BEGIN_

        _GET_(self%id_temp, temp)
        _GET_(self%id_salt, salt)
        _GET_(self%id_dens, dens)
        _GET_(self%id_oxy, oxy)

        tempk = (temp + 273.15) / 100.0
        osat = exp(-173.9894 + 255.5907 / tempk + 146.4813 * log(tempk) - 22.2040 * tempk + &
            salt * (-0.037376 + 0.016504 * tempk - 0.0020564 * tempk * tempk)) ! mol kg-1
        osat = osat * dens ! mol m-3
        osat = osat * 32.0 * 1e-6 ! mg m-3
        doxy = (pvel / 86400.0) * (osat - oxy)

        _ADD_SURFACE_FLUX_(self%id_oxy, doxy)

        _SURFACE_LOOP_END_

    end subroutine do_surface

    !> Runs pelagic processes of the biogeochemical model
    subroutine do(self, _ARGUMENTS_DO_)
        class(type_imr_norwecom), intent(in) :: self !! NORWECOM FABM model class
        _DECLARE_ARGUMENTS_DO_

        ! Temporary variables
        real(rk) :: temp, nit, pho, sil, sis, det, detp, oxy, dia, fla, par, rad
        real(rk) :: umax, rad_lim, nit_lim, pho_lim, sil_lim
        real(rk) :: prod_dia, resp_dia, mort_dia, prod_fla, resp_fla, mort_fla
        real(rk) :: dnit, dpho, dsil, dsis, ddet, doxy, ddia, dfla
        real(rk) :: gpp, npp, chla

        _LOOP_BEGIN_

        ! Get local variables
        _GET_(self%id_temp, temp)
        _GET_(self%id_par, par)
        _GET_(self%id_nit, nit)
        _GET_(self%id_pho, pho)
        _GET_(self%id_sil, sil)
        _GET_(self%id_sis, sis)
        _GET_(self%id_det, det)
        _GET_(self%id_detp, detp)
        _GET_(self%id_oxy, oxy)
        _GET_(self%id_dia, dia)
        _GET_(self%id_fla, fla)

        ! Convert watts to micro einstein
        rad = par / 0.217_rk ! W m-2 -> uE m-2 s-1

        !----- Phytoplankton terms -----!
        
        ! Diatoms growth
        umax = self%a1 * exp(self%a2 * temp)
        rad_lim = slim(self%dia_ar, umax, rad)
        nit_lim = slim(self%dia_an, umax, nit)
        pho_lim = slim(self%dia_ap, umax, pho)
        sil_lim = slim(self%dia_as, umax, sil)
        prod_dia = umax * min(rad_lim, nit_lim, pho_lim, sil_lim) * dia
        resp_dia = self%a5 * dia * exp(self%a6 * temp)
        mort_dia = self%cc3 * dia
        
        ! Constrain diatom losses
        if (dia < self%diamin) then
            resp_dia = 0.0_rk
            mort_dia = 0.0_rk
        end if

        ! Flagellates growth
        umax = self%a3 * exp(self%a4 * temp)
        rad_lim = slim(self%fla_ar, umax, rad)
        nit_lim = slim(self%fla_an, umax, nit)
        pho_lim = slim(self%fla_ap, umax, pho)
        prod_fla = umax * min(rad_lim, nit_lim, pho_lim) * fla
        resp_fla = self%a5 * fla * exp(self%a6 * temp)
        mort_fla = self%cc3 * fla

        ! Constrain flagellate losses
        if (fla < self%flamin) then
            resp_fla = 0.0_rk
            mort_fla = 0.0_rk
        end if

        ! Primary production
        gpp = self%scc2 * (prod_dia + prod_fla)
        npp = self%scc2 * (prod_dia + prod_fla - (resp_dia + resp_fla))
        chla = (dia + fla) / self%n2chla
     
        !----- Fluxes -----!
        dnit = resp_dia + resp_fla + self%cc4 * det - (prod_dia + prod_fla)
        dpho = self%cc1 * dnit
        dsil = self%scc4 * sis - self%cc2 * prod_dia
        dsis = self%cc2 * (resp_dia + mort_dia) - self%scc4 * sis
        ddet = mort_dia + mort_fla - self%cc4 * det
        doxy = (self%scc1 * (prod_dia + prod_fla - (resp_dia + resp_fla + self%cc4 * det))) * 1e-3 ! mg m-3 -> mg l
        ddia = prod_dia - (resp_dia + mort_dia)
        dfla = prod_fla - (resp_fla + mort_fla)

        !----- Update FABM -----!
        _ADD_SOURCE_(self%id_nit, dnit)
        _ADD_SOURCE_(self%id_pho, dpho)
        _ADD_SOURCE_(self%id_sil, dsil)
        _ADD_SOURCE_(self%id_sis, dsis)
        _ADD_SOURCE_(self%id_det, ddet)
        _ADD_SOURCE_(self%id_oxy, doxy)
        _ADD_SOURCE_(self%id_dia, ddia)
        _ADD_SOURCE_(self%id_fla, dfla)

        _SET_DIAGNOSTIC_(self%id_totsil, sil+sis)
        _SET_DIAGNOSTIC_(self%id_gpp, gpp)
        _SET_DIAGNOSTIC_(self%id_npp, npp)
        _SET_DIAGNOSTIC_(self%id_chla, chla)

        _LOOP_END_

    contains

        !> Calculates the substrate limitation term
        real(rk) function slim(alpha, pmax, s)
            real(rk), intent(in) :: alpha !! Substrate affinity
            real(rk), intent(in) :: pmax !! Production max
            real(rk), intent(in) :: s !! Substrate concentration
            if (s < 0.0) then
                slim = 0.0
            else
                slim = s / (s + (pmax / alpha))
            end if
        end function slim

    end subroutine do



    !> Sets the sinking speed of diatoms
    subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
        class(type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

        real(rk) :: dia, sil, vdia, csil

        csil = 28.09_rk

        _LOOP_BEGIN_

        _GET_(self%id_dia, dia)
        _GET_(self%id_sil, sil)

        if (sil / csil < self%sib) then
            vdia = self%srdia_max
        else
            vdia = self%srdia_min + (self%srdia_max - self%srdia_min) / (sil / csil)
        end if
        
        ! Stop sinking if concentration is too low
        if (dia < self%diamin) vdia = 0.0_rk

        _ADD_VERTICAL_VELOCITY_(self%id_dia, -1.0 * vdia)

        _LOOP_END_
    end subroutine get_vertical_movement

end module imr_norwecom