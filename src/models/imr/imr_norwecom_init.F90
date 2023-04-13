submodule (imr_norwecom) imr_norwecom_init
    !! Submodule for initializing the FABM implementation of NORWECOM

    implicit none

contains

    module subroutine initialize(self, configunit)
        !! Initialize the NORWECOM FABM model
        class(type_imr_norwecom), intent(inout), target :: self !! NORWECOM FABM model class
        integer, intent(in) :: configunit

        real(rk) :: pmax, dia_ar, dia_an, dia_ap, dia_as, fla_ar, fla_an, fla_ap

        !----- Initialize state variables -----!
        call self%register_state_variable(self%id_nit, "nit", "mgN m-3", "Nitrate concentration", &
            minimum = 0.0_rk, initial_value = 168.0_rk)
        call self%register_state_variable(self%id_pho, "pho", "mgP m-3", "Phosphate concentration", &
            minimum = 0.0_rk, initial_value = 25.0_rk)
        call self%register_state_variable(self%id_sil, "sil", "mgSi m-3", "Silicate concentraton", &
            minimum = 0.0_rk, initial_value = 155.0_rk)
        call self%register_state_variable(self%id_sis, "sis", "mgSi m-3", "Biogenic silica concentration", &
            minimum = 0.0_rk, initial_value = 3.0_rk, vertical_movement = -3.47e-5_rk)
        call self%register_state_variable(self%id_det, "det", "mgN m-3", "Nitrogen detritus concentration", &
            minimum = 0.0_rk, initial_value = 0.1_rk)
        call self%register_state_variable(self%id_detp, "detp", "mgP m-3", "Phosphorus detritus concentration", &
            minimum = 0.0_rk, initial_value = 0.1_rk, vertical_movement = -3.47e-5_rk)
        call self%register_state_variable(self%id_oxy, "oxy", "mg m-3", "Dissolved oxygen concentration", &
            minimum = 0.0_rk, initial_value = 10.0_rk)
        call self%register_state_variable(self%id_dia, "dia", "mgN m-3", "Diatoms concentration", &
            minimum = 0.0001_rk, initial_value = 0.1_rk)
        call self%register_state_variable(self%id_fla, "fla", "mgN m-3", "Flagellates concentration", &
            minimum = 0.0001_rk, initial_value = 0.1_rk, vertical_movement = -2.89e-6_rk)
        call self%register_state_variable(self%id_mic, "mic", "mgN m-3", "Microzooplankton concentraton", &
            minimum = 0.0001_rk, initial_value = 0.1_rk)
        call self%register_state_variable(self%id_mes, "mes", "mgN m-3", "Mesozooplankton concentration", &
            minimum = 0.0001_rk, initial_value = 0.1_rk)
        call self%register_bottom_state_variable(self%id_botdet, "botdet", "mgN m-2", "Bottom nitrogen detritus", &
            & minimum = 0.0001_rk, initial_value = 0.1_rk)
        call self%register_bottom_state_variable(self%id_botdetp, "botdetp", "mgP m-2", "Bottom phosphorus detritus", &
            & minimum = 0.0001_rk, initial_value = 0.1_rk)
        call self%register_bottom_state_variable(self%id_botsis, "botsis", "mgSi m-2", "Bottom biogenic silica", &
            & minimum = 0.0001_rk, initial_value = 0.1_rk)
        call self%register_bottom_state_variable(self%id_burdet, "burdet", "mgN m-2", "Burried nitrogen detritus", &
            & minimum = 0.0_rk, initial_value = 0.0_rk)
        call self%register_bottom_state_variable(self%id_burdetp, "butdetp", "mgP m-2", "Burried phosphorus detritus", &
            & minimum = 0.0_rk, initial_value = 0.0_rk)
        call self%register_bottom_state_variable(self%id_bursis, "bursis", "mgSi m-2", "Burried biogenic silica", &
            & minimum = 0.0_rk, initial_value = 0.0_rk)

        !----- Initialize diagnostic variables -----!
        call self%register_diagnostic_variable(self%id_chla, "chla", "mgChla m-3", "Chlorophyll a concentration")
        call self%register_diagnostic_variable(self%id_gpp, "gpp", "mgC m-3 s-1", "Gross primary production")
        call self%register_diagnostic_variable(self%id_npp, "npp", "mgC m-3 s-1", "Net primary production")
        call self%register_diagnostic_variable(self%id_gsp, "gsp", "mgC m-3 s-1", "Gross secondary production")
        call self%register_diagnostic_variable(self%id_nsp, "nsp", "mgC m-3 s-1", "Net secondary production")
        call self%register_diagnostic_variable(self%id_totsil, "totsis", "mgSi m-3", "Total silicate concentration")
        call self%register_diagnostic_variable(self%id_totpho, "totpho", "mgP m-3", "Total phosphate concentration")
        call self%register_diagnostic_variable(self%id_radlim, "radlim", "[0-1]", "Light limitation")
        call self%register_diagnostic_variable(self%id_nitlim, "nitlim", "[0-1]", "Nitrate limitation")
        call self%register_diagnostic_variable(self%id_pholim, "pholim", "[0-1]", "Phosphate limitation")
        call self%register_diagnostic_variable(self%id_sillim, "sillim", "[0-1]", "Silicate limitation")

        !----- Initialize dependencies -----!
        call self%register_dependency(self%id_temp, standard_variables%temperature)
        call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
        call self%register_dependency(self%id_dens, standard_variables%density)
        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_horizontal_dependency(self%id_bstress, standard_variables%bottom_stress)

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
        if (zooplankton) then
            call self%get_parameter(self%cc3, "cc3", "s-1", "Phytoplankton mortality rate", default = 1.6e-7_rk)
        else
            call self%get_parameter(self%cc3, "cc3", "s-1", "Phytoplankton mortality rate", default = 1.6e-6_rk)
        end if
        call self%get_parameter(self%cc4, "cc4", "s-1", "Detritus decomposition rate", default = 1.52e-7_rk)
        call self%get_parameter(self%diamin, "diamin", "mgN m-3", "Minimum diatoms concentration", default = 0.1_rk)
        call self%get_parameter(self%flamin, "flamin", "mgN m-3", "Minimum flagellates concentration", default = 0.1_rk)
        call self%get_parameter(self%n2chla, "n2chla", "mgN mgChla-1", "Cellular fraction of nitrate and chlorophyll a", default = 11.0_rk)
        call self%get_parameter(self%srdia_min, "srdia_min", "mgN m-3", "Minimum diatoms sinking rate", default = 3.47e-6_rk)
        call self%get_parameter(self%srdia_max, "srdia_max", "mgN m-3", "Maximum diatoms sinking rate", default = 3.47e-5_rk)
        call self%get_parameter(self%sib, "sib", "uM", "Si concentration with diatom max sinking rate", default = 1.0_rk)
        call self%get_parameter(self%srdet, "srdet", "s-1", "Detritus sinking rate", default = -3.47e-5_rk)
        call self%get_parameter(self%scc1, "scc1", "mgO mgN-1", "O/N consumption rate", default = 19.71_rk)
        call self%get_parameter(self%scc2, "scc2", "mgC mgN-1", "Intercellular C/N ratio", default = 5.68_rk)
        call self%get_parameter(self%scc4, "scc4", "s-1", "Biogenic silica decomposition rate", default = 6.41e-8_rk)

        ! Light
        call self%get_parameter(self%v, "v", "m mgChla-1", "Chlorophyll a extinction coefficient", default = 1.38e-2_rk)
        call self%get_parameter(self%pvel, "pvel", "-", "Air-water oxygen exchange factor", default = 5.0_rk)

        ! Zooplankton
        call self%get_parameter(self%pi11, "pi11", "[0-1]", "Mesozooplankton prey preference for diatoms", default = 0.333_rk)
        call self%get_parameter(self%pi12, "pi12", "[0-1]", "Mesozooplankton prey preference for microzooplankton", default = 0.333_rk)
        call self%get_parameter(self%pi13, "pi13", "[0-1]", "Mesozooplankton prey preference for detritus", default = 0.333_rk)
        call self%get_parameter(self%pi21, "pi21", "[0-1]", "Microzooplankton prey preference for flagellates", default = 0.5_rk)
        call self%get_parameter(self%pi22, "pi22", "[0-1]", "Microzooplankton prey preference for detritus", default = 0.5_rk)
        call self%get_parameter(self%beta, "beta", "[0-1]", "Zooplankton assimilation efficiency", default = 0.75_rk)
        call self%get_parameter(self%mju2, "mju2", "d-1", "Maximum loss rate of zooplankton", default = 0.2_rk)
        call self%get_parameter(self%delta, "delta", "[0-1]", "Fraction of zooplankton losses to detritus", default = 0.6_rk)
        call self%get_parameter(self%eps, "eps", "[0-1]", "Fraction of zooplankton losses to nitrate", default = 0.4_rk)
        call self%get_parameter(self%mes_g, "mes_g", "d-1", "Mesozooplankton maximum growth rate", default = 0.4_rk)
        call self%get_parameter(self%mic_g, "mic_g", "d-1", "Microzooplankton maximum growth rate", default = 0.5_rk)
        call self%get_parameter(self%k3, "k3", "mmol m-3", "Half-saturation constant for zooplankton ingestion", default = 1.0_rk)
        call self%get_parameter(self%k6, "k6", "mmol m-3", "Half-saturation constant for zooplankton loss", default = 0.2_rk)
        call self%get_parameter(self%q10, "q10", "", "Temperature dependence on zooplankton growth", default = 1.5_rk)

        ! Sediment
        call self%get_parameter(self%tau1, "tau1", "Pa", "Bottom stress threshold for sedimentation", default = 0.064_rk)
        call self%get_parameter(self%tau2, "tau2", "Pa", "Bottom stress threshold for resuspension", default = 0.78_rk)
        call self%get_parameter(self%c2, "c2", "s m-1", "Slope of the linear increase in bottom flux", default = 100.0_rk)

        !----- Initialize aggregated variables -----!
        call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_nit)
        call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_det)
        call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_dia)
        call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_fla)

    contains
        
        subroutine get_affinities(dia_ar, dia_an, dia_ap, dia_as, fla_ar, fla_an, fla_ap)
            !! Calculate substrate affinities for diatoms and flagellates
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

end submodule imr_norwecom_init