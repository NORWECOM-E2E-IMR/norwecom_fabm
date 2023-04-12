#include "fabm_driver.h"

submodule (imr_norwecom) imr_norwecom_npzd
    !! Submodule for performing biogeochemical processes in the FABM implementation
    !! of the NORWECOM model.

    implicit none

contains

    module subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
        !! Perform surface-water exchange processes
        class(type_imr_norwecom), intent(in) :: self !! NORWECOM FABM model class
        _DECLARE_ARGUMENTS_DO_SURFACE_

        ! Temporary variables
        real(rk) :: temp, salt, dens, oxy
        real(rk) :: tempk, pvel, osat, doxy

        _SURFACE_LOOP_BEGIN_

        _GET_(self%id_temp, temp)
        _GET_(self%id_salt, salt)
        _GET_(self%id_dens, dens)
        _GET_(self%id_oxy, oxy)

        ! TODO: move parameters to model definition
        tempk = (temp + 273.15) / 100.0
        osat = exp(-173.9894 + 255.5907 / tempk + 146.4813 * log(tempk) - 22.2040 * tempk + &
            salt * (-0.037376 + 0.016504 * tempk - 0.0020564 * tempk * tempk)) ! mol kg-1
        osat = osat * dens ! mol m-3
        osat = osat * 32.0 * 1e-3_rk ! mg m-3
        doxy = (self%pvel / 86400.0) * (osat - oxy)

        _ADD_SURFACE_FLUX_(self%id_oxy, doxy)

        _SURFACE_LOOP_END_
    end subroutine do_surface

    module subroutine do(self, _ARGUMENTS_DO_)
        !! Perform pelagic processes of the biogeochemical model
        class(type_imr_norwecom), intent(in) :: self !! NORWECOM FABM model class
        _DECLARE_ARGUMENTS_DO_

        ! Temporary variables
        real(rk) :: temp, nit, pho, sil, sis, det, detp, oxy, dia, fla, par, rad
        real(rk) :: nit2dia, dia2nit, dia2det, nit2fla, fla2nit, fla2det
        real(rk) :: dnit, dpho, dsil, dsis, ddet, doxy, ddia, dfla, ddetp
        real(rk) :: gpp, npp, chla, gsp, nsp
        real(rk) :: mes, mic, dmes, dmic
        real(rk) :: g11, g12, g13, g21, g22
        real(rk) :: dia2mes, mic2mes, det2mes, mes2det, mes2nit
        real(rk) :: fla2mic, det2mic, mic2det, mic2nit

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
        _GET_(self%id_mes, mes)
        _GET_(self%id_mic, mic)

        ! Convert watts to micro einstein
        rad = par / 0.217_rk ! W m-2 -> uE m-2 s-1

        !----- Phytoplankton terms -----!
        call do_phytoplankton(self, temp, dia, fla, rad, nit, pho, sil, nit2dia, &
            & dia2nit, dia2det, nit2fla, fla2nit, fla2det, gpp, npp, chla)

        !----- Zooplankton terms -----!
        call do_zooplankton(self, temp, dia, fla, det, mic, mes, dia2mes, mic2mes, &
            & det2mes, fla2mic, det2mic, mes2det, mes2nit, mic2det, mic2nit, gsp, nsp)

        !----- Fluxes -----!
        dnit = dia2nit + fla2nit + 0.10_rk*(dia2det + fla2det) + self%cc4*det &
            & + mes2nit + mic2nit - (nit2dia + nit2fla)
        dpho = self%cc1*(dia2nit + fla2nit + 0.25_rk*(dia2det + fla2det) &
            & + mes2nit + mic2nit - (nit2dia + nit2fla)) + 1.3_rk*self%cc4*detp
        dsil = self%scc4*sis - self%cc2*nit2dia
        dsis = self%cc2*(dia2nit + dia2det + dia2mes) - self%scc4*sis
        ddet = 0.9_rk*(dia2det + fla2det) + mes2det + mic2det &
            & - (self%cc4*det + det2mes + det2mic)
        ddetp = self%cc1*(0.75_rk*(dia2det + fla2det) + mes2det + mic2det &
            & - (det2mes + det2mic)) - 1.3_rk*self%cc4*detp
        doxy = -1.0_rk*self%scc1*dnit
        ddia = nit2dia - (dia2nit + dia2det + dia2mes)
        dfla = nit2fla - (fla2nit + fla2det + fla2mic)
        dmes = dia2mes + mic2mes + det2mes - (mes2det + mes2nit)
        dmic = fla2mic + det2mic - (mic2mes + mic2det + mic2nit)

        !----- Update FABM -----!
        _ADD_SOURCE_(self%id_nit, dnit)
        _ADD_SOURCE_(self%id_pho, dpho)
        _ADD_SOURCE_(self%id_sil, dsil)
        _ADD_SOURCE_(self%id_sis, dsis)
        _ADD_SOURCE_(self%id_det, ddet)
        _ADD_SOURCE_(self%id_detp, ddetp)
        _ADD_SOURCE_(self%id_oxy, doxy)
        _ADD_SOURCE_(self%id_dia, ddia)
        _ADD_SOURCE_(self%id_fla, dfla)
        _ADD_SOURCE_(self%id_mes, dmes)
        _ADD_SOURCE_(self%id_mic, dmic)

        _SET_DIAGNOSTIC_(self%id_totsil, sil+sis)
        _SET_DIAGNOSTIC_(self%id_totpho, pho+detp)
        _SET_DIAGNOSTIC_(self%id_gpp, gpp)
        _SET_DIAGNOSTIC_(self%id_npp, npp)
        _SET_DIAGNOSTIC_(self%id_gsp, gsp)
        _SET_DIAGNOSTIC_(self%id_nsp, nsp)
        _SET_DIAGNOSTIC_(self%id_chla, chla)

        _LOOP_END_

    end subroutine do

    subroutine do_phytoplankton(self, temp, dia, fla, rad, nit, pho, sil, pdia, rdia, mdia, pfla, rfla, mfla, gpp, npp, chla)
        class(type_imr_norwecom), intent(in) :: self
        real(rk), intent(in) :: temp !! Temperature
        real(rk), intent(in) :: dia !! Diatom concentration
        real(rk), intent(in) :: fla !! Flagellate concentration
        real(rk), intent(in) :: rad !! Photoactive radiation
        real(rk), intent(in) :: nit !! Nitrate concentration
        real(rk), intent(in) :: pho !! Phosphate concentration
        real(rk), intent(in) :: sil !! Silicate concentration
        real(rk), intent(out) :: pdia !! Diatom production
        real(rk), intent(out) :: rdia !! Diatom respiration
        real(rk), intent(out) :: mdia !! Diatom background mortality
        real(rk), intent(out) :: pfla !! Flagellate production
        real(rk), intent(out) :: rfla !! Flagellate respiration
        real(rk), intent(out) :: mfla !! Flagellate mortality
        real(rk), intent(out) :: gpp !! Gross primary production
        real(rk), intent(out) :: npp !! Net primary production
        real(rk), intent(out) :: chla !! Chlorophyll a concentration

        ! Local variables
        real(rk) :: umax, rad_lim, nit_lim, pho_lim, sil_lim

        ! Diatoms
        umax = self%a1*exp(self%a2*temp)
        rad_lim = slim(self%dia_ar, umax, rad)
        nit_lim = slim(self%dia_an, umax, nit)
        pho_lim = slim(self%dia_ap, umax, pho)
        sil_lim = slim(self%dia_as, umax, sil)
        pdia = umax*min(rad_lim, nit_lim, pho_lim, sil_lim)*dia
        rdia = self%a5*dia*exp(self%a6*temp)
        mdia = self%cc3*dia

        ! Constrain phytoplankton losses
        if (dia < self%diamin) then
            rdia = 0.0_rk
            mdia = 0.0_rk
        end if

        ! Flagellates growth
        umax = self%a3*exp(self%a4*temp)
        rad_lim = slim(self%fla_ar, umax, rad)
        nit_lim = slim(self%fla_an, umax, nit)
        pho_lim = slim(self%fla_ap, umax, pho)
        pfla = umax*min(rad_lim, nit_lim, pho_lim)*fla
        rfla = self%a5*fla*exp(self%a6*temp)
        mfla = self%cc3*fla

        ! Constrain flagellate losses
        if (fla < self%flamin) then
            rfla = 0.0_rk
            mfla = 0.0_rk
        end if

        ! Primary production and chlorophyll a
        gpp = self%scc2*(pdia + pfla)
        npp = gpp - self%scc2*(rdia + rfla)
        chla = (dia + fla)/self%n2chla
    
    contains

        real(rk) function slim(alpha, pmax, s)
            !! Calculates the substrate limitation term
            real(rk), intent(in) :: alpha !! Substrate affinity
            real(rk), intent(in) :: pmax !! Production max
            real(rk), intent(in) :: s !! Substrate concentration
            if (s < 0.0) then
                slim = 0.0
            else
                slim = s/(s + (pmax/alpha))
            end if
        end function slim

    end subroutine do_phytoplankton

    subroutine do_zooplankton(self, temp, dia, fla, det, mic, mes, dia2mes, mic2mes, det2mes, fla2mic, det2mic, &
            mes2det, mes2nit, mic2det, mic2nit, gsp, nsp)
        !! Do zooplankton
        class(type_imr_norwecom), intent(in) :: self
        real(rk), intent(in) :: temp
        real(rk), intent(in) :: dia !! Diatom concentration
        real(rk), intent(in) :: fla !! Flagellates concentration
        real(rk), intent(in) :: det !! Nitrogen detritus concentration
        real(rk), intent(in) :: mic !! Microzooplankton concentration
        real(rk), intent(in) :: mes !! Mesozooplankton concentration
        real(rk), intent(out) :: dia2mes !! Grazing by mesozooplankton on diatoms [d-1]
        real(rk), intent(out) :: mic2mes !! Grazing by mesozooplankton on microzooplankton [d-1]
        real(rk), intent(out) :: det2mes !! Grazing by mesozooplankton on detritus [d-1]
        real(rk), intent(out) :: fla2mic !! Grazing by microzooplankton on flagellates [d-1]
        real(rk), intent(out) :: det2mic !! Grazing by microzooplankton on detritus [d-1]
        real(rk), intent(out) :: mes2det
        real(rk), intent(out) :: mes2nit
        real(rk), intent(out) :: mic2det
        real(rk), intent(out) :: mic2nit
        real(rk), intent(out) :: gsp
        real(rk), intent(out) :: nsp

        ! Local variables
        real :: denum, p11, p12, p13, p21, p22, tmp

        ! Mesozooplankton
        denum = self%pi11*dia + self%pi12*mic + self%pi13*det + eps
        p11 = self%pi11*dia/denum
        p12 = self%pi12*mic/denum
        p13 = self%pi13*det/denum
        tmp = (tfac(temp)*self%mes_g)/(self%cnit*self%k3 + p11*dia + p12*mic + p13*det)
        dia2mes = (tmp*p11*dia*mes)/day_sec
        mic2mes = (tmp*p12*mic*mes)/day_sec
        det2mes = (tmp*p13*det*mes)/day_sec

        ! Microzooplankton
        denum = self%pi21*fla + self%pi22*det + eps
        p21 = self%pi21*fla/denum
        p22 = self%pi22*det/denum
        tmp = (tfac(temp)*self%mic_g)/(self%cnit*self%k3 + p21*fla + p22*det)
        fla2mic = (tmp*p21*fla*mic)/day_sec
        det2mic = (tmp*p22*det*mic)/day_sec

        ! Zooplankton specific fluxes
        mes2det = (self%delta*self%mju2*(mes/(mes + self%cnit*self%k6))*mes + &
            (1.0_rk - self%beta)*(dia2mes + mic2mes + det2mes))/day_sec
        mes2nit = (self%eps*self%mju2*(mes/(mes + self%cnit*self%k6))*mes)/day_sec
        mic2det = (self%delta*self%mju2*(mic/(mic + self%cnit*self%k6))*mic + &
            (1.0_rk - self%beta)*(fla2mic + det2mic))/day_sec
        mic2nit = (self%eps*self%mju2*(mic/(mic + self%cnit*self%k6))*mic)/day_sec

        ! Zooplankton production
        gsp = self%scc2*self%beta*(dia2mes + mic2mes + det2mes + fla2mic + det2mic)
        nsp = gsp - self%scc2*(mes2nit + mic2nit)
    
    contains

        real(rk) function tfac(te)
            !! Calculates temperature dependence on zooplankton growth
            real(rk), intent(in) :: te !! Temperature [degC]
            tfac = self%q10**((te - 10.0_rk)/10.0_rk)
        end function tfac

    end subroutine do_zooplankton

    

end submodule imr_norwecom_npzd