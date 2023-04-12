#include "fabm_driver.h"

submodule (imr_norwecom) imr_norwecom_utils
    !! Submodule for utiliy procedures used in the FABM implementation of NORWECOM.

    implicit none

contains

    module subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
        !! Sets the sinking speed of diatoms
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
        _ADD_VERTICAL_VELOCITY_(self%id_det, self%srdet)

        _LOOP_END_
    end subroutine get_vertical_movement

    module subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
        !! Sets the light extinction coefficient
        class (type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_GET_EXTINCTION_

        real(rk) :: dia, fla
        real(rk) :: my_extinction
    
        _LOOP_BEGIN_
    
        _GET_(self%id_dia, dia)
        _GET_(self%id_fla, fla)

        my_extinction = self%v * ((dia + fla)/self%n2chla)

        _SET_EXTINCTION_( my_extinction )
    
        _LOOP_END_
    end subroutine get_light_extinction

end submodule imr_norwecom_utils