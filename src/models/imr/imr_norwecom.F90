#include "fabm_driver.h"

module imr_norwecom

    use fabm_types

    implicit none

    private

    type, extends(type_base_model), public :: type_imr_norwecom

    contains
        procedure :: initialize
    end type

contains

    subroutine initialize(self, configunit)
        class(type_imr_norwecom), intent(inout), target :: self
        integer, intent(in) :: configunit
    end subroutine initialize

end module imr_norwecom