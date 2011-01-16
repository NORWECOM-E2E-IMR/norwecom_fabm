#define REALTYPE double precision
#define _ZERO_ 0.0d0
#define _ONE_  1.0d0

! Older Fortran compilers do not allow derived types to contain allocatable members
! (A Fortran >95 feature, defined in ISO Technical Report TR 15581 and part of the Fortran 2003 specification).
! As a workaround, they can be declared with the pointer attribute, which does bring a slight performance penalty.
! By using the below preprocessor definitions, the allocatable attribute is automatically replaced by the pointer
! attribute where needed, and related function calls are changed as well.
#ifdef _ISO_TR_15581_
#define _ALLOCATABLE_ allocatable
#define _NULL_ 
#define _ALLOCATED_ allocated
#else
#define _ALLOCATABLE_ pointer
#define _NULL_ =>null()
#define _ALLOCATED_ associated
#endif

! Data type for location variable(s)
#define _LOCATION_TYPE_ integer

! Define remaining dimensions in a 1D loop.
! NB if there is only one spatial dimension, there are no remaining dimensions!
#if _FABM_DIMENSION_COUNT_>1
#define _ARG_LOCATION_1DLOOP_ ,_LOCATION_1DLOOP_
#else
#define _ARG_LOCATION_1DLOOP_
#define _LOCATION_1DLOOP_
#define _VARIABLE_1DLOOP_ LOCATION
#endif

! Define dimension attribute and index specifyer for horizontal (2D) fields.
#ifdef _FABM_HORIZONTAL_IS_SCALAR_
#define _INDEX_LOCATION_HZ_
#define _ATTR_LOCATION_DIMENSIONS_HZ_
#define _ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE_ ,dimension(:)
#define _ARG_LOCATION_HZ_
#define _ARG_LOCATION_DIMENSIONS_HZ_
#else
#define _INDEX_LOCATION_HZ_ (LOCATION_HZ)
#define _ATTR_LOCATION_DIMENSIONS_HZ_ ,dimension(LOCATION_DIMENSIONS_HZ)
#define _ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE_ ,dimension(:,LOCATION_DIMENSIONS_HZ)
#define _ARG_LOCATION_HZ_ ,LOCATION_HZ
#define _ARG_LOCATION_DIMENSIONS_HZ_ ,LOCATION_DIMENSIONS_HZ
#endif

! Define dimension attribute and index specifyer for full 3D fields.
#define _INDEX_LOCATION_ (LOCATION)
#define _ATTR_LOCATION_DIMENSIONS_ ,dimension(LOCATION_DIMENSIONS)
#define _ATTR_LOCATION_DIMENSIONS_PLUS_ONE_ ,dimension(:,LOCATION_DIMENSIONS)
#define _ARG_LOCATION_ ,LOCATION
#define _ARG_LOCATION_DIMENSIONS_ ,:

#define _FABM_ARGS_0D_ environment _ARG_LOCATION_
#define _DECLARE_FABM_ARGS_0D_ type (type_environment),intent(inout) :: environment;_LOCATION_TYPE_,intent(in) :: LOCATION
#define _FABM_ARGS_IN_0D_ root%environment _ARG_LOCATION_

! Spatial loop for quantities defined on hortizontal slice of the full spatial domain.
#define _FABM_ENTER_HZ_
#define _FABM_LEAVE_HZ_

! Expressions for indexing space-dependent FABM variables defined on horizontal slices of the domain.
#define _INDEX_SURFACE_EXCHANGE_(index) (index)

#ifdef _FABM_USE_1D_LOOP_

! 1D vectorized: FABM subroutines operate on one spatial dimension.

! Dummy argument and argument declaration for location specification.
#define _LOCATION_ND_ fabm_loop_start,fabm_loop_stop _ARG_LOCATION_1DLOOP_
#define _DECLARE_LOCATION_ARG_ND_ _LOCATION_TYPE_,intent(in) :: fabm_loop_start,fabm_loop_stop _ARG_LOCATION_1DLOOP_; _LOCATION_TYPE_ :: _VARIABLE_1DLOOP_

! Beginning and end of spatial loop
#define _FABM_LOOP_BEGIN_ do _VARIABLE_1DLOOP_=fabm_loop_start,fabm_loop_stop
#define _FABM_LOOP_END_ end do

! Dimensionality of generic space-dependent arguments.
#define _ATTR_DIMENSIONS_0_ ,dimension(:)
#define _ATTR_DIMENSIONS_1_ ,dimension(:,:)
#define _ATTR_DIMENSIONS_2_ ,dimension(:,:,:)

! Expressions for indexing space-dependent FABM variables defined on the full spatial domain.
! These may be overridden by the host-specific driver (if it needs another order of dimensions).
! In that case, do not redefine the expressions here.
#ifndef _INDEX_ODE_
#define _INDEX_ODE_(variable) (_VARIABLE_1DLOOP_-fabm_loop_start+1,variable)
#endif
#ifndef _INDEX_PPDD_
#define _INDEX_PPDD_(variable1,variable2) (_VARIABLE_1DLOOP_-fabm_loop_start+1,variable1,variable2)
#endif
#ifndef _INDEX_CONSERVED_QUANTITY_
#define _INDEX_CONSERVED_QUANTITY_(variable) (_VARIABLE_1DLOOP_-fabm_loop_start+1,variable)
#endif
#ifndef _INDEX_VERTICAL_MOVEMENT_
#define _INDEX_VERTICAL_MOVEMENT_(variable) (_VARIABLE_1DLOOP_-fabm_loop_start+1,variable)
#endif
#define _INDEX_EXTINCTION_ (_VARIABLE_1DLOOP_-fabm_loop_start+1)

#else

! Not vectorized: FABM subroutines operate one the local state only.

! Dummy argument and argument declaration for location specification.
#define _LOCATION_ND_ LOCATION
#define _DECLARE_LOCATION_ARG_ND_ _LOCATION_TYPE_,intent(in) :: LOCATION

! Beginning and end of spatial loop
#define _FABM_LOOP_BEGIN_
#define _FABM_LOOP_END_

! Dimensionality of generic space-dependent arguments.
#define _ATTR_DIMENSIONS_0_
#define _ATTR_DIMENSIONS_1_ ,dimension(:)
#define _ATTR_DIMENSIONS_2_ ,dimension(:,:)

! Expressions for indexing space-dependent FABM variables defined on the full spatial domain.
#define _INDEX_ODE_(variable) (variable)
#define _INDEX_PPDD_(variable1,variable2) (variable1,variable2)
#define _INDEX_EXTINCTION_
#define _INDEX_CONSERVED_QUANTITY_(variable) (variable)
#define _INDEX_VERTICAL_MOVEMENT_(variable) (variable)

#endif

! For FABM: standard arguments used in calling biogeochemical routines.
#define _FABM_ARGS_ND_IN_ root%environment,_LOCATION_ND_

! For BGC models: FABM arguments to routines implemented by biogeochemical models.
#define _FABM_ARGS_ND_ environment,_LOCATION_ND_
#define _FABM_ARGS_DO_RHS_ _FABM_ARGS_ND_,rhs
#define _FABM_ARGS_DO_PPDD_ _FABM_ARGS_ND_,pp,dd
#define _FABM_ARGS_GET_EXTINCTION_ _FABM_ARGS_ND_,extinction
#define _FABM_ARGS_GET_CONSERVED_QUANTITIES_ _FABM_ARGS_ND_,sums
#define _FABM_ARGS_GET_SURFACE_EXCHANGE_ _FABM_ARGS_0D_,flux

! For BGC models: Declaration of FABM arguments to routines implemented by biogeochemical models.
#define _DECLARE_FABM_ARGS_ND_ type (type_environment),intent(inout) :: environment;_DECLARE_LOCATION_ARG_ND_
#define _DECLARE_FABM_ARGS_DO_RHS_  _DECLARE_FABM_ARGS_ND_;REALTYPE _ATTR_DIMENSIONS_1_,intent(inout) :: rhs
#define _DECLARE_FABM_ARGS_DO_PPDD_ _DECLARE_FABM_ARGS_ND_;REALTYPE _ATTR_DIMENSIONS_2_,intent(inout) :: pp,dd
#define _DECLARE_FABM_ARGS_GET_EXTINCTION_ _DECLARE_FABM_ARGS_ND_;REALTYPE _ATTR_DIMENSIONS_0_,intent(inout) :: extinction
#define _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_ _DECLARE_FABM_ARGS_ND_;REALTYPE _ATTR_DIMENSIONS_1_,intent(inout) :: sums
#define _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_ _DECLARE_FABM_ARGS_0D_;REALTYPE,dimension(:),intent(inout) :: flux

! Macros for declaring/accessing variable identifiers of arbitrary type.
#define _TYPE_STATE_VARIABLE_ID_ type (type_state_variable_id)
#define _TYPE_DIAGNOSTIC_VARIABLE_ID_ integer
#define _TYPE_DEPENDENCY_ID_ integer
#define _TYPE_CONSERVED_QUANTITY_ID_ integer

! For BGC models: Expressions for setting space-dependent FABM variables defined on the full spatial domain.
#define _SET_ODE_(variable,value) rhs _INDEX_ODE_(variable%id) = rhs _INDEX_ODE_(variable%id) + value
#define _SET_DD_(variable1,variable2,value) dd _INDEX_PPDD_(variable1%id,variable2%id) = dd _INDEX_PPDD_(variable1%id,variable2%id) + value
#define _SET_PP_(variable1,variable2,value) pp _INDEX_PPDD_(variable1%id,variable2%id) = pp _INDEX_PPDD_(variable1%id,variable2%id) + value
#define _SET_EXTINCTION_(value) extinction _INDEX_EXTINCTION_ = extinction _INDEX_EXTINCTION_ + value
#define _SET_CONSERVED_QUANTITY_(variable,value) sums _INDEX_CONSERVED_QUANTITY_(variable) = sums _INDEX_CONSERVED_QUANTITY_(variable) + value
#define _SET_VERTICAL_MOVEMENT_(variable,value) vertical_movement _INDEX_VERTICAL_MOVEMENT_(variable%id) = value
#define _SET_SURFACE_EXCHANGE_(variable,value) flux _INDEX_SURFACE_EXCHANGE_(variable%id) = value

! For BGC models: quick expressions for setting a single element in both the destruction and production matrix.
#define _SET_DD_SYM_(variable1,variable2,value) _SET_DD_(variable1,variable2,value);_SET_PP_(variable1,variable2,value)
#define _SET_PP_SYM_(variable1,variable2,value) _SET_PP_(variable1,variable2,value);_SET_DD_(variable1,variable2,value)

! For BGC models: read-only access to values of external dependencies
#define _GET_DEPENDENCY_(variable,target) target = environment%var(variable)%data _INDEX_LOCATION_
#define _GET_DEPENDENCY_HZ_(variable,target) target = environment%var_hz(variable)%data _INDEX_LOCATION_HZ_

! For FABM: read/write access to state variables
#define _GET_STATE_EX_(env,variable,target) target = env%var(variable%dependencyid)%data _INDEX_LOCATION_
#define _SET_STATE_EX_(env,variable,value) env%var(variable%dependencyid)%data _INDEX_LOCATION_ = value

! For BGC models: read-only access to state variable values
#define _GET_STATE_(variable,target) _GET_STATE_EX_(environment,variable,target)

! For BGC models: write access to diagnostic variables
#ifdef _FABM_MANAGE_DIAGNOSTICS_
#define _SET_DIAG_(index,value) if (index.ne.id_not_used) environment%diag(index,LOCATION) = value
#define _SET_DIAG_HZ_(index,value) if (index.ne.id_not_used) environment%diag_hz(index) = value
#else
#define _SET_DIAG_(index,value) if (index.ne.id_not_used) environment%var(index)%data _INDEX_LOCATION_ = value
#define _SET_DIAG_HZ_(index,value) if (index.ne.id_not_used) environment%var_hz(index)%data _INDEX_LOCATION_HZ_ = value
#endif

!#define _SET_ODE_(variable,value) if (variable%id.ne.id_not_used) rhs _INDEX_ODE_(variable%id) = rhs _INDEX_ODE_(variable%id) + value