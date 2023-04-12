add_library(fabm_models_imr OBJECT
    imr_model_library.F90
    imr_norwecom.F90
    )

add_dependencies(fabm_models_imr fabm_base)
