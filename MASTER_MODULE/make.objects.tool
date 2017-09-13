# ------------------------------------------------------------
# environment & interface modules external to the main package
# ------------------------------------------------------------
OBJECTS_GEOCAPE_MAIN  = $(OBJ)/geocape_tool.o

OBJECTS_GEOCAPE_PREP  = $(OBJ)/geocape_profile_prep.o               \
                           $(OBJ)/generate_plume.o                  \
                           $(OBJ)/geocape_surface_prep.o            \
                           $(OBJ)/solarspec_prep.o                  \
                           $(OBJ)/geocape_xsecs_prep.o              \
                           $(OBJ)/get_hitran_crs.o                  \
                           $(OBJ)/bspline.o                         \
                           $(OBJ)/prepare_aercld_optical_property.o \
                           $(OBJ)/index_water.o                     \
                           $(OBJ)/gauss.o                           \
                           $(OBJ)/julday.o                          \
                           $(OBJ)/caldat.o                          \
                           $(OBJ)/view_angles.o                     \
                           $(OBJ)/solar_angles.o                    \
                           $(OBJ)/satellite_angles.o            

OBJECTS_NETCDF        = $(OBJ)/netcdf_wrt.o

OBJECTS_MODULE        = $(GC_OBJ)/GC_parameters_module.o            \
                        $(GC_OBJ)/GC_control_types_module.o         \
                        $(GC_OBJ)/GC_variables_module.o             \
                        $(GC_OBJ)/GC_error_module.o                 \
                        $(GC_OBJ)/GC_read_input_module.o            \
                        $(GC_OBJ)/GC_aerosols_module.o              \
                        $(GC_OBJ)/GC_clouds_module.o                \
                        $(GC_OBJ)/GC_solar_module.o                 \
                        $(GC_OBJ)/GC_xsections_module.o             \
                        $(GC_OBJ)/GC_surface_module.o               \
                        $(GC_OBJ)/GC_profiles_module.o              \
                        $(GC_OBJ)/GC_Vlidort_module.o               \
                        $(GC_OBJ)/GC_convolution_module.o           \
                        $(GC_OBJ)/GC_netcdf_module.o