# ------------------------------------------------------------
# environment & interface modules external to the main package
# ------------------------------------------------------------
OBJECTS_GEOCAPE_V2P6  = $(OBJ)/geocape_tool_v2p6.o
OBJECTS_GEOCAPE_V2P6G = $(OBJ)/geocape_tool_v2p6G.o

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

OBJECTS_MODULEG       = $(GC_OBJG)/GC_parameters_module.o            \
                        $(GC_OBJG)/GC_variables_module.o             \
                        $(GC_OBJG)/GC_error_module.o                 \
                        $(GC_OBJG)/GC_read_input_module.o            \
                        $(GC_OBJG)/GC_aerosols_module.o              \
                        $(GC_OBJG)/GC_clouds_module.o                \
                        $(GC_OBJG)/GC_solar_module.o                 \
                        $(GC_OBJG)/GC_xsections_module.o             \
                        $(GC_OBJG)/GC_surface_module.o               \
                        $(GC_OBJG)/GC_profiles_module.o              \
                        $(GC_OBJG)/GC_Vlidort_module.o               \
                        $(GC_OBJG)/GC_convolution_module.o           \
                        $(GC_OBJG)/GC_netcdf_module.o
